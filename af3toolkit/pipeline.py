"""
Top-level pipeline: analyze_complex() ties together contact detection,
position re-mapping, and domain/validation annotation for both chains
of an AF3 receptor-Reelin complex.
"""

from dataclasses import dataclass

import pandas as pd

from .receptors import get_receptor_config
from .reelin import get_reelin_repeat_coords, is_validated, REELIN_REPEATS
from .sequence_mapping import build_position_maps
from .contacts import find_contacts, format_contacts_as_dataframe


@dataclass
class AnalysisResult:
    contacts_table: pd.DataFrame      # full per-residue contact table, all columns
    domain_summary: pd.DataFrame      # contact counts per receptor domain
    excluded_domains: list            # receptor domains not present in this construct
    calcium_contacts: list            # (residue_name, full_seq_pos, exon, domain), if requested


def _annotate_receptor_side(df, residue_col, pdb_pos_col, receptor_config, selected_exons):
    """Add Full-Seq Pos / Exon / Domain columns for the receptor chain."""
    pos_maps = build_position_maps(receptor_config["exon_coords"], selected_exons)
    exon_to_domain = receptor_config["exon_to_domain"]

    def map_row(pdb_pos):
        if pdb_pos == "" or pd.isna(pdb_pos):
            return pd.Series([None, None, None])
        pdb_pos = int(pdb_pos)
        full_pos = pos_maps["to_full"].get(pdb_pos)
        exon = pos_maps["to_segment"].get(pdb_pos)
        domain = exon_to_domain.get(exon) if exon is not None else None
        return pd.Series([full_pos, exon, domain])

    df[["Receptor Full-Seq Pos", "Receptor Exon", "Receptor Domain"]] = df[pdb_pos_col].apply(map_row)
    return df, pos_maps


def _annotate_reelin_side(df, pdb_pos_col, selected_repeats):
    """Add Full-Seq Pos / Repeat / Validated columns for the Reelin chain."""
    repeat_coords = get_reelin_repeat_coords()
    missing = [r for r in selected_repeats if r not in repeat_coords]
    if missing:
        raise ValueError(
            f"Reelin repeat(s) {missing} have no full_seq_range set in reelin.py. "
            "Fill in REELIN_REPEATS before running this construct."
        )

    pos_maps = build_position_maps(repeat_coords, selected_repeats)

    def map_row(pdb_pos):
        if pdb_pos == "" or pd.isna(pdb_pos):
            return pd.Series([None, None, None])
        pdb_pos = int(pdb_pos)
        full_pos = pos_maps["to_full"].get(pdb_pos)
        repeat = pos_maps["to_segment"].get(pdb_pos)
        validated = is_validated(repeat) if repeat is not None else None
        return pd.Series([full_pos, repeat, validated])

    df[["Reelin Full-Seq Pos", "Reelin Repeat", "Reelin Validated"]] = df[pdb_pos_col].apply(map_row)
    return df


def _build_domain_summary(df, receptor_config, selected_exons):
    """Count contacts per receptor domain, and list domains excluded from this construct."""
    counts = (
        df[df["Receptor Domain"].notna()]
        .groupby("Receptor Domain")
        .size()
        .reset_index(name="Contact Count")
        .sort_values("Contact Count", ascending=False)
    )

    all_domains = set(receptor_config["exon_to_domain"].values())
    modeled_domains = {
        receptor_config["exon_to_domain"][e]
        for e in selected_exons
        if e in receptor_config["exon_to_domain"]
    }
    excluded_domains = sorted(all_domains - modeled_domains)

    return counts, excluded_domains


def analyze_complex(
    pdb_path: str,
    receptor: str,
    receptor_selected_exons: list,
    reelin_selected_repeats: list,
    receptor_chain: str = "A",
    reelin_chain: str = "B",
    distance_cutoff: float = 3.5,
    include_calcium: bool = True,
    calcium_cutoff: float = 3.5,
) -> AnalysisResult:
    """
    Run the full post-AF3 contact analysis for one receptor-Reelin complex.

    Parameters
    ----------
    pdb_path : str
        Path to the AF3-output PDB file for this construct.
    receptor : str
        One of "ApoER2", "LDLR", "VLDLR" (see receptors.RECEPTOR_CONFIGS).
    receptor_selected_exons : list[int]
        Exon numbers spliced into this construct, IN ORDER.
    reelin_selected_repeats : list[str]
        Reelin repeat PAIR labels spliced into this construct, IN ORDER
        (e.g. ["5-6"]). Pair labels match the keys in reelin.REELIN_REPEATS
        (currently resolved at pair-level, not individual repeat number).
    receptor_chain, reelin_chain : str
        Chain IDs in the PDB file.
    distance_cutoff : float
        Angstrom cutoff for residue-residue contacts.
    include_calcium : bool
        Whether to also run calcium-coordination analysis on the receptor chain.
    calcium_cutoff : float
        Angstrom cutoff for ion-coordination contacts.

    Returns
    -------
    AnalysisResult
    """
    receptor_config = get_receptor_config(receptor)

    chain1_contacts, chain2_contacts = find_contacts(
        pdb_path, chain1_id=receptor_chain, chain2_id=reelin_chain, distance_cutoff=distance_cutoff
    )
    df = format_contacts_as_dataframe(
        chain1_contacts, chain2_contacts, chain1_id=receptor_chain, chain2_id=reelin_chain
    )
    df = df.rename(columns={
        f"Chain {receptor_chain} Residue": "Receptor Residue",
        f"Chain {receptor_chain} PDB Pos": "Receptor PDB Pos",
        f"Chain {reelin_chain} Residue": "Reelin Residue",
        f"Chain {reelin_chain} PDB Pos": "Reelin PDB Pos",
    })

    df, _ = _annotate_receptor_side(
        df, "Receptor Residue", "Receptor PDB Pos", receptor_config, receptor_selected_exons
    )
    df = _annotate_reelin_side(df, "Reelin PDB Pos", reelin_selected_repeats)

    column_order = [
        "Receptor Residue", "Receptor PDB Pos", "Receptor Full-Seq Pos", "Receptor Exon", "Receptor Domain",
        "Reelin Residue", "Reelin PDB Pos", "Reelin Full-Seq Pos", "Reelin Repeat", "Reelin Validated",
    ]
    df = df[column_order]

    domain_summary, excluded_domains = _build_domain_summary(df, receptor_config, receptor_selected_exons)

    calcium_contacts = []
    if include_calcium:
        from .contacts import find_chain_ion_contacts
        ca_raw = find_chain_ion_contacts(
            pdb_path, chain_id=receptor_chain, ion_resname="CA", distance_cutoff=calcium_cutoff
        )
        pos_maps = build_position_maps(receptor_config["exon_coords"], receptor_selected_exons)
        for resname, pdb_pos in ca_raw:
            full_pos = pos_maps["to_full"].get(pdb_pos)
            exon = pos_maps["to_segment"].get(pdb_pos)
            domain = receptor_config["exon_to_domain"].get(exon) if exon is not None else None
            calcium_contacts.append((resname, full_pos, exon, domain))

    return AnalysisResult(
        contacts_table=df,
        domain_summary=domain_summary,
        excluded_domains=excluded_domains,
        calcium_contacts=calcium_contacts,
    )
