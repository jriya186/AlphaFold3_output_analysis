# af3toolkit

Post-AlphaFold3 structural contact analysis for LDL-receptor-family / Reelin complexes (ApoER2, LDLR, VLDLR).

## Background

This toolkit was developed based on research conducted during an internship in the Ho-Beffert Lab at Boston University, under the supervision of Dr. Uwe Beffert and Dr. Angela Ho. The original exploratory analysis was conducted in the form of Jupyter notebooks during the internship. This package is a subsequent independent repackaging of that work into a structured, reusable Python toolkit.

The goal of the study was to analyze the binding interface between ApoER2 (LRP8) splice variants and the extracellular glycoprotein Reelin. ApoER2 is a member of the LDL receptor superfamily and a key receptor in the Reelin signaling pathway, which governs neuronal migration during brain development and synaptic plasticity in the adult brain. Reelin signals through ApoER2 (and VLDLR) by binding to the receptor's extracellular LA repeat domains, triggering phosphorylation of the adaptor protein DAB1 and downstream cytoskeletal reorganization. Disruption of this pathway has been implicated in neurodevelopmental disorders and Alzheimer's disease, where ApoE, another ApoER2 ligand, competes for overlapping binding sites.

ApoER2 undergoes extensive alternative splicing, producing isoforms with different combinations of extracellular exons. Understanding which splice variant-specific domains mediate Reelin binding, and with what geometry, has direct implications for understanding isoform-specific receptor function in neurodevelopment and disease. AlphaFold3 was used to model binding between specific ApoER2 exon combinations and Reelin repeat fragments; this toolkit provides the post-AF3 analysis pipeline for interpreting those structural predictions.

## What this does

AF3 inputs for this kind of analysis are spliced fragments, not full-length proteins (e.g. only exons 2-8 of a receptor, or only repeats 5-6 of Reelin). That means raw PDB residue numbers do NOT correspond to real positions in the full-length protein. This toolkit re-maps every contact back to:

- the correct full-sequence position
- the exon (receptor) or repeat (Reelin) it falls in
- the structural domain it belongs to (receptor side)
- whether that region has experimental structural support or is AF3 prediction only (Reelin side)

and explicitly flags which domains were not part of the modeled construct, so "no contacts in domain X" is never confused with "no binding in domain X."

Note: only Reelin repeats 5-6 have a published crystal structure of the receptor-bound interface. Contacts in all other repeats are AF3 predictions with no experimental structure to validate them against, and are flagged accordingly in all output.

## Install

```
pip install -e .

# optional, for nicer CLI tables:
pip install -e ".[cli-pretty]"
```

## CLI usage

```
af3toolkit analyze \
  --pdb apoer2_reelin56_complex.pdb \
  --receptor ApoER2 \
  --receptor-exons 2,3,4,5,6,7,8 \
  --reelin-repeats 5-6 \
  --save-csv results/apoer2_repeat56_contacts.csv
```

`--save-csv` is optional -- omit it and the full table prints to terminal without writing a file.

Discover what is already configured:

```
af3toolkit list-receptors
af3toolkit list-reelin-repeats
```

## Output

Every contact row includes both the receptor and Reelin side, fully re-mapped:

| Receptor Residue | Receptor PDB Pos | Receptor Full-Seq Pos | Receptor Exon | Receptor Domain | Reelin Residue | Reelin PDB Pos | Reelin Full-Seq Pos | Reelin Repeat | Reelin Validated |
|---|---|---|---|---|---|---|---|---|---|

Above the table, the CLI also prints:

- a domain-level contact summary (counts per receptor domain)
- a list of domains excluded from the construct (not modeled, not evidence of no binding)
- calcium-coordinating residues, since LA repeats are calcium-dependent

## Pre-AF3 sequence splicing (Python only, no CLI)

Before running AF3, you often need a spliced construct -- specific exons or repeats concatenated together, not the full-length protein. `sequence_splicing.py` handles this, reusing the same exon/repeat coordinate dicts as the rest of the package.

```python
from af3toolkit.sequence_splicing import splice_sequence, clean_sequence
from af3toolkit.receptors import RECEPTOR_CONFIGS

full_seq = clean_sequence("""
MGPWGWKLRWTVALLLAAAGTAV...  # paste full UniProt sequence, newlines are fine
""")

spliced = splice_sequence(
    full_seq,
    RECEPTOR_CONFIGS["LDLR"]["exon_coords"],
    selected_segments=[2, 3],
)
print(spliced)  # ready to paste into AF3
```

For Reelin, swap in its repeat coordinates:

```python
from af3toolkit.reelin import get_reelin_repeat_coords

reelin_spliced = splice_sequence(
    reelin_full_seq,
    get_reelin_repeat_coords(),
    selected_segments=["5-6"],
)
```

`segment_breakdown()` is also available if you want to print every exon's sequence before splicing, as a sanity check.

## Why receptor and Reelin sides are treated differently

Only Reelin repeats 5-6 have a published crystal structure of the actual receptor-bound interface. Every other repeat's contacts in an AF3 model are predictions with no experimental structure to check them against. The `Reelin Validated` column makes that distinction explicit in the output, rather than implying uniform confidence across the whole protein.

## Before using this on your own data

- `af3toolkit/receptors.py`: exon-to-domain labels are approximate (based on general LDL-receptor-family architecture) -- verify against each receptor's UniProt feature table before citing exact domain boundaries.
- `af3toolkit/reelin.py`: repeat coordinates are set for repeat pairs 1-2, 3-4, 5-6, and 7-8 based on human Reelin (UniProt P78509) Ensembl coordinates. Only repeats 5-6 are experimentally validated; all others are AF3-prediction only.

## Repo structure

```
af3toolkit/
  receptors.py         # exon coordinates + domain labels (ApoER2, LDLR, VLDLR)
  reelin.py            # repeat coordinates + validated/prediction-only tags
  sequence_mapping.py  # generic local-PDB-position -> full-sequence mapper
  sequence_splicing.py # pre-AF3 input generation: splice exons/repeats into a sequence string
  contacts.py          # Biopython-based contact + ion-coordination detection
  confidence.py        # best-AF3-model selection by iPTM/pTM
  pipeline.py          # ties it all together: analyze_complex()
  cli.py               # command-line interface
notebooks/archive/     # original exploratory notebooks this package was built from
                       # (kept for reference, not maintained)
                       # sequences.yml documents the conda environment they were run in
tests/                 # pytest tests for the position-mapping logic
```

## Development

```
pip install -e ".[dev]"
pytest tests/
```
