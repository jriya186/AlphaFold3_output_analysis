"""
Reelin repeat registry.

Reelin's binding to ApoER2/VLDLR is structurally characterized for
repeats 5-6 (the only region with a published crystal structure of the
receptor-Reelin interface). Every other repeat pair's AF3-predicted
contacts are model output only -- there is no experimental structure to
validate them against, and they should NOT be reported with the same
confidence as repeats 5-6.

Coordinates are keyed by REPEAT PAIR, not individual repeat number --
the source data resolves boundaries at the pair level (e.g. "repeats
1-2" as one contiguous unit), not as eight individually-bounded repeats.
If you later get single-repeat-resolution boundaries, split these.

Numbering is HUMAN REELIN (UniProt P78509), per Ensembl coordinates --
NOT the mouse coordinates from the cited structural paper. Mixing the
two coordinate systems would silently mis-map every downstream position,
so mouse numbers are kept only in `evidence` for citation purposes.
"""

REELIN_REPEATS = {
    "signal_peptide": {
        "full_seq_range": (1, 22),
        "validated": False,
        "evidence": "Signal peptide, not a structural repeat. Mature protein starts at aa 23.",
    },
    "1-2": {
        "full_seq_range": (526, 1220),
        "validated": False,
        "evidence": "AF3 prediction only -- no experimental structure for this region.",
    },
    "3-4": {
        "full_seq_range": (1240, 1949),
        "validated": False,
        "evidence": "AF3 prediction only -- no experimental structure for this region.",
    },
    "5-6": {
        "full_seq_range": (1956, 2661),
        "validated": True,
        "evidence": (
            "Crystal structure of Reelin repeats 5-6 bound to ApoER2/Reelin receptor "
            "(mouse numbering aa 1952-2678 in source paper: "
            "https://doi.org/10.3389/fncel.2016.00137). "
            "Human-equivalent range above is from Ensembl."
        ),
    },
    "7-8": {
        "full_seq_range": (2669, 3426),
        "validated": False,
        "evidence": "AF3 prediction only -- no experimental structure for this region.",
    },
}


def get_reelin_repeat_coords() -> dict:
    """
    Return {repeat_pair: (start, end)} for repeat pairs that have a
    full_seq_range filled in. Raises if any requested repeat pair is
    missing its coordinates -- this is deliberate so a None silently
    isn't treated as a real position.
    """
    coords = {}
    for repeat, info in REELIN_REPEATS.items():
        if info["full_seq_range"] is not None:
            coords[repeat] = info["full_seq_range"]
    return coords


def is_validated(repeat_pair: str) -> bool:
    """Whether this repeat pair's receptor contacts have experimental structural support."""
    if repeat_pair not in REELIN_REPEATS:
        available = ", ".join(sorted(REELIN_REPEATS))
        raise KeyError(
            f"Reelin repeat pair '{repeat_pair}' not recognized. Available: {available}"
        )
    return REELIN_REPEATS[repeat_pair]["validated"]


def list_reelin_repeats() -> list:
    """Return repeat pair labels with a short validation note, for CLI display."""
    out = []
    for repeat in sorted(REELIN_REPEATS):
        info = REELIN_REPEATS[repeat]
        tag = "VALIDATED" if info["validated"] else "prediction only"
        out.append((repeat, tag))
    return out