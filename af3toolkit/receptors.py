"""
Receptor configuration registry.

Each receptor entry holds:
  - exon_coords: {exon_number: (full_sequence_start, full_sequence_end)}
                 1-indexed, inclusive, based on the full-length UniProt sequence.
  - exon_to_domain: {exon_number: domain_name} -- structural/functional domain
                 each exon falls into. LDL-receptor-family proteins are
                 modular (LA repeats, EGF-like / beta-propeller domains,
                 O-linked sugar domain, transmembrane, cytoplasmic tail),
                 and exon boundaries roughly track domain boundaries.

IMPORTANT: exon_coords below were transcribed from UniProt sequences used
in the original analysis notebooks and should be re-verified against the
current UniProt entry if you update the sequence. exon_to_domain labels
are based on general LDL-receptor-family domain architecture and are
APPROXIMATE -- verify against UniProt's "Family & Domains" feature table
for each receptor before citing exact boundaries in a paper or talk.
"""

RECEPTOR_CONFIGS = {
    "ApoER2": {
        "uniprot_id": "Q14114",  # LRP8 / ApoER2 -- verify
        "exon_coords": {
            1: (1, 42),
            2: (43, 82),
            3: (83, 123),
            4: (124, 166),
            5: (167, 295),
            6: (296, 336),
            7: (337, 376),
            8: (377, 418),
            9: (419, 476),
            10: (477, 552),
            11: (553, 592),
            12: (593, 638),
            13: (639, 686),
            14: (687, 737),
            15: (738, 812),
            16: (813, 835),
            17: (836, 892),
            18: (893, 951),
            19: (952, 963),
        },
        # verify each label against UniProt Q14114 feature table.
        "exon_to_domain": {
            1: "Signal peptide",
            2: "LA repeat 1",
            3: "LA repeat 2",
            4: "LA repeat 3",
            5: "LA repeat 4",
            6: "LA repeat 5",
            7: "LA repeat 6",
            8: "LA repeat 7",
            9: "EGF-precursor homology domain A",
            10: "EGF-precursor homology domain A",
            11: "YWTD beta-propeller",
            12: "YWTD beta-propeller",
            13: "EGF-precursor homology domain B",
            14: "EGF-precursor homology domain B",
            15: "O-linked sugar domain / proline-rich insert",
            16: "O-linked sugar domain / proline-rich insert",
            17: "Transmembrane domain",
            18: "Cytoplasmic tail (NPxY motif)",
            19: "Cytoplasmic tail (NPxY motif)",
        },
    },
    "LDLR": {
        "uniprot_id": "P01130",  # verify
        "exon_coords": {
            1: (1, 23),
            2: (23, 64),
            3: (64, 105),
            4: (105, 232),
            5: (232, 273),
            6: (273, 314),
            7: (314, 354),
            8: (354, 396),
            9: (396, 453),
            10: (453, 529),
            11: (529, 569),
            12: (569, 615),
            13: (616, 663),
            14: (663, 714),
            15: (714, 771),
            16: (771, 797),
            17: (797, 849),
            18: (850, 860),
        },
        # verify each label against UniProt P01130 feature table.
        "exon_to_domain": {
            1: "Signal peptide",
            2: "LA repeat 1",
            3: "LA repeat 2",
            4: "LA repeat 3",
            5: "LA repeat 4",
            6: "LA repeat 5",
            7: "LA repeat 6",
            8: "LA repeat 7",
            9: "EGF-precursor homology domain A",
            10: "YWTD beta-propeller",
            11: "YWTD beta-propeller",
            12: "YWTD beta-propeller",
            13: "EGF-precursor homology domain B",
            14: "EGF-precursor homology domain B",
            15: "EGF-precursor homology domain B",
            16: "O-linked sugar domain",
            17: "Transmembrane domain",
            18: "Cytoplasmic tail (NPxY motif)",
        },
    },
    "VLDLR": {
        "uniprot_id": "P98155",  # verify
        "exon_coords": {
            1: (1, 28),
            2: (29, 68),
            3: (69, 109),
            4: (110, 150),
            5: (151, 274),
            6: (275, 315),
            7: (316, 356),
            8: (357, 396),
            9: (397, 438),
            10: (439, 495),
            11: (496, 568),
            12: (569, 608),
            13: (609, 654),
            14: (655, 702),
            15: (703, 751),
            16: (752, 779),
            17: (780, 806),
            18: (807, 862),
            19: (863, 873),
        },
        # verify each label against UniProt P98155 feature table.
        "exon_to_domain": {
            1: "Signal peptide",
            2: "LA repeat 1",
            3: "LA repeat 2",
            4: "LA repeat 3",
            5: "LA repeat 4",
            6: "LA repeat 5",
            7: "LA repeat 6",
            8: "LA repeat 7",
            9: "LA repeat 8",
            10: "EGF-precursor homology domain A",
            11: "YWTD beta-propeller",
            12: "YWTD beta-propeller",
            13: "EGF-precursor homology domain B",
            14: "EGF-precursor homology domain B",
            15: "O-linked sugar domain",
            16: "O-linked sugar domain",
            17: "Transmembrane domain",
            18: "Cytoplasmic tail (NPxY motif)",
            19: "Cytoplasmic tail (NPxY motif)",
        },
    },
}


def get_receptor_config(receptor_name: str) -> dict:
    """Look up a receptor config by name, with a helpful error if not found."""
    try:
        return RECEPTOR_CONFIGS[receptor_name]
    except KeyError:
        available = ", ".join(sorted(RECEPTOR_CONFIGS))
        raise KeyError(
            f"Receptor '{receptor_name}' not found. Available receptors: {available}"
        )


def list_receptors() -> list:
    """Return the list of receptor names currently in the registry."""
    return sorted(RECEPTOR_CONFIGS)
