"""
Generic local-PDB-position -> full-sequence-position mapping.

Used identically for receptor exons and Reelin repeats: given a dict of
segment_id -> (full_seq_start, full_seq_end) and the ordered list of
segments actually spliced into a given AF3 construct, reconstruct what
each raw PDB residue number actually corresponds to in the full-length
protein.

This deliberately does NOT assume PDB numbering matches full-sequence
numbering -- AF3 constructs almost never start at the true residue 1
unless the very first exon/repeat happens to be included.
"""


def build_position_maps(segment_coords: dict, selected_segments: list) -> dict:
    """
    Parameters
    ----------
    segment_coords : dict
        {segment_id: (full_seq_start, full_seq_end)}, 1-indexed inclusive,
        covering the full-length protein (e.g. RECEPTOR_CONFIGS[...]["exon_coords"]
        or a dict built from REELIN_REPEATS).
    selected_segments : list
        Segment IDs, IN THE ORDER they were spliced together to build the
        AF3 input construct. Order matters -- it determines local numbering.

    Returns
    -------
    dict with two lookups:
        "to_full":    {local_pdb_position: full_sequence_position}
        "to_segment": {local_pdb_position: segment_id}

    Raises
    ------
    KeyError if a selected segment isn't in segment_coords.
    """
    pdb_to_full = {}
    pdb_to_segment = {}
    pdb_pos = 1

    for seg in selected_segments:
        if seg not in segment_coords:
            available = ", ".join(str(s) for s in sorted(segment_coords))
            raise KeyError(
                f"Segment {seg} not found in segment_coords. Available: {available}"
            )
        start, end = segment_coords[seg]
        for full_pos in range(start, end + 1):
            pdb_to_full[pdb_pos] = full_pos
            pdb_to_segment[pdb_pos] = seg
            pdb_pos += 1

    return {"to_full": pdb_to_full, "to_segment": pdb_to_segment}
