"""
Pre-AF3 sequence splicing.

Given a full-length protein sequence and a dict of segment (exon/repeat)
coordinates on that sequence, build the spliced amino acid string for a
chosen set of segments, in the order given. This is the input-generation
counterpart to sequence_mapping.build_position_maps(), which does the
reverse operation on already-folded AF3/PDB output.

Ported from the splicing logic in sequence_edit.ipynb (get_spliced_ldlr,
get_spliced_vldlr), generalized to work for any receptor or Reelin
construct using the same segment_coords pattern as the rest of the
package.
"""

import re


def clean_sequence(sequence: str) -> str:
    """
    Strip whitespace/newlines from a pasted sequence (e.g. copied straight
    out of UniProt across multiple lines) so it's just amino acid letters.
    """
    return re.sub(r"\s+", "", sequence)


def splice_sequence(full_sequence: str, segment_coords: dict, selected_segments: list) -> str:
    """
    Build a spliced amino acid sequence by concatenating selected segments,
    in order, from the full-length sequence.

    Parameters
    ----------
    full_sequence : str
        The full-length protein sequence (e.g. from UniProt). Whitespace is
        stripped automatically, but it should otherwise contain only amino
        acid letters -- no headers, FASTA ">" lines, or other characters.
    segment_coords : dict
        {segment_id: (full_seq_start, full_seq_end)}, 1-indexed inclusive.
        This is the same dict shape used elsewhere in the package (e.g.
        RECEPTOR_CONFIGS[...]["exon_coords"] or Reelin repeat coordinates).
    selected_segments : list
        Segment IDs to splice together, IN ORDER. Order determines the
        resulting sequence's order -- this is not auto-sorted, since you
        may deliberately want a non-sequential splice.

    Returns
    -------
    str
        The spliced amino acid sequence.

    Raises
    ------
    KeyError
        If a selected segment isn't a key in segment_coords.
    ValueError
        If a segment's coordinates fall outside the bounds of full_sequence
        (this usually means segment_coords and full_sequence are mismatched
        -- e.g. coordinates from a different UniProt version/isoform than
        the sequence string you're using).
    """
    full_sequence = clean_sequence(full_sequence)

    pieces = []
    for seg in selected_segments:
        if seg not in segment_coords:
            available = ", ".join(str(s) for s in sorted(segment_coords, key=str))
            raise KeyError(
                f"Segment {seg} not found in segment_coords. Available: {available}"
            )
        start, end = segment_coords[seg]
        if start < 1 or end > len(full_sequence):
            raise ValueError(
                f"Segment {seg} coordinates ({start}-{end}) fall outside the full "
                f"sequence length ({len(full_sequence)}). segment_coords and "
                f"full_sequence may be mismatched (e.g. different UniProt versions)."
            )
        pieces.append(full_sequence[start - 1:end])

    return "".join(pieces)


def segment_breakdown(full_sequence: str, segment_coords: dict) -> dict:
    """
    Return {segment_id: sequence_substring} for every segment defined in
    segment_coords. Useful as a sanity check before splicing -- e.g. to
    print each exon's sequence and length and confirm coordinates look
    right, the way the original notebooks printed exon-by-exon breakdowns.
    """
    full_sequence = clean_sequence(full_sequence)
    breakdown = {}
    for seg, (start, end) in segment_coords.items():
        breakdown[seg] = full_sequence[start - 1:end]
    return breakdown
