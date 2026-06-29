"""
Tests for the core position-mapping logic. This is the piece most worth
protecting with tests, since a silent off-by-one here would mislabel
every downstream domain/repeat call without throwing an error.
"""

from af3toolkit.sequence_mapping import build_position_maps


def test_single_segment_starts_at_one():
    coords = {1: (1, 42), 2: (43, 82)}
    maps = build_position_maps(coords, selected_segments=[2])
    # First residue of the spliced construct (local pos 1) should map
    # to the first residue of exon 2 in the full sequence (43), not 1.
    assert maps["to_full"][1] == 43
    assert maps["to_segment"][1] == 2


def test_multi_segment_continuity():
    coords = {2: (43, 82), 3: (83, 123)}
    maps = build_position_maps(coords, selected_segments=[2, 3])
    # Exon 2 spans 40 residues (43-82 inclusive), so local position 41
    # should be the first residue of exon 3 (full-seq pos 83).
    assert maps["to_full"][40] == 82
    assert maps["to_full"][41] == 83
    assert maps["to_segment"][40] == 2
    assert maps["to_segment"][41] == 3


def test_missing_segment_raises():
    coords = {1: (1, 42)}
    try:
        build_position_maps(coords, selected_segments=[2])
        assert False, "Expected KeyError for missing segment"
    except KeyError:
        pass


def test_excluded_first_exon_shifts_numbering():
    # This is the exact scenario from the conversation: exon 1 (1-42) is
    # excluded, exon 2 (43-82) is first in the construct. PDB position 1
    # should NOT be full-seq position 1.
    coords = {1: (1, 42), 2: (43, 82)}
    maps = build_position_maps(coords, selected_segments=[2])
    assert maps["to_full"][1] != 1
    assert maps["to_full"][1] == 43
