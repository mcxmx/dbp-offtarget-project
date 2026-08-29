from src.utils import gc_content, hamming_distance, reverse_complement, reverse_complement_canonical, sequence_identity


def test_reverse_complement():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAGC") == "GCTT"


def test_gc_content():
    assert gc_content("AAAA") == 0.0
    assert gc_content("GGCC") == 1.0


def test_hamming_distance_and_identity():
    assert hamming_distance("AAAA", "AAAT") == 1
    assert sequence_identity("AAAA", "AAAT") == 0.75


def test_reverse_complement_canonical():
    assert reverse_complement_canonical("ACGT") == "ACGT"
    assert reverse_complement_canonical("ATGC") == "ATGC"
    assert reverse_complement_canonical("AAGC") == "AAGC"
