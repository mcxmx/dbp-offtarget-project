import pytest

from src.sequence_equivalence import canonical_rc, has_rc_split_leakage, rc_equivalent, reverse_complement


def test_reverse_complement_and_canonical_rc_are_deterministic():
    assert reverse_complement("AAGTCGA") == "TCGACTT"
    assert canonical_rc("AAGTCGA") == canonical_rc("TCGACTT")
    assert canonical_rc("tcgactt") == canonical_rc("AAGTCGA")


def test_rc_equivalent_pairs_share_canonical_class():
    seq = "CTGACG"
    rc = reverse_complement(seq)
    assert rc_equivalent(seq, rc)
    assert canonical_rc(seq) == canonical_rc(rc)


def test_invalid_dna_rejected():
    with pytest.raises(ValueError):
        reverse_complement("ACNGT")


def test_split_helper_detects_reverse_complement_leakage():
    train = {"AAAAAAA", "CTGACGA"}
    test_with_leak = {"TTTTTTT", "CCCCCCC"}
    test_without_leak = {"CCCCCCC", "GGGGAAT"}
    assert has_rc_split_leakage(train, test_with_leak)
    assert not has_rc_split_leakage(train, test_without_leak)
