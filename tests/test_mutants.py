from src.generate_mutants import double_mutant_rows, single_mutant_rows


class DummyPair:
    pass


def test_single_mutant_count():
    import pandas as pd

    pair = pd.Series({"pair_id": "x", "protein_name": "p", "protein_sequence": "AAAA", "target_dna": "ACGT"})
    rows = single_mutant_rows(pair)
    assert len(rows) == 12


def test_double_mutant_count():
    import pandas as pd

    pair = pd.Series({"pair_id": "x", "protein_name": "p", "protein_sequence": "AAAA", "target_dna": "ACGT"})
    rows = double_mutant_rows(pair)
    assert len(rows) == 54

