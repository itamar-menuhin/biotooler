"""Test to probe pyChimera API contract for return shapes and indexing.

This test establishes a test-backed API contract for pyChimera's calc_cARS function
before implementing positional windowing. We verify:
1. What calc_cARS returns in codon/nt/aa modes
2. Whether return values are scalars or arrays
3. Array lengths and their relationship to input sequence length
4. Behavior with and without win_params

These tests skip cleanly when pychimera is not installed.
"""

import pytest

# Skip entire module if pychimera is not installed
pychimera = pytest.importorskip("pychimera")


def test_calc_cars_codon_mode_returns_scalar_without_win_params():
    """Test calc_cARS in codon mode without win_params returns a scalar.

    Expected: calc_cARS returns a single float score when return_vec is not specified.
    The score represents the average maximal common substring length.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS, nt2codon

    # Create tiny reference sequences (9 nt each = 3 codons)
    ref_seqs_nt = [
        "ATGAAATAA",  # ATG AAA TAA (M K *)
        "ATGGCATAA",  # ATG GCA TAA (M A *)
    ]

    # Create target sequence (9 nt = 3 codons)
    target_nt = "ATGAAATAA"  # Same as first reference

    # Convert to codon representation
    ref_codons = nt2codon(ref_seqs_nt)
    target_codons = nt2codon([target_nt])

    # Build suffix array from reference
    suffix_array = build_suffix_array(ref_codons)

    # Call calc_cARS without win_params (global mode)
    result = calc_cARS(target_codons[0], suffix_array, max_len=40, max_pos=0.5)

    # Assert it returns a scalar (float)
    assert isinstance(result, (float, int)), f"Expected scalar, got {type(result)}"
    assert result >= 0, "cARS score should be non-negative"


def test_calc_cars_codon_mode_returns_vector_with_return_vec():
    """Test calc_cARS in codon mode with return_vec=True returns an array of codon length.

    Expected: With return_vec=True, calc_cARS returns an array of per-codon maximal
    common substring lengths. Array length should match the number of codons.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS, nt2codon

    # Create tiny reference sequences (9 nt each = 3 codons)
    ref_seqs_nt = [
        "ATGAAATAA",  # ATG AAA TAA (M K *)
        "ATGGCATAA",  # ATG GCA TAA (M A *)
    ]

    # Create target sequence (9 nt = 3 codons)
    target_nt = "ATGAAATAA"

    # Convert to codon representation
    ref_codons = nt2codon(ref_seqs_nt)
    target_codons = nt2codon([target_nt])

    # Build suffix array from reference
    suffix_array = build_suffix_array(ref_codons)

    # Call calc_cARS with return_vec=True
    result = calc_cARS(
        target_codons[0], suffix_array, max_len=40, max_pos=0.5, return_vec=True
    )

    # Assert it returns an array-like object
    assert hasattr(result, "__len__"), "Expected array-like object with return_vec=True"
    assert hasattr(result, "__getitem__"), "Expected indexable array-like object"

    # Assert length matches codon count (9 nt / 3 = 3 codons)
    expected_codon_count = len(target_nt) // 3
    assert len(result) == expected_codon_count, (
        f"Expected vector length {expected_codon_count} (codon count), "
        f"got {len(result)}"
    )


def test_calc_cars_codon_mode_with_win_params_returns_scalar():
    """Test calc_cARS in codon mode with win_params returns a scalar (position-specific score).

    Expected: Even with win_params (position-specific mode), calc_cARS returns a single
    scalar score when return_vec is not specified. The win_params affect how matches
    are weighted by position but the output is still aggregated to a scalar.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS, nt2codon

    # Create tiny reference sequences
    ref_seqs_nt = [
        "ATGAAATAA",
        "ATGGCATAA",
    ]

    # Create target sequence
    target_nt = "ATGAAATAA"

    # Convert to codon representation
    ref_codons = nt2codon(ref_seqs_nt)
    target_codons = nt2codon([target_nt])

    # Build suffix array
    suffix_array = build_suffix_array(ref_codons)

    # Define window parameters for position-specific mode
    win_params = {"size": 2, "center": 0, "by_start": True, "by_stop": True}

    # Call calc_cARS with win_params (position-specific mode)
    result = calc_cARS(
        target_codons[0],
        suffix_array,
        win_params=win_params,
        max_len=40,
        max_pos=0.5,
    )

    # Assert it returns a scalar
    assert isinstance(result, (float, int)), f"Expected scalar, got {type(result)}"
    assert result >= 0, "PScARS score should be non-negative"


def test_calc_cars_codon_mode_with_win_params_and_return_vec():
    """Test calc_cARS in codon mode with win_params and return_vec=True returns array.

    Expected: With both win_params and return_vec=True, calc_cARS returns an array
    of per-codon scores. Array length matches codon count.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS, nt2codon

    # Create tiny reference sequences
    ref_seqs_nt = [
        "ATGAAATAA",
        "ATGGCATAA",
    ]

    # Create target sequence
    target_nt = "ATGAAATAA"

    # Convert to codon representation
    ref_codons = nt2codon(ref_seqs_nt)
    target_codons = nt2codon([target_nt])

    # Build suffix array
    suffix_array = build_suffix_array(ref_codons)

    # Define window parameters
    win_params = {"size": 2, "center": 0, "by_start": True, "by_stop": True}

    # Call calc_cARS with win_params and return_vec=True
    result = calc_cARS(
        target_codons[0],
        suffix_array,
        win_params=win_params,
        max_len=40,
        max_pos=0.5,
        return_vec=True,
    )

    # Assert it returns an array-like object
    assert hasattr(result, "__len__"), "Expected array-like object"
    assert hasattr(result, "__getitem__"), "Expected indexable array-like object"

    # Assert length matches codon count
    expected_codon_count = len(target_nt) // 3
    assert len(result) == expected_codon_count, (
        f"Expected vector length {expected_codon_count} (codon count), "
        f"got {len(result)}"
    )


def test_calc_cars_nt_mode_returns_scalar():
    """Test calc_cARS in nucleotide mode without win_params returns a scalar.

    Expected: When working directly with nucleotide sequences (not converting to codons),
    calc_cARS still returns a scalar score.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS

    # Use nucleotide sequences directly (no nt2codon conversion)
    ref_seqs_nt = [
        "ATGAAATAA",
        "ATGGCATAA",
    ]

    target_nt = "ATGAAATAA"

    # Build suffix array from raw nucleotide sequences
    suffix_array = build_suffix_array(ref_seqs_nt)

    # Call calc_cARS with nucleotide sequence
    result = calc_cARS(target_nt, suffix_array, max_len=40, max_pos=0.5)

    # Assert it returns a scalar
    assert isinstance(result, (float, int)), f"Expected scalar, got {type(result)}"
    assert result >= 0, "cARS score should be non-negative"


def test_calc_cars_nt_mode_with_return_vec():
    """Test calc_cARS in nucleotide mode with return_vec=True returns array of nt length.

    Expected: With return_vec=True in nucleotide mode, calc_cARS returns an array
    of per-nucleotide scores. Array length should match nucleotide count.
    """
    # Import required functions from pychimera
    from pychimera import build_suffix_array, calc_cARS

    # Use nucleotide sequences directly
    ref_seqs_nt = [
        "ATGAAATAA",
        "ATGGCATAA",
    ]

    target_nt = "ATGAAATAA"

    # Build suffix array
    suffix_array = build_suffix_array(ref_seqs_nt)

    # Call calc_cARS with return_vec=True
    result = calc_cARS(target_nt, suffix_array, max_len=40, max_pos=0.5, return_vec=True)

    # Assert it returns an array-like object
    assert hasattr(result, "__len__"), "Expected array-like object"
    assert hasattr(result, "__getitem__"), "Expected indexable array-like object"

    # Assert length matches nucleotide count
    expected_nt_count = len(target_nt)
    assert len(result) == expected_nt_count, (
        f"Expected vector length {expected_nt_count} (nucleotide count), "
        f"got {len(result)}"
    )


def test_nt2codon_converts_list_of_sequences():
    """Test nt2codon converts a list of nucleotide sequences to codon representation.

    Expected: nt2codon takes a list of nucleotide strings and returns a list-like
    object of codon-encoded sequences.
    """
    from pychimera import nt2codon

    # Create nucleotide sequences (9 nt = 3 codons each)
    nt_seqs = [
        "ATGAAATAA",
        "ATGGCATAA",
    ]

    # Convert to codons
    codon_seqs = nt2codon(nt_seqs)

    # Assert returns list-like object with same number of sequences
    assert hasattr(codon_seqs, "__len__"), "Expected list-like object"
    assert hasattr(codon_seqs, "__getitem__"), "Expected indexable object"
    assert len(codon_seqs) == len(nt_seqs), (
        f"Expected {len(nt_seqs)} sequences, got {len(codon_seqs)}"
    )


def test_build_suffix_array_accepts_codon_sequences():
    """Test build_suffix_array accepts codon-encoded sequences.

    Expected: build_suffix_array can build a suffix array from codon-encoded
    sequences returned by nt2codon.
    """
    from pychimera import build_suffix_array, nt2codon

    # Create and convert sequences
    ref_seqs_nt = ["ATGAAATAA", "ATGGCATAA"]
    ref_codons = nt2codon(ref_seqs_nt)

    # Build suffix array
    suffix_array = build_suffix_array(ref_codons)

    # Assert it returns something (exact type may vary)
    assert suffix_array is not None, "build_suffix_array should return a value"


def test_build_suffix_array_accepts_nt_sequences():
    """Test build_suffix_array accepts raw nucleotide sequences.

    Expected: build_suffix_array can also build a suffix array directly from
    nucleotide sequences without codon conversion.
    """
    from pychimera import build_suffix_array

    # Use raw nucleotide sequences
    ref_seqs_nt = ["ATGAAATAA", "ATGGCATAA"]

    # Build suffix array
    suffix_array = build_suffix_array(ref_seqs_nt)

    # Assert it returns something
    assert suffix_array is not None, "build_suffix_array should return a value"
