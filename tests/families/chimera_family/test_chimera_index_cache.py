"""Tests for chimera index caching functionality.

Verifies that ChimeraIndexKey, ChimeraReferenceIndex, and get_reference_index
provide deterministic caching with alphabet-specific keys.
"""

import pytest

from biotooler.core.reference_sequences import ReferenceSequenceSet


def test_chimera_index_key_creation():
    """Test that ChimeraIndexKey can be created and is hashable."""
    try:
        from biotooler.families.chimera.index import ChimeraIndexKey
    except ImportError:
        pytest.skip("pychimera not installed")

    key1 = ChimeraIndexKey(alphabet="codon", genetic_code_table=1, ref_hash="abc123")
    key2 = ChimeraIndexKey(alphabet="codon", genetic_code_table=1, ref_hash="abc123")
    key3 = ChimeraIndexKey(alphabet="aa", genetic_code_table=1, ref_hash="abc123")

    # Same keys should be equal
    assert key1 == key2

    # Different alphabet should create different key
    assert key1 != key3

    # Should be hashable (for use as dict key)
    assert hash(key1) == hash(key2)
    assert hash(key1) != hash(key3)


def test_chimera_index_key_string_representation():
    """Test that ChimeraIndexKey has readable string representation."""
    try:
        from biotooler.families.chimera.index import ChimeraIndexKey
    except ImportError:
        pytest.skip("pychimera not installed")

    key = ChimeraIndexKey(
        alphabet="codon",
        genetic_code_table=1,
        ref_hash="abcdef1234567890abcdef1234567890",
    )
    key_str = str(key)

    assert "codon" in key_str
    assert "gc1" in key_str
    assert "abcdef1234567890" in key_str


def test_get_reference_index_without_cache():
    """Test that get_reference_index builds index without cache."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a simple reference set
    cds = {
        "gene1": "ATGATGATG",  # M M M
        "gene2": "ATGGCATAA",  # M A *
    }
    ref_set = ReferenceSequenceSet(cds)

    # Build index without cache
    try:
        index = get_reference_index(ref_set, alphabet="codon")
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should return a ChimeraReferenceIndex
    assert index is not None
    assert index.suffix_array is not None
    assert index.key is not None
    assert index.key.alphabet == "codon"
    assert index.key.genetic_code_table == 1


def test_get_reference_index_with_cache_hit():
    """Test that get_reference_index reuses cached index."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a simple reference set
    cds = {
        "gene1": "ATGATGATG",
        "gene2": "ATGGCATAA",
    }
    ref_set = ReferenceSequenceSet(cds)

    # Build index with cache
    cache = {}
    try:
        index1 = get_reference_index(ref_set, alphabet="codon", cache=cache)
    except ImportError:
        pytest.skip("pychimera not installed")

    # Cache should now contain one entry
    assert len(cache) == 1

    # Build again with same parameters - should hit cache
    index2 = get_reference_index(ref_set, alphabet="codon", cache=cache)

    # Should be the exact same object (cache hit)
    assert index1 is index2
    assert len(cache) == 1


def test_get_reference_index_with_cache_miss_different_alphabet():
    """Test that different alphabets create different cache entries."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set with both CDS and proteins
    cds = {
        "gene1": "ATGATGATG",
        "gene2": "ATGGCATAA",
    }
    proteins = {
        "gene1": "MMM",
        "gene2": "MA",
    }
    ref_set = ReferenceSequenceSet(cds, proteins)

    # Build indices with different alphabets
    cache = {}
    try:
        index_codon = get_reference_index(ref_set, alphabet="codon", cache=cache)
        index_aa = get_reference_index(ref_set, alphabet="aa", cache=cache)
        index_nt = get_reference_index(ref_set, alphabet="nt", cache=cache)
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should create three different cache entries
    assert len(cache) == 3

    # Should be different objects
    assert index_codon is not index_aa
    assert index_codon is not index_nt
    assert index_aa is not index_nt

    # Keys should have different alphabets
    assert index_codon.key.alphabet == "codon"
    assert index_aa.key.alphabet == "aa"
    assert index_nt.key.alphabet == "nt"


def test_get_reference_index_with_cache_miss_different_genetic_code():
    """Test that different genetic code tables create different cache entries."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create a reference set
    cds = {
        "gene1": "ATGATGATG",
        "gene2": "ATGGCATAA",
    }
    ref_set = ReferenceSequenceSet(cds, genetic_code_table=1)

    # Build indices with different genetic code tables
    cache = {}
    try:
        index_gc1 = get_reference_index(
            ref_set, alphabet="codon", genetic_code_table=1, cache=cache
        )
        index_gc11 = get_reference_index(
            ref_set, alphabet="codon", genetic_code_table=11, cache=cache
        )
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should create two different cache entries
    assert len(cache) == 2

    # Should be different objects
    assert index_gc1 is not index_gc11

    # Keys should have different genetic code tables
    assert index_gc1.key.genetic_code_table == 1
    assert index_gc11.key.genetic_code_table == 11


def test_get_reference_index_deterministic_hash():
    """Test that reference hash is deterministic for same sequences."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create two reference sets with same sequences but different order
    cds1 = {
        "gene1": "ATGATGATG",
        "gene2": "ATGGCATAA",
    }
    cds2 = {
        "gene2": "ATGGCATAA",
        "gene1": "ATGATGATG",
    }
    ref_set1 = ReferenceSequenceSet(cds1)
    ref_set2 = ReferenceSequenceSet(cds2)

    # Build indices
    cache1 = {}
    cache2 = {}
    try:
        index1 = get_reference_index(ref_set1, alphabet="codon", cache=cache1)
        index2 = get_reference_index(ref_set2, alphabet="codon", cache=cache2)
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should have same hash (deterministic despite different order)
    assert index1.key.ref_hash == index2.key.ref_hash

    # Should be able to hit cache with either reference set
    shared_cache = {}
    index_a = get_reference_index(ref_set1, alphabet="codon", cache=shared_cache)
    index_b = get_reference_index(ref_set2, alphabet="codon", cache=shared_cache)

    # Should hit cache (same object)
    assert index_a is index_b
    assert len(shared_cache) == 1


def test_get_reference_index_different_sequences_different_hash():
    """Test that different sequences produce different hashes."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    # Create two reference sets with different sequences
    cds1 = {
        "gene1": "ATGATGATG",
    }
    cds2 = {
        "gene1": "ATGATGAAA",
    }
    ref_set1 = ReferenceSequenceSet(cds1)
    ref_set2 = ReferenceSequenceSet(cds2)

    # Build indices in shared cache
    cache = {}
    try:
        index1 = get_reference_index(ref_set1, alphabet="codon", cache=cache)
        index2 = get_reference_index(ref_set2, alphabet="codon", cache=cache)
    except ImportError:
        pytest.skip("pychimera not installed")

    # Should have different hashes
    assert index1.key.ref_hash != index2.key.ref_hash

    # Should create two cache entries
    assert len(cache) == 2


def test_get_reference_index_invalid_alphabet():
    """Test that invalid alphabet raises ValueError."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    cds = {"gene1": "ATGATGATG"}
    ref_set = ReferenceSequenceSet(cds)

    # Should raise ValueError for invalid alphabet
    try:
        with pytest.raises(ValueError, match="alphabet must be"):
            get_reference_index(ref_set, alphabet="invalid")
    except ImportError:
        pytest.skip("pychimera not installed")


def test_get_reference_index_cache_is_optional():
    """Test that cache parameter is optional."""
    try:
        from biotooler.families.chimera.index import get_reference_index
    except ImportError:
        pytest.skip("pychimera not installed")

    cds = {"gene1": "ATGATGATG"}
    ref_set = ReferenceSequenceSet(cds)

    # Should work without cache parameter
    try:
        index1 = get_reference_index(ref_set, alphabet="codon")
        index2 = get_reference_index(ref_set, alphabet="codon")
    except ImportError:
        pytest.skip("pychimera not installed")

    # Without cache, should build new indices each time (different objects)
    # Note: Can't use 'is' comparison as implementation detail, but can verify they work
    assert index1 is not None
    assert index2 is not None
