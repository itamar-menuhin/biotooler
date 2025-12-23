"""Tests for ViennaRNA window fold caching functionality."""

from unittest.mock import MagicMock, patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.families.viennarna.cache import (
    FoldKey,
    ViennaFoldCache,
    _hash_string,
    get_context_mfe_cached,
)


class TestFoldKey:
    """Tests for FoldKey dataclass."""

    def test_fold_key_frozen(self):
        """Test that FoldKey is immutable (frozen)."""
        key = FoldKey(
            seq_hash="abc123",
            flank_left=5,
            flank_right=5,
        )

        # Should not be able to modify
        with pytest.raises(AttributeError):
            key.flank_left = 10  # type: ignore[misc]

    def test_fold_key_equality(self):
        """Test that FoldKeys with same values are equal."""
        key1 = FoldKey(
            seq_hash="abc",
            flank_left=5,
            flank_right=5,
        )
        key2 = FoldKey(
            seq_hash="abc",
            flank_left=5,
            flank_right=5,
        )

        assert key1 == key2
        assert hash(key1) == hash(key2)

    def test_fold_key_different_values_not_equal(self):
        """Test that FoldKeys with different values are not equal."""
        key1 = FoldKey(
            seq_hash="abc",
            flank_left=5,
            flank_right=5,
        )
        key2 = FoldKey(
            seq_hash="xyz",  # Different sequence
            flank_left=5,
            flank_right=5,
        )

        assert key1 != key2


class TestViennaFoldCache:
    """Tests for ViennaFoldCache class."""

    def test_cache_init_empty(self):
        """Test that cache initializes empty."""
        cache = ViennaFoldCache()
        assert len(cache) == 0

    def test_cache_get_miss(self):
        """Test that get returns None for cache miss."""
        cache = ViennaFoldCache()
        key = FoldKey(
            seq_hash="abc",
            flank_left=0,
            flank_right=0,
        )

        assert cache.get(key) is None

    def test_cache_set_and_get(self):
        """Test that set stores value and get retrieves it."""
        cache = ViennaFoldCache()
        key = FoldKey(
            seq_hash="abc",
            flank_left=0,
            flank_right=0,
        )
        structure = "(((...)))"
        mfe = -5.2

        cache.set(key, structure, mfe)

        result = cache.get(key)
        assert result is not None
        assert result == (structure, mfe)
        assert len(cache) == 1

    def test_cache_clear(self):
        """Test that clear removes all entries."""
        cache = ViennaFoldCache()
        key = FoldKey(
            seq_hash="abc",
            flank_left=0,
            flank_right=0,
        )

        cache.set(key, "(((...)))", -5.2)
        assert len(cache) == 1

        cache.clear()
        assert len(cache) == 0
        assert cache.get(key) is None

    def test_cache_multiple_entries(self):
        """Test that cache can store multiple entries."""
        cache = ViennaFoldCache()

        key1 = FoldKey(
            seq_hash="abc",
            flank_left=0,
            flank_right=0,
        )
        key2 = FoldKey(
            seq_hash="xyz",
            flank_left=0,
            flank_right=0,
        )

        cache.set(key1, "(((...)))", -5.2)
        cache.set(key2, ".........", -0.5)

        assert len(cache) == 2
        assert cache.get(key1) == ("(((...)))", -5.2)
        assert cache.get(key2) == (".........", -0.5)


class TestHashFunctions:
    """Tests for hash utility functions."""

    def test_hash_string_deterministic(self):
        """Test that hash is deterministic."""
        s = "ACGUACGU"
        hash1 = _hash_string(s)
        hash2 = _hash_string(s)

        assert hash1 == hash2
        assert len(hash1) == 64  # SHA256 produces 64 hex chars

    def test_hash_string_different_inputs(self):
        """Test that different strings produce different hashes."""
        hash1 = _hash_string("ACGU")
        hash2 = _hash_string("UGCA")

        assert hash1 != hash2


class TestGetContextMfeCached:
    """Tests for get_context_mfe_cached function."""

    def test_caching_prevents_redundant_folds(self):
        """Test that caching prevents redundant fold computations."""
        # Create mock RNA module
        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()
            seq_str = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bp

            # First call - should compute
            structure1, mfe1 = get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=cache,
            )

            assert len(fold_compound_calls) == 1
            assert structure1 == "." * 20
            assert mfe1 == -10.0

            # Second call with same parameters - should use cache
            structure2, mfe2 = get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=cache,
            )

            # Should NOT have made another fold_compound call
            assert len(fold_compound_calls) == 1, (
                "Cache should have prevented second fold computation"
            )
            assert structure2 == structure1
            assert mfe2 == mfe1

    def test_different_windows_not_cached_together(self):
        """Test that different windows are not confused in cache."""
        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()
            seq_str = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bp

            # First window
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=cache,
            )

            # Different window (different start position)
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=10,
                window_size=20,
                cache=cache,
            )

            # Should have made two fold_compound calls (different sequences)
            assert len(fold_compound_calls) == 2

    def test_identical_sequences_use_cache(self):
        """Test that identical sequences from different positions use cache."""
        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()
            # Create sequence with repeating pattern
            pattern = "ACGTACGTACGTACGTACGT"  # 20 bp
            seq_str = pattern * 2  # 40 bp

            # First window 0-20
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=cache,
            )

            assert len(fold_compound_calls) == 1

            # Second window 20-40 (identical sequence)
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=20,
                window_size=20,
                cache=cache,
            )

            # Should use cache since sequences are identical
            assert len(fold_compound_calls) == 1, (
                "Identical sequences should use cache"
            )

    def test_flanks_affect_cache_key(self):
        """Test that different flank configurations produce different cache keys."""
        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()
            seq_str = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bp

            # First call with no flanks
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=10,
                window_size=10,
                flank_left=0,
                flank_right=0,
                cache=cache,
            )

            assert len(fold_compound_calls) == 1

            # Second call with flanks - should recompute
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=10,
                window_size=10,
                flank_left=5,
                flank_right=5,
                cache=cache,
            )

            # Should have made another fold_compound call due to different flanks
            assert len(fold_compound_calls) == 2

    def test_dna_to_rna_conversion(self):
        """Test that DNA (T) is converted to RNA (U) before folding."""
        rna_mock = MagicMock()
        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            seq_str = "ACGTACGTACGTACGTACGT"  # DNA with T

            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
            )

            # Verify that sequence passed to fold_compound has U not T
            assert len(sequences_folded) == 1
            assert "T" not in sequences_folded[0]
            assert "U" in sequences_folded[0]
            assert sequences_folded[0] == "ACGUACGUACGUACGUACGU"

    def test_window_beyond_sequence_raises_error(self):
        """Test that window extending beyond sequence raises error."""
        rna_mock = MagicMock()

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            seq_str = "ACGTACGTACGT"  # 12 bp

            with pytest.raises(ValueError, match="Window extends beyond sequence"):
                get_context_mfe_cached(
                    seq_str=seq_str,
                    window_start=5,
                    window_size=20,  # Extends beyond sequence
                )

    def test_boundary_normalization_left_flank(self):
        """Test that left flank is normalized when window starts near beginning."""
        rna_mock = MagicMock()
        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            seq_str = "ACGTACGTACGTACGTACGT"  # 20 bp

            # Window at start with requested left flank that extends before 0
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=2,
                window_size=10,
                flank_left=5,  # Would go to position -3, should clip to 0
            )

            # Should have folded from position 0, not -3
            # Context: positions 0-12 (start=2-5=0 clipped, end=2+10=12)
            assert len(sequences_folded[0]) == 12

    def test_boundary_normalization_right_flank(self):
        """Test that right flank is normalized when extending beyond sequence."""
        rna_mock = MagicMock()
        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            seq_str = "ACGTACGTACGTACGTACGT"  # 20 bp

            # Window ending near sequence end with flank extending beyond
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=5,
                window_size=10,
                flank_right=10,  # Would go to position 25, should clip to 20
            )

            # Should have folded to position 20, not 25
            # Context: positions 5-20 (start=5, end=5+10+10=25 clipped to 20)
            assert len(sequences_folded[0]) == 15

    def test_no_cache_always_recomputes(self):
        """Test that without cache, every call recomputes."""
        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            seq_str = "ACGTACGTACGTACGTACGTACGTACGTACGT"  # 32 bp

            # Two calls with same parameters but no cache
            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=None,  # No cache
            )

            get_context_mfe_cached(
                seq_str=seq_str,
                window_start=0,
                window_size=20,
                cache=None,  # No cache
            )

            # Should have computed both times
            assert len(fold_compound_calls) == 2


class TestWindowMFEFeatureWithCache:
    """Test WindowMFEFeature with caching."""

    def test_window_mfe_uses_cache(self):
        """Test that WindowMFEFeature uses cache when provided."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()

            # Create feature with cache
            feature = WindowMFEFeature(
                window_starts=[0, 10, 20],
                window_size=10,
                cache=cache,
            )

            # Create sequence with repeating pattern that would result in identical substrings
            seq_str = "ACGTACGTAC" * 4  # 40 bp with repeating pattern
            record = SeqRecord(Seq(seq_str), id="test")

            result = feature(record)

            # All three windows should have been computed
            assert "MFE_0" in result
            assert "MFE_10" in result
            assert "MFE_20" in result

            # But the repeating pattern means some are identical
            # Windows at 0, 10, 20 with size 10 from "ACGTACGTAC" repeated:
            # Window 0-10: ACGTACGTAC
            # Window 10-20: ACGTACGTAC (identical!)
            # Window 20-30: ACGTACGTAC (identical!)
            # So only 1 unique fold computation should occur
            assert len(fold_compound_calls) == 1, (
                f"Expected 1 fold computation due to cache, got {len(fold_compound_calls)}"
            )

    def test_window_mfe_without_cache(self):
        """Test that WindowMFEFeature works without cache (always recomputes)."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            # Create feature without cache
            feature = WindowMFEFeature(
                window_starts=[0, 10, 20],
                window_size=10,
                cache=None,  # No cache
            )

            seq_str = "ACGTACGTAC" * 4  # Same repeating pattern
            record = SeqRecord(Seq(seq_str), id="test")

            feature(record)

            # Without cache, all 3 windows should be computed even if identical
            assert len(fold_compound_calls) == 3


class TestContextWindowFoldFeatureWithCache:
    """Test ContextWindowFoldFeature with caching."""

    def test_context_window_fold_uses_cache(self):
        """Test that ContextWindowFoldFeature uses cache when provided."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            # Create simple structure for testing
            structure = "." * len(seq)
            fc_mock.mfe.return_value = (structure, mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()

            # Create feature with cache and flanks
            feature = ContextWindowFoldFeature(
                starts_nt=[0, 15, 30],
                window_size_nt=10,
                flank_left_nt=5,
                flank_right_nt=5,
                cache=cache,
            )

            # Create sequence with repeating pattern
            seq_str = "ACGTACGTACGTACGTAC" * 3  # 54 bp
            record = SeqRecord(Seq(seq_str), id="test")

            result = feature(record)

            # All three windows should have been computed
            assert "CTX_MFE_0" in result
            assert "CTX_MFE_15" in result
            assert "CTX_MFE_30" in result

            # Context slices:
            # Window 0: context 0-15 (no left flank at start)
            # Window 15: context 10-30
            # Window 30: context 25-45
            # All different sequences, so 3 fold computations
            assert len(fold_compound_calls) == 3

    def test_context_window_fold_cache_with_identical_contexts(self):
        """Test caching when multiple windows have identical context slices."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        rna_mock = MagicMock()
        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            structure = "." * len(seq)
            fc_mock.mfe.return_value = (structure, mfe_value)
            return fc_mock

        rna_mock.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=rna_mock,
        ):
            cache = ViennaFoldCache()

            # Create feature with NO flanks - simpler cache key matching
            feature = ContextWindowFoldFeature(
                starts_nt=[0, 20, 40],
                window_size_nt=20,
                flank_left_nt=0,
                flank_right_nt=0,
                cache=cache,
            )

            # Create sequence with exact repeating pattern
            pattern = "ACGTACGTACGTACGTACGT"  # 20 bp
            seq_str = pattern * 3  # 60 bp, three identical 20bp windows
            record = SeqRecord(Seq(seq_str), id="test")

            result = feature(record)

            # All three windows should be in result
            assert "CTX_MFE_0" in result
            assert "CTX_MFE_20" in result
            assert "CTX_MFE_40" in result

            # But only 1 fold computation due to identical sequences and cache
            assert len(fold_compound_calls) == 1, (
                f"Expected 1 fold computation due to cache and identical sequences, "
                f"got {len(fold_compound_calls)}"
            )
