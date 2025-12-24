"""Tests for FeatureSet ORF window computation (wide format)."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_store import attach_orf
from biotooler.features import FeatureSet


class TestFeatureSetOrfWindowsWideFormat:
    """Tests for FeatureSet.compute_orf_windows with wide-format output."""

    def test_returns_single_row_dataframe(self):
        """Test that output is a single-row DataFrame (wide format)."""

        def simple_feature(record):
            return {"length": len(record.seq)}

        fs = FeatureSet(simple_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        assert result.shape[0] == 1  # Single row (wide format)
        assert result.shape[1] > 3  # Has metadata + feature columns

    def test_wide_format_column_names_with_suffixes(self):
        """Test that columns are named with window_start suffix (_0, _3, _6)."""

        def count_a(record):
            return {"a_count": str(record.seq).count("A")}

        fs = FeatureSet(count_a, name="count")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        # Should have columns like count.a_count_0, count.a_count_3, count.a_count_6
        assert "count.a_count_0" in result.columns
        assert "count.a_count_3" in result.columns
        assert "count.a_count_6" in result.columns

    def test_metadata_columns_present(self):
        """Test that metadata columns (record_id, orf_start, orf_end) are included."""

        def dummy_feature(record):
            return {"dummy": 0}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test_id")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns
        assert result["record_id"].iloc[0] == "test_id"
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15

    def test_column_ordering_metadata_first(self):
        """Test that metadata columns come before feature columns."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        cols = list(result.columns)
        # First three should be metadata
        assert cols[0] == "record_id"
        assert cols[1] == "orf_start"
        assert cols[2] == "orf_end"

    def test_column_ordering_features_sorted_by_key_then_window(self):
        """Test that feature columns are sorted by feature key, then window_start."""

        def multi_feature(record):
            return {"a_feat": 1, "z_feat": 2}

        fs = FeatureSet(multi_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        feature_cols = [c for c in result.columns if c.startswith("test.")]

        # Should be sorted by feature key, then window_start
        # Expected: a_feat_0, a_feat_3, a_feat_6, z_feat_0, z_feat_3, z_feat_6
        assert feature_cols[0] == "test.a_feat_0"
        assert feature_cols[1] == "test.a_feat_3"
        assert feature_cols[2] == "test.a_feat_6"
        assert feature_cols[3] == "test.z_feat_0"
        assert feature_cols[4] == "test.z_feat_3"
        assert feature_cols[5] == "test.z_feat_6"

    def test_feature_values_computed_correctly(self):
        """Test that feature values are computed correctly for each window."""

        def count_g(record):
            return {"g_count": str(record.seq).count("G")}

        fs = FeatureSet(count_g, name="test")
        # Windows: ATGAAACCC (1 G), AAACCCGGG (3 G), CCCGGGTTT (3 G)
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        assert result["test.g_count_0"].iloc[0] == 1
        assert result["test.g_count_3"].iloc[0] == 3
        assert result["test.g_count_6"].iloc[0] == 3

    def test_orf_resolution_explicit_orf_parameter(self):
        """Test ORF resolution using explicit orf parameter."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("NNNNATGAAACCCGGGTTTNNNN"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(4, 19), window_nt=9, step_nt=3)

        assert result["orf_start"].iloc[0] == 4
        assert result["orf_end"].iloc[0] == 19

    def test_orf_resolution_orf_index_parameter(self):
        """Test ORF resolution using orf_index to select from candidates."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        # Sequence with multiple ORFs: ATGAAATAA (0,9), longer ones...
        record = SeqRecord(Seq("ATGAAATAGATGCCCTAGTAA"), id="seq1")

        # Select first ORF candidate (index 0)
        result = fs.compute_orf_windows(record, orf_index=0, window_nt=6, step_nt=3)

        # First candidate should be (0, 9)
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 9

    def test_orf_resolution_attached_orf(self):
        """Test ORF resolution using attached ORF from record."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
        attach_orf(record, (0, 15))

        # No orf or orf_index provided - should use attached ORF
        result = fs.compute_orf_windows(record, window_nt=9, step_nt=3)

        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15

    def test_missing_orf_raises_clear_error(self):
        """Test that missing ORF info raises ValueError with clear message."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
        # No ORF attached, no orf or orf_index provided

        with pytest.raises(ValueError, match="No ORF information provided"):
            fs.compute_orf_windows(record, window_nt=9, step_nt=3)

    def test_protein_sequence_raises_clear_error(self):
        """Test that protein sequences raise ValueError with clear message."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("MKLVLS"), id="protein1")
        record.annotations["molecule_type"] = "protein"

        with pytest.raises(ValueError, match="only supported for DNA/RNA sequences, not protein"):
            fs.compute_orf_windows(record, orf=(0, 6), window_nt=3, step_nt=3)

    def test_step_not_multiple_of_3_raises_error(self):
        """Test that step_nt not multiple of 3 raises ValueError."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=4)

    def test_no_windows_returns_metadata_only(self):
        """Test that when no windows fit, returns row with metadata only."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        record = SeqRecord(Seq("ATGAAA"), id="seq1")

        # Window size larger than ORF with drop_partial=True
        result = fs.compute_orf_windows(
            record, orf=(0, 6), window_nt=12, step_nt=3, drop_partial=True
        )

        # Should have only metadata columns
        assert result.shape[0] == 1
        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns
        assert len([c for c in result.columns if c.startswith("test.")]) == 0

    def test_drop_partial_true_excludes_partial_windows(self):
        """Test that partial windows are excluded when drop_partial=True."""

        def count_a(record):
            return {"a_count": str(record.seq).count("A")}

        fs = FeatureSet(count_a, name="test")
        # 14 nt, windows of 9 nt with step 3 nt -> only 2 full windows fit
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        result = fs.compute_orf_windows(
            record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=True
        )

        # Should only have columns for windows 0 and 3 (full windows)
        feature_cols = [c for c in result.columns if c.startswith("test.")]
        assert len(feature_cols) == 2
        assert "test.a_count_0" in result.columns
        assert "test.a_count_3" in result.columns
        assert "test.a_count_6" not in result.columns

    def test_drop_partial_false_includes_partial_windows(self):
        """Test that partial windows are included when drop_partial=False."""

        def count_a(record):
            return {"a_count": str(record.seq).count("A")}

        fs = FeatureSet(count_a, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        result = fs.compute_orf_windows(
            record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=False
        )

        # Should have columns for all windows including partial ones
        feature_cols = [c for c in result.columns if c.startswith("test.")]
        assert len(feature_cols) == 5
        assert "test.a_count_0" in result.columns
        assert "test.a_count_3" in result.columns
        assert "test.a_count_6" in result.columns
        assert "test.a_count_9" in result.columns
        assert "test.a_count_12" in result.columns

    def test_multiple_features_per_window(self):
        """Test computing multiple features per window."""

        def multi_features(record):
            seq = str(record.seq)
            return {
                "length": len(seq),
                "gc_count": seq.count("G") + seq.count("C"),
                "a_count": seq.count("A"),
            }

        fs = FeatureSet(multi_features, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        # Should have 3 features × 3 windows = 9 feature columns
        feature_cols = [c for c in result.columns if c.startswith("test.")]
        assert len(feature_cols) == 9

        # Check that all feature types are present for each window
        for window_start in [0, 3, 6]:
            assert f"test.a_count_{window_start}" in result.columns
            assert f"test.gc_count_{window_start}" in result.columns
            assert f"test.length_{window_start}" in result.columns

    def test_numeric_sorting_of_window_starts(self):
        """Test that window starts are sorted numerically (3, 6, 9, 12 not 12, 3, 6, 9)."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="test")
        # Long ORF to get window_start >= 10
        long_seq = "ATG" + "AAA" * 20  # 63 nt
        record = SeqRecord(Seq(long_seq), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 63), window_nt=9, step_nt=3)

        feature_cols = [c for c in result.columns if c.startswith("test.")]

        # Extract window_starts
        window_starts = [int(c.rsplit("_", 1)[1]) for c in feature_cols]

        # Should be sorted numerically
        assert window_starts == sorted(window_starts)
        # Verify we have double-digit window_starts
        assert any(ws >= 10 for ws in window_starts)

    def test_featureset_name_in_column_names(self):
        """Test that FeatureSet name appears in column names."""

        def dummy_feature(record):
            return {"feat": 1}

        fs = FeatureSet(dummy_feature, name="custom_name")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

        feature_cols = [c for c in result.columns if c.startswith("custom_name.")]
        assert len(feature_cols) == 3
        assert "custom_name.feat_0" in result.columns
        assert "custom_name.feat_3" in result.columns
        assert "custom_name.feat_6" in result.columns


class ToyPositionalCodonFeature:
    """Toy positional feature that returns codon indices [0, 1, 2, ...] for testing."""

    @property
    def position_space(self):
        """Return CODON position space."""
        from biotooler.features.aggregation import PositionSpace

        return PositionSpace.CODON

    @property
    def vector_keys(self):
        """Return aggregation specs - using mean aggregation."""
        import numpy as np

        from biotooler.features.aggregation import AggregationSpec

        return {"codon_idx": AggregationSpec(aggregation_fn=np.mean)}

    def compute_vector(self, record, **kwargs):
        """Return array [0, 1, 2, ...] for each codon."""
        import numpy as np

        seq_len = len(record.seq)
        num_codons = seq_len // 3
        return {"codon_idx": np.arange(num_codons, dtype=float)}


def test_step_nt_does_not_subsample_codons():
    """Regression test: step_nt moves window but aggregates all codons in range.

    This test proves that step_nt does NOT subsample codons/residues; it only
    moves the window start position. Each window aggregates all positions in
    its range [window_start, window_end).

    Given:
    - 18 nt sequence (6 codons: indices 0, 1, 2, 3, 4, 5)
    - window_nt=12 (4 codons)
    - step_nt=6 (1 codon step)

    Expected windows:
    - Window at 0: codons 0-4 (nt 0-12), mean = (0+1+2+3)/4 = 1.5
    - Window at 6: codons 2-6 (nt 6-18), mean = (2+3+4+5)/4 = 3.5

    If step_nt incorrectly subsampled, window 6 might only aggregate
    codons 2 and 4 (every other codon), giving mean=3.0 instead of 3.5.
    """
    feature = ToyPositionalCodonFeature()
    fs = FeatureSet(feature, name="test")

    # 18 nt = 6 codons
    record = SeqRecord(Seq("ATGATGATGATGATGATG"), id="test")

    result = fs.compute_orf_windows(record, orf=(0, 18), window_nt=12, step_nt=6, drop_partial=True)

    # Window 0: codons 0, 1, 2, 3 -> mean = (0+1+2+3)/4 = 1.5
    assert "test.codon_idx_0" in result.columns
    assert abs(result["test.codon_idx_0"].iloc[0] - 1.5) < 1e-10

    # Window 6: codons 2, 3, 4, 5 -> mean = (2+3+4+5)/4 = 3.5
    assert "test.codon_idx_6" in result.columns
    assert abs(result["test.codon_idx_6"].iloc[0] - 3.5) < 1e-10
