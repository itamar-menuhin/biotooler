"""Tests for FeatureSet.compute_orf_windows_v2 with positional features."""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_store import attach_orf
from biotooler.features import FeatureSet
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class ToyResidueFeature:
    """A toy positional feature for testing - computes per-residue values."""

    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space."""
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specs for each feature key."""
        return {
            "gc": AggregationSpec(name="MEAN", aggregation_fn=np.mean),
            "count": AggregationSpec(name="SUM", aggregation_fn=np.sum),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-residue GC indicator and count."""
        seq = str(record.seq).upper()
        positions = kwargs.get("positions", None)

        if positions is not None:
            # Compute only for requested positions
            gc_vector = np.array([1.0 if seq[i] in "GC" else 0.0 for i in positions])
            count_vector = np.ones(len(positions))
        else:
            # Compute for all positions (backward compatibility)
            gc_vector = np.array([1.0 if b in "GC" else 0.0 for b in seq])
            count_vector = np.ones(len(seq))

        return {
            "gc": gc_vector,
            "count": count_vector,
        }


class ToyCodonFeature:
    """A toy positional feature for testing - computes per-codon values."""

    @property
    def position_space(self) -> PositionSpace:
        """Return CODON position space."""
        return PositionSpace.CODON

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specs for each feature key."""
        return {
            "start_a": AggregationSpec(name="MEAN", aggregation_fn=np.mean),
            "codon_index": AggregationSpec(name="SUM", aggregation_fn=np.sum),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-codon features."""
        seq = str(record.seq).upper()
        num_codons = len(seq) // 3
        positions = kwargs.get("positions", None)

        if positions is not None:
            # Compute only for requested codon positions
            start_a = np.array([1.0 if seq[i * 3] == "A" else 0.0 for i in positions])
            codon_index = np.array(positions, dtype=float)
        else:
            # Compute for all codons (backward compatibility)
            start_a = np.array([1.0 if seq[i * 3] == "A" else 0.0 for i in range(num_codons)])
            codon_index = np.arange(num_codons, dtype=float)

        return {
            "start_a": start_a,
            "codon_index": codon_index,
        }


class TestFeatureSetOrfWindowsV2Basic:
    """Basic tests for compute_orf_windows_v2."""

    def test_returns_single_row_dataframe(self):
        """Test that output is a single-row DataFrame (wide format)."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        assert result.shape[0] == 1  # Single row
        assert result.shape[1] > 3  # Has metadata + feature columns

    def test_metadata_columns_present(self):
        """Test that metadata columns are included."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test_id")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns
        assert result["record_id"].iloc[0] == "test_id"
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15

    def test_column_names_with_window_start_suffix(self):
        """Test that columns are named with window_start suffix."""
        fs = FeatureSet(ToyResidueFeature(), name="gc_feat")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        # Should have columns like GC_FEAT_gc_MEAN_0, GC_FEAT_gc_MEAN_3, GC_FEAT_gc_MEAN_6
        assert "GC_FEAT_gc_MEAN_0" in result.columns
        assert "GC_FEAT_gc_MEAN_3" in result.columns
        assert "GC_FEAT_gc_MEAN_6" in result.columns
        assert "GC_FEAT_count_SUM_0" in result.columns
        assert "GC_FEAT_count_SUM_3" in result.columns
        assert "GC_FEAT_count_SUM_6" in result.columns

    def test_column_ordering_metadata_first(self):
        """Test that metadata columns come before feature columns."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        cols = list(result.columns)
        # First three should be metadata
        assert cols[0] == "record_id"
        assert cols[1] == "orf_start"
        assert cols[2] == "orf_end"

    def test_non_positional_feature_raises_error(self):
        """Test that non-positional features raise ValueError."""

        def simple_feature(record):
            return {"length": len(record.seq)}

        fs = FeatureSet(simple_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="does not implement PositionalFeature protocol"):
            fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)


class TestFeatureSetOrfWindowsV2ResidueSpace:
    """Tests for v2 with residue-space features."""

    def test_residue_feature_aggregation(self):
        """Test that residue features are aggregated correctly."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        # Sequence: ATGAAACCCGGGTTT (15 nt)
        # GC content: [0,0,1,0,0,0,1,1,1,1,1,1,0,0,0]
        # Windows (9nt, step 3):
        #   0-9: ATGAAACCC -> [0,0,1,0,0,0,1,1,1] -> mean = 4/9 ≈ 0.444
        #   3-12: AAACCCGGG -> [0,0,0,1,1,1,1,1,1] -> mean = 6/9 ≈ 0.667
        #   6-15: CCCGGGTTT -> [1,1,1,1,1,1,0,0,0] -> mean = 6/9 ≈ 0.667
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        # Check GC mean aggregation
        assert abs(result["TEST_gc_MEAN_0"].iloc[0] - 4 / 9) < 1e-10
        assert abs(result["TEST_gc_MEAN_3"].iloc[0] - 6 / 9) < 1e-10
        assert abs(result["TEST_gc_MEAN_6"].iloc[0] - 6 / 9) < 1e-10

        # Check count sum aggregation (should be window size)
        assert result["TEST_count_SUM_0"].iloc[0] == 9.0
        assert result["TEST_count_SUM_3"].iloc[0] == 9.0
        assert result["TEST_count_SUM_6"].iloc[0] == 9.0

    def test_residue_feature_partial_windows(self):
        """Test residue features with partial windows."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        # 14 nt sequence, windows of 9 nt with step 3
        # Full windows: 0-9, 3-12
        # Partial window: 6-14 (8 nt, should be dropped with drop_partial=True)
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        result = fs.compute_orf_windows_v2(
            record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=True
        )

        # Should only have columns for windows 0 and 3
        feature_cols = [c for c in result.columns if c.startswith("TEST_")]
        assert len(feature_cols) == 4  # 2 windows × 2 keys
        assert "TEST_gc_MEAN_0" in result.columns
        assert "TEST_gc_MEAN_3" in result.columns
        assert "TEST_gc_MEAN_6" not in result.columns

    def test_residue_feature_with_drop_partial_false(self):
        """Test residue features include partial windows when drop_partial=False."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTT"), id="seq1")

        result = fs.compute_orf_windows_v2(
            record, orf=(0, 14), window_nt=9, step_nt=3, drop_partial=False
        )

        # Should have windows 0, 3, 6, 9, 12
        assert "TEST_gc_MEAN_0" in result.columns
        assert "TEST_gc_MEAN_3" in result.columns
        assert "TEST_gc_MEAN_6" in result.columns
        assert "TEST_gc_MEAN_9" in result.columns
        assert "TEST_gc_MEAN_12" in result.columns

        # Window at 12 has only 2 nt (14-12=2), so count should be 2
        assert result["TEST_count_SUM_12"].iloc[0] == 2.0


class TestFeatureSetOrfWindowsV2CodonSpace:
    """Tests for v2 with codon-space features."""

    def test_codon_feature_aggregation(self):
        """Test that codon features are aggregated correctly."""
        fs = FeatureSet(ToyCodonFeature(), name="test")
        # Sequence: ATGAAACCCGGGTTT (15 nt = 5 codons)
        # Codons: ATG, AAA, CCC, GGG, TTT
        # Start with A: [1, 1, 0, 0, 0]
        # Codon indices: [0, 1, 2, 3, 4]
        # Windows (9nt = 3 codons, step 3nt = 1 codon):
        #   0-9 (codons 0-3): ATG,AAA,CCC -> start_a mean = (1+1+0)/3 = 0.667
        #                                  -> codon_index sum = 0+1+2 = 3
        #   3-12 (codons 1-4): AAA,CCC,GGG -> start_a mean = (1+0+0)/3 = 0.333
        #                                   -> codon_index sum = 1+2+3 = 6
        #   6-15 (codons 2-5): CCC,GGG,TTT -> start_a mean = (0+0+0)/3 = 0.0
        #                                   -> codon_index sum = 2+3+4 = 9
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        # Check start_a mean aggregation
        assert abs(result["TEST_start_a_MEAN_0"].iloc[0] - 2 / 3) < 1e-10
        assert abs(result["TEST_start_a_MEAN_3"].iloc[0] - 1 / 3) < 1e-10
        assert abs(result["TEST_start_a_MEAN_6"].iloc[0] - 0.0) < 1e-10

        # Check codon_index sum aggregation
        assert result["TEST_codon_index_SUM_0"].iloc[0] == 3.0
        assert result["TEST_codon_index_SUM_3"].iloc[0] == 6.0
        assert result["TEST_codon_index_SUM_6"].iloc[0] == 9.0

    def test_codon_feature_partial_windows(self):
        """Test codon features handle partial windows correctly."""
        fs = FeatureSet(ToyCodonFeature(), name="test")
        # 18 nt = 6 codons, windows of 9 nt (3 codons), step 3 nt (1 codon)
        # Windows: 0-9 (codons 0-3), 3-12 (codons 1-4), 6-15 (codons 2-5), 9-18 (codons 3-6)
        # With drop_partial=True, all windows fit
        record = SeqRecord(Seq("ATGAAACCCGGGTTTAAA"), id="seq1")

        result = fs.compute_orf_windows_v2(
            record, orf=(0, 18), window_nt=9, step_nt=3, drop_partial=True
        )

        # Should have 4 windows
        feature_cols = [c for c in result.columns if c.startswith("TEST_")]
        assert len(feature_cols) == 8  # 4 windows × 2 keys
        assert "TEST_start_a_MEAN_0" in result.columns
        assert "TEST_start_a_MEAN_3" in result.columns
        assert "TEST_start_a_MEAN_6" in result.columns
        assert "TEST_start_a_MEAN_9" in result.columns


class TestFeatureSetOrfWindowsV2OrfResolution:
    """Tests for ORF resolution in v2."""

    def test_explicit_orf_parameter(self):
        """Test ORF resolution using explicit orf parameter."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("NNNNATGAAACCCGGGTTTNNNN"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(4, 19), window_nt=9, step_nt=3)

        assert result["orf_start"].iloc[0] == 4
        assert result["orf_end"].iloc[0] == 19

    def test_orf_index_parameter(self):
        """Test ORF resolution using orf_index."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAATAGATGCCCTAGTAA"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf_index=0, window_nt=6, step_nt=3)

        # First candidate should be (0, 9)
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 9

    def test_attached_orf(self):
        """Test ORF resolution using attached ORF."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")
        attach_orf(record, (0, 15))

        result = fs.compute_orf_windows_v2(record, window_nt=9, step_nt=3)

        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15

    def test_missing_orf_raises_error(self):
        """Test that missing ORF raises ValueError."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="No ORF information provided"):
            fs.compute_orf_windows_v2(record, window_nt=9, step_nt=3)


class TestFeatureSetOrfWindowsV2Validation:
    """Tests for validation in v2."""

    def test_protein_sequence_raises_error(self):
        """Test that protein sequences raise ValueError."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("MKLVLS"), id="protein1")
        record.annotations["molecule_type"] = "protein"

        with pytest.raises(ValueError, match="only supported for DNA/RNA sequences, not protein"):
            fs.compute_orf_windows_v2(record, orf=(0, 6), window_nt=3, step_nt=3)

    def test_step_not_multiple_of_3_raises_error(self):
        """Test that step_nt not multiple of 3 raises ValueError."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=4)


class TestFeatureSetOrfWindowsV2EmptyCases:
    """Tests for edge cases with no windows."""

    def test_no_windows_returns_metadata_only(self):
        """Test that when no windows fit, returns row with metadata only."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAA"), id="seq1")

        # Window size larger than ORF with drop_partial=True
        result = fs.compute_orf_windows_v2(
            record, orf=(0, 6), window_nt=12, step_nt=3, drop_partial=True
        )

        # Should have only metadata columns
        assert result.shape[0] == 1
        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns
        assert len([c for c in result.columns if c.startswith("TEST_")]) == 0


class TestFeatureSetOrfWindowsV2Sorting:
    """Tests for column sorting in v2."""

    def test_numeric_sorting_of_window_starts(self):
        """Test that window starts are sorted numerically."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        # Long ORF to get window_start >= 10
        long_seq = "ATG" + "AAA" * 20  # 63 nt
        record = SeqRecord(Seq(long_seq), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 63), window_nt=9, step_nt=3)

        feature_cols = [c for c in result.columns if c.startswith("TEST_")]

        # The columns should be sorted by feature key, then by window_start
        # Since ToyResidueFeature has "count" and "gc" keys (alphabetically sorted),
        # we expect: count_0, count_3, ..., count_54, gc_0, gc_3, ..., gc_54

        # Verify we have double-digit window_starts
        window_starts = [int(c.rsplit("_", 1)[1]) for c in feature_cols]
        assert any(ws >= 10 for ws in window_starts)

        # Verify overall column ordering is correct
        # Extract (key, window_start) tuples
        key_window_pairs = []
        for col in feature_cols:
            parts = col.split(".")[-1].rsplit("_", 1)  # Get "key_window" and split
            key = parts[0]
            window = int(parts[1])
            key_window_pairs.append((key, window))

        # Check that they're sorted by key first, then window
        assert key_window_pairs == sorted(key_window_pairs)

    def test_features_sorted_by_key_then_window(self):
        """Test that features are sorted by key name then window_start."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=3)

        feature_cols = [c for c in result.columns if c.startswith("TEST_")]

        # Features: "count" and "gc" (alphabetically: count < gc)
        # Windows: 0, 3, 6
        # Expected order: count_0, count_3, count_6, gc_0, gc_3, gc_6
        expected = [
            "TEST_count_SUM_0",
            "TEST_count_SUM_3",
            "TEST_count_SUM_6",
            "TEST_gc_MEAN_0",
            "TEST_gc_MEAN_3",
            "TEST_gc_MEAN_6",
        ]
        assert feature_cols == expected
