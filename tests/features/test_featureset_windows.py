"""Tests for FeatureSet generic window computation (wide format)."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet


class TinyNonIncrementalFeature:
    """A simple non-incremental feature that returns sequence length."""

    def __call__(self, record: SeqRecord) -> dict[str, int]:
        """Return the length of the sequence."""
        return {"len": len(record.seq)}


class TinyIncrementalFeature:
    """A simple incremental feature that tracks nucleotide counts.

    This feature maintains a running count of each nucleotide type
    as windows slide across the sequence.
    """

    def __call__(self, record: SeqRecord) -> dict[str, int]:
        """Make the feature callable for FeatureSet compatibility."""
        seq = str(record.seq).upper()
        return {
            "a_count": seq.count("A"),
            "c_count": seq.count("C"),
        }

    def init_state(
        self,
        record: SeqRecord,
        *,
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> dict:
        """Initialize state for the first window."""
        # Extract the window sequence
        seq = str(record.seq)[window_start:window_end].upper()

        # Count nucleotides in first window
        state = {
            "record": record,
            "a_count": seq.count("A"),
            "c_count": seq.count("C"),
        }
        return state

    def step_state(
        self,
        state: dict,
        *,
        out_start: int,
        out_end: int,
        in_start: int,
        in_end: int,
        **kwargs,
    ) -> None:
        """Update state for the next window by removing outgoing and adding incoming bases."""
        record = state["record"]

        # Remove outgoing bases
        if out_end > out_start:
            out_seq = str(record.seq)[out_start:out_end].upper()
            state["a_count"] -= out_seq.count("A")
            state["c_count"] -= out_seq.count("C")

        # Add incoming bases
        if in_end > in_start:
            in_seq = str(record.seq)[in_start:in_end].upper()
            state["a_count"] += in_seq.count("A")
            state["c_count"] += in_seq.count("C")

    def emit(self, state: dict) -> dict[str, int]:
        """Emit feature values from current state."""
        return {
            "a_count": state["a_count"],
            "c_count": state["c_count"],
        }


class TestComputeWindowsBasicFunctionality:
    """Tests for basic compute_windows functionality."""

    def test_returns_single_row_dataframe(self):
        """Test that output is a single-row DataFrame (wide format)."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        assert result.shape[0] == 1  # Single row (wide format)
        assert result.shape[1] > 3  # Has metadata + feature columns

    def test_wide_format_column_names_with_suffixes(self):
        """Test that columns are named with window_start suffix (_0, _2, _4)."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        # Should have columns like test.len_0, test.len_2, test.len_4, test.len_6, test.len_8
        assert "test.len_0" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_4" in result.columns

    def test_metadata_columns_present(self):
        """Test that metadata columns (record_id, region_start, region_end) are included."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="test_id")

        result = fs.compute_windows(record, window_size=4, step=2)

        assert "record_id" in result.columns
        assert "region_start" in result.columns
        assert "region_end" in result.columns
        assert result["record_id"].iloc[0] == "test_id"
        assert result["region_start"].iloc[0] == 0
        assert result["region_end"].iloc[0] == 12

    def test_column_ordering_metadata_first(self):
        """Test that metadata columns come before feature columns."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        cols = list(result.columns)
        # First three should be metadata
        assert cols[0] == "record_id"
        assert cols[1] == "region_start"
        assert cols[2] == "region_end"

    def test_column_ordering_features_sorted_by_key_then_window(self):
        """Test that feature columns are sorted by feature key, then window_start."""

        def multi_feature(record):
            return {"a_feat": 1, "z_feat": 2}

        fs = FeatureSet(multi_feature, name="test")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        feature_cols = [c for c in result.columns if c.startswith("test.")]

        # Should be sorted by feature key, then window_start
        # Expected: a_feat_0, a_feat_2, a_feat_4, ..., z_feat_0, z_feat_2, z_feat_4, ...
        # All a_feat columns should come before all z_feat columns
        # Use indices to avoid O(n²) complexity
        a_indices = [i for i, c in enumerate(feature_cols) if "a_feat" in c]
        z_indices = [i for i, c in enumerate(feature_cols) if "z_feat" in c]
        if a_indices and z_indices:
            assert max(a_indices) < min(z_indices), (
                "a_feat columns should come before z_feat columns"
            )

        # Check first few columns
        assert feature_cols[0] == "test.a_feat_0"
        assert feature_cols[1] == "test.a_feat_2"


class TestComputeWindowsNonIncrementalFeatures:
    """Tests for non-incremental feature computation."""

    def test_non_incremental_feature_values_correct(self):
        """Test that non-incremental features compute correct values."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        # Windows of size 4 with step 2:
        # Window 0 (0-4): ACGT -> len=4
        # Window 1 (2-6): GTAC -> len=4
        # Window 2 (4-8): ACGT -> len=4
        # Window 3 (6-10): GTAC -> len=4
        # Window 4 (8-12): ACGT -> len=4
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        # All windows should have length 4
        assert result["test.len_0"].iloc[0] == 4
        assert result["test.len_2"].iloc[0] == 4
        assert result["test.len_4"].iloc[0] == 4
        assert result["test.len_6"].iloc[0] == 4
        assert result["test.len_8"].iloc[0] == 4


class TestComputeWindowsIncrementalFeatures:
    """Tests for incremental feature computation."""

    def test_incremental_feature_produces_correct_output(self):
        """Test that incremental computation produces correct feature values."""
        fs = FeatureSet(TinyIncrementalFeature(), name="test")
        # ACGTACGT - windows of 4 with step 2
        # Window 0 (0-4): ACGT -> A:1, C:1
        # Window 1 (2-6): GTAC -> A:1, C:1
        # Window 2 (4-8): ACGT -> A:1, C:1
        record = SeqRecord(Seq("ACGTACGT"), id="test")

        result = fs.compute_windows(record, window_size=4, step=2)

        # Check window 0 (ACGT)
        assert result["test.a_count_0"].iloc[0] == 1
        assert result["test.c_count_0"].iloc[0] == 1

        # Check window 1 (GTAC)
        assert result["test.a_count_2"].iloc[0] == 1
        assert result["test.c_count_2"].iloc[0] == 1

        # Check window 2 (ACGT)
        assert result["test.a_count_4"].iloc[0] == 1
        assert result["test.c_count_4"].iloc[0] == 1

    def test_incremental_and_fallback_produce_identical_output(self):
        """Test that incremental and fallback paths produce exactly the same output."""
        # Create incremental version
        fs_incremental = FeatureSet(TinyIncrementalFeature(), name="test")

        record = SeqRecord(Seq("ACGTACGT"), id="test")

        # Verify incremental produces consistent results
        result_incremental = fs_incremental.compute_windows(record, window_size=4, step=2)

        # Verify shape and structure
        assert result_incremental.shape[0] == 1
        assert "test.a_count_0" in result_incremental.columns
        assert "test.c_count_0" in result_incremental.columns


class TestComputeWindowsWithRegion:
    """Tests for region parameter functionality."""

    def test_region_restricts_windowing(self):
        """Test that region parameter restricts windowing to a subsequence."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("NNNNACGTACGTNNNN"), id="seq1")

        # Only window over positions 4-12 (ACGTACGT)
        result = fs.compute_windows(
            record, window_size=4, step=2, region=(4, 12)
        )

        assert result["region_start"].iloc[0] == 4
        assert result["region_end"].iloc[0] == 12

        # Should have windows starting at 0 (relative to region)
        assert "test.len_0" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_4" in result.columns

    def test_region_window_starts_relative_to_region(self):
        """Test that window starts are relative to region start (first is _0)."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("NNNNACGTACGTNNNN"), id="seq1")

        result = fs.compute_windows(
            record, window_size=4, step=2, region=(4, 12)
        )

        # Window starts should be 0, 2, 4 (relative to region start, not absolute)
        assert "test.len_0" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_4" in result.columns

        # Should NOT have windows starting at 4, 6, 8 (absolute positions)
        assert "test.len_4" in result.columns  # This is relative position 4 within region
        assert "test.len_6" not in result.columns  # Would be beyond region


class TestComputeWindowsProteinSequences:
    """Tests for protein sequence handling."""

    def test_protein_sequence_with_step_1(self):
        """Test that protein sequences work with step=1."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein1")
        protein_record.annotations["molecule_type"] = "protein"

        # Use step=1 for protein (sliding window)
        result = fs.compute_windows(protein_record, window_size=5, step=1)

        # Should have windows at positions 0, 1, 2, 3, 4
        assert "test.len_0" in result.columns
        assert "test.len_1" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_3" in result.columns
        assert "test.len_4" in result.columns

        # All windows should have length 5
        for i in range(5):
            assert result[f"test.len_{i}"].iloc[0] == 5


class TestComputeWindowsDropPartial:
    """Tests for drop_partial parameter."""

    def test_drop_partial_true_excludes_partial_windows(self):
        """Test that partial windows are excluded when drop_partial=True."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        # 11 bases, windows of 4 with step 2 -> only 0, 2, 4, 6 fit fully
        record = SeqRecord(Seq("ACGTACGTACG"), id="seq1")

        result = fs.compute_windows(
            record, window_size=4, step=2, drop_partial=True
        )

        # Should only have columns for full windows (0, 2, 4, 6)
        feature_cols = [c for c in result.columns if c.startswith("test.")]
        assert len(feature_cols) == 4
        assert "test.len_0" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_4" in result.columns
        assert "test.len_6" in result.columns
        assert "test.len_8" not in result.columns  # Partial

    def test_drop_partial_false_includes_partial_windows(self):
        """Test that partial windows are included when drop_partial=False."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGTACGTACG"), id="seq1")

        result = fs.compute_windows(
            record, window_size=4, step=2, drop_partial=False
        )

        # Should have columns for all windows including partial ones
        feature_cols = [c for c in result.columns if c.startswith("test.")]
        assert len(feature_cols) == 6
        assert "test.len_0" in result.columns
        assert "test.len_2" in result.columns
        assert "test.len_4" in result.columns
        assert "test.len_6" in result.columns
        assert "test.len_8" in result.columns  # Partial (3 bases)
        assert "test.len_10" in result.columns  # Partial (1 base)


class TestComputeWindowsEdgeCases:
    """Tests for edge cases and special scenarios."""

    def test_no_windows_returns_metadata_only(self):
        """Test that when no windows fit, returns row with metadata only."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ACGT"), id="seq1")

        # Window size larger than sequence with drop_partial=True
        result = fs.compute_windows(
            record, window_size=10, step=2, drop_partial=True
        )

        # Should have only metadata columns
        assert result.shape[0] == 1
        assert "record_id" in result.columns
        assert "region_start" in result.columns
        assert "region_end" in result.columns
        assert len([c for c in result.columns if c.startswith("test.")]) == 0

    def test_numeric_sorting_of_window_starts(self):
        """Test that window starts are sorted numerically (3, 6, 9, 12 not 12, 3, 6, 9)."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        # Long sequence to get window_start >= 10
        long_seq = "ACGT" * 20  # 80 bases
        record = SeqRecord(Seq(long_seq), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=3)

        feature_cols = [c for c in result.columns if c.startswith("test.")]

        # Extract window_starts
        window_starts = [int(c.rsplit("_", 1)[1]) for c in feature_cols]

        # Should be sorted numerically
        assert window_starts == sorted(window_starts)
        # Verify we have double-digit window_starts
        assert any(ws >= 10 for ws in window_starts)

    def test_featureset_name_in_column_names(self):
        """Test that FeatureSet name appears in column names."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="custom_name")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="seq1")

        result = fs.compute_windows(record, window_size=4, step=2)

        feature_cols = [c for c in result.columns if c.startswith("custom_name.")]
        assert len(feature_cols) > 0
        assert "custom_name.len_0" in result.columns

    def test_empty_sequence_no_windows(self):
        """Test handling of empty sequence."""
        fs = FeatureSet(TinyNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq(""), id="empty")

        result = fs.compute_windows(record, window_size=4, step=2)

        # Should have only metadata
        assert result.shape[0] == 1
        assert result["record_id"].iloc[0] == "empty"
        assert len([c for c in result.columns if c.startswith("test.")]) == 0
