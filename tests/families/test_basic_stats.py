"""Tests for basic_stats family feature."""

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.windowing import iter_windows
from biotooler.families.basic_stats import BasicStatsFeature
from biotooler.features import FeatureSet


class TestBasicStatsBaseline:
    """Tests for BasicStatsFeature baseline computation."""

    def test_dna_sequence_basic(self):
        """Test basic DNA sequence statistics."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGT"), id="test")
        result = feature(record)

        assert result["length"] == 8
        assert result["count_a"] == 2
        assert result["count_c"] == 2
        assert result["count_g"] == 2
        assert result["count_t"] == 2
        assert result["fraction_a"] == 0.25
        assert result["fraction_c"] == 0.25
        assert result["fraction_g"] == 0.25
        assert result["fraction_t"] == 0.25
        assert result["gc_fraction"] == 0.5

    def test_rna_sequence_basic(self):
        """Test basic RNA sequence statistics (U instead of T)."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGUACGU"), id="test")
        result = feature(record)

        assert result["length"] == 8
        assert result["count_a"] == 2
        assert result["count_c"] == 2
        assert result["count_g"] == 2
        assert result["count_u"] == 2
        assert "count_t" not in result  # RNA should not have T
        assert result["gc_fraction"] == 0.5

    def test_protein_sequence_basic(self):
        """Test basic protein sequence statistics."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        result = feature(record)

        assert result["length"] == 9
        assert result["count_m"] == 1
        assert result["count_k"] == 1
        assert result["count_a"] == 1
        assert "gc_fraction" not in result  # Protein should not have GC fraction

    def test_dna_with_n(self):
        """Test DNA sequence with N (ambiguous base)."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTNACGTN"), id="test")
        result = feature(record)

        assert result["length"] == 10
        assert result["count_a"] == 2
        assert result["count_c"] == 2
        assert result["count_g"] == 2
        assert result["count_t"] == 2
        assert result["count_n"] == 2
        # GC fraction is (G+C)/total_length, N is counted in length but not as G or C
        assert result["gc_fraction"] == 0.4  # (2+2)/10

    def test_high_gc_content(self):
        """Test sequence with high GC content."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("GGGCCCGGGCCC"), id="test")
        result = feature(record)

        assert result["gc_fraction"] == 1.0
        assert result["count_a"] == 0
        assert result["count_t"] == 0

    def test_low_gc_content(self):
        """Test sequence with low GC content."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("AAAATTTTAAAA"), id="test")
        result = feature(record)

        assert result["gc_fraction"] == 0.0
        assert result["count_g"] == 0
        assert result["count_c"] == 0

    def test_empty_sequence(self):
        """Test empty sequence edge case."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq(""), id="test")
        result = feature(record)

        assert result["length"] == 0
        # Empty sequence should have 0 fractions, not divide by zero
        assert result["fraction_a"] == 0.0
        assert result["fraction_c"] == 0.0
        assert result["gc_fraction"] == 0.0

    def test_single_nucleotide(self):
        """Test single nucleotide sequence."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("A"), id="test")
        result = feature(record)

        assert result["length"] == 1
        assert result["count_a"] == 1
        assert result["fraction_a"] == 1.0
        assert result["fraction_c"] == 0.0
        assert result["gc_fraction"] == 0.0


class TestBasicStatsIncremental:
    """Tests for BasicStatsFeature incremental computation."""

    def test_init_state_creates_correct_structure(self):
        """Test that init_state creates proper state structure."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGT"), id="test")

        state = feature.init_state(record, window_start=0, window_end=4)

        assert "record" in state
        assert "counts" in state
        assert "alphabet" in state
        assert "char_to_idx" in state
        assert isinstance(state["counts"], np.ndarray)
        assert state["counts"].dtype == np.int64

    def test_step_state_updates_counts(self):
        """Test that step_state correctly updates counts."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGT"), id="test")

        # Initialize with first window [0:4] = "ACGT"
        state = feature.init_state(record, window_start=0, window_end=4)
        initial_emit = feature.emit(state)
        assert initial_emit["count_a"] == 1
        assert initial_emit["count_c"] == 1
        assert initial_emit["count_g"] == 1
        assert initial_emit["count_t"] == 1

        # Step to window [2:6] = "GTAC"
        # Remove [0:2] = "AC", Add [4:6] = "AC"
        feature.step_state(state, out_start=0, out_end=2, in_start=4, in_end=6)
        stepped_emit = feature.emit(state)

        # After step: still have A, C, G, T each once
        assert stepped_emit["count_a"] == 1
        assert stepped_emit["count_c"] == 1
        assert stepped_emit["count_g"] == 1
        assert stepped_emit["count_t"] == 1

    def test_emit_produces_same_format_as_baseline(self):
        """Test that emit produces same format as __call__."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGT"), id="test")

        # Baseline
        baseline = feature(record)

        # Incremental
        state = feature.init_state(record, window_start=0, window_end=8)
        incremental = feature.emit(state)

        # Should have same keys
        assert set(baseline.keys()) == set(incremental.keys())

        # Should have same values
        for key in baseline.keys():
            assert baseline[key] == incremental[key]


class TestBasicStatsIncrementalVsBaseline:
    """Tests comparing incremental and baseline computation."""

    def test_single_window_exact_match(self):
        """Test that incremental matches baseline for a single window."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGTGGGCCC"), id="test")

        # Baseline
        baseline = feature(record)

        # Incremental
        state = feature.init_state(record, window_start=0, window_end=14)
        incremental = feature.emit(state)

        # Exact match
        for key in baseline.keys():
            baseline_value = baseline[key]
            incremental_value = incremental[key]
            if isinstance(baseline_value, float):
                assert abs(baseline_value - incremental_value) < 1e-10  # type: ignore[operator]
            else:
                assert baseline_value == incremental_value

    def test_sliding_windows_match_baseline(self):
        """Test that incremental matches baseline for sliding windows."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGTACGTACGT"), id="test")
        window_size = 8
        step = 4

        # Baseline: compute each window independently
        baseline_results = []
        for window in iter_windows(record, window_size=window_size, step=step):
            baseline_results.append(feature(window))

        # Incremental: use manual state stepping
        incremental_results = []
        state = None
        for _i, window in enumerate(iter_windows(record, window_size=window_size, step=step)):
            window_start = int(window.annotations["start"])
            window_end = int(window.annotations["end"])

            if state is None:
                state = feature.init_state(record, window_start=window_start, window_end=window_end)
            else:
                prev_start = window_start - step
                prev_end = prev_start + window_size
                feature.step_state(
                    state,
                    out_start=prev_start,
                    out_end=window_start,
                    in_start=prev_end,
                    in_end=window_end
                )

            incremental_results.append(feature.emit(state))

        # Compare all windows
        assert len(baseline_results) == len(incremental_results)
        for i, (baseline, incremental) in enumerate(
            zip(baseline_results, incremental_results, strict=True)
        ):
            for key in baseline.keys():
                baseline_value = baseline[key]
                incremental_value = incremental[key]
                if isinstance(baseline_value, float):
                    assert abs(baseline_value - incremental_value) < 1e-10, (  # type: ignore[operator]
                        f"Window {i}, key {key}: "
                        f"baseline={baseline_value}, incremental={incremental_value}"
                    )
                else:
                    assert baseline_value == incremental_value, (
                        f"Window {i}, key {key}: "
                        f"baseline={baseline_value}, incremental={incremental_value}"
                    )

    def test_featureset_compute_windows_matches_baseline(self):
        """Test that FeatureSet.compute_windows matches baseline computation."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGTACGTACGT"), id="test")
        window_size = 8
        step = 4

        # Baseline: compute each window independently
        baseline_results = []
        for window in iter_windows(record, window_size=window_size, step=step):
            baseline_results.append(feature(window))

        # Incremental: use FeatureSet.compute_windows
        fs = FeatureSet(feature, name="stats")
        incremental_df = fs.compute_windows(
            record,
            window_size=window_size,
            step=step
        )

        # Extract window starts from dataframe columns
        window_starts = []
        for col in incremental_df.columns:
            if col.startswith("stats.length_"):
                window_start = int(col.split("_")[-1])
                window_starts.append(window_start)
        window_starts.sort()

        # Compare each window
        assert len(baseline_results) == len(window_starts)
        for i, (baseline, window_start) in enumerate(
            zip(baseline_results, window_starts, strict=True)
        ):
            for key, baseline_value in baseline.items():
                col_name = f"stats.{key}_{window_start}"
                incremental_value = incremental_df[col_name].values[0]

                if isinstance(baseline_value, float):
                    assert abs(baseline_value - incremental_value) < 1e-10, (
                        f"Window {i} (start={window_start}), key {key}: "
                        f"baseline={baseline_value}, incremental={incremental_value}"
                    )
                else:
                    assert baseline_value == incremental_value, (
                        f"Window {i} (start={window_start}), key {key}: "
                        f"baseline={baseline_value}, incremental={incremental_value}"
                    )

    def test_overlapping_windows(self):
        """Test with highly overlapping windows (step=1)."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGTACGTAC"), id="test")
        window_size = 6
        step = 1

        # Baseline
        baseline_results = []
        for window in iter_windows(record, window_size=window_size, step=step):
            baseline_results.append(feature(window))

        # Incremental via FeatureSet
        fs = FeatureSet(feature, name="stats")
        incremental_df = fs.compute_windows(record, window_size=window_size, step=step)

        # Should have 5 windows (0-5, 1-6, 2-7, 3-8, 4-9)
        assert len(baseline_results) == 5

        # Verify each window
        for i, baseline in enumerate(baseline_results):
            for key, baseline_value in baseline.items():
                col_name = f"stats.{key}_{i}"
                incremental_value = incremental_df[col_name].values[0]
                if isinstance(baseline_value, float):
                    assert abs(baseline_value - incremental_value) < 1e-10
                else:
                    assert baseline_value == incremental_value

    def test_rna_incremental_vs_baseline(self):
        """Test RNA sequence incremental vs baseline."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ACGUACGUACGU"), id="test")
        window_size = 6
        step = 3

        # Baseline
        baseline_results = []
        for window in iter_windows(record, window_size=window_size, step=step):
            baseline_results.append(feature(window))

        # Incremental
        fs = FeatureSet(feature, name="stats")
        incremental_df = fs.compute_windows(record, window_size=window_size, step=step)

        # Compare
        for i, baseline in enumerate(baseline_results):
            window_start = i * step
            for key, baseline_value in baseline.items():
                col_name = f"stats.{key}_{window_start}"
                incremental_value = incremental_df[col_name].values[0]
                if isinstance(baseline_value, float):
                    assert abs(baseline_value - incremental_value) < 1e-10
                else:
                    assert baseline_value == incremental_value

    def test_protein_incremental_vs_baseline(self):
        """Test protein sequence incremental vs baseline."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("MKALVSWGRDEPQFIH"), id="test")
        window_size = 8
        step = 4

        # Baseline
        baseline_results = []
        for window in iter_windows(record, window_size=window_size, step=step):
            baseline_results.append(feature(window))

        # Incremental
        fs = FeatureSet(feature, name="stats")
        incremental_df = fs.compute_windows(record, window_size=window_size, step=step)

        # Compare
        for i, baseline in enumerate(baseline_results):
            window_start = i * step
            for key, baseline_value in baseline.items():
                col_name = f"stats.{key}_{window_start}"
                incremental_value = incremental_df[col_name].values[0]
                if isinstance(baseline_value, float):
                    assert abs(baseline_value - incremental_value) < 1e-10
                else:
                    assert baseline_value == incremental_value


class TestBasicStatsEdgeCases:
    """Tests for edge cases and special sequences."""

    def test_uniform_sequence(self):
        """Test sequence with only one nucleotide type."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("AAAAAAA"), id="test")
        result = feature(record)

        assert result["count_a"] == 7
        assert result["count_c"] == 0
        assert result["count_g"] == 0
        assert result["count_t"] == 0
        assert result["fraction_a"] == 1.0
        assert result["gc_fraction"] == 0.0

    def test_alternating_sequence(self):
        """Test alternating sequence pattern."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("ATATATAT"), id="test")
        result = feature(record)

        assert result["count_a"] == 4
        assert result["count_t"] == 4
        assert result["count_c"] == 0
        assert result["count_g"] == 0
        assert result["gc_fraction"] == 0.0

    def test_all_gc_sequence(self):
        """Test sequence with only G and C."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("GCGCGCGC"), id="test")
        result = feature(record)

        assert result["count_g"] == 4
        assert result["count_c"] == 4
        assert result["count_a"] == 0
        assert result["count_t"] == 0
        assert result["gc_fraction"] == 1.0

    def test_sequence_with_x_and_stop(self):
        """Test protein sequence with X (unknown) and * (stop)."""
        feature = BasicStatsFeature()
        record = SeqRecord(Seq("MKALX*SWGR"), id="test")
        result = feature(record)

        assert result["length"] == 10
        assert result["count_x"] == 1
        assert result["count_*"] == 1

    def test_alphabet_consistency_across_windows(self):
        """Test that alphabet is consistent across windows even if some chars are absent."""
        feature = BasicStatsFeature()
        # Full sequence has all bases, but individual windows may not
        record = SeqRecord(Seq("AAAACCCCGGGGTTTT"), id="test")

        # Window 1: only A's
        state1 = feature.init_state(record, window_start=0, window_end=4)
        result1 = feature.emit(state1)
        assert "count_a" in result1
        assert "count_c" in result1  # Should be present even if count is 0
        assert "count_g" in result1
        assert "count_t" in result1

        # Window 2: only C's
        state2 = feature.init_state(record, window_start=4, window_end=8)
        result2 = feature.emit(state2)
        assert "count_a" in result2  # Should be present even if count is 0
        assert "count_c" in result2
        assert "count_g" in result2
        assert "count_t" in result2


class TestBasicStatsFeatureSetIntegration:
    """Integration tests with FeatureSet."""

    def test_compute_windows_returns_wide_format(self):
        """Test that compute_windows returns wide format DataFrame."""
        feature = BasicStatsFeature()
        fs = FeatureSet(feature, name="stats")
        record = SeqRecord(Seq("ACGTACGTACGT"), id="test")

        result = fs.compute_windows(record, window_size=4, step=2)

        # Should be single row (wide format)
        assert result.shape[0] == 1

        # Should have metadata columns
        assert "record_id" in result.columns
        assert "region_start" in result.columns
        assert "region_end" in result.columns

        # Should have feature columns with window suffixes
        assert "stats.length_0" in result.columns
        assert "stats.count_a_0" in result.columns
        assert "stats.gc_fraction_0" in result.columns

    def test_compute_windows_with_region(self):
        """Test compute_windows with region parameter."""
        feature = BasicStatsFeature()
        fs = FeatureSet(feature, name="stats")
        record = SeqRecord(Seq("AAAACGTACGTCCCC"), id="test")

        # Compute only on middle region [4:11] = "CGTACGT"
        result = fs.compute_windows(
            record,
            window_size=4,
            step=2,
            region=(4, 11)
        )

        assert result["region_start"].values[0] == 4
        assert result["region_end"].values[0] == 11

        # First window should be "CGTA"
        assert result["stats.count_c_0"].values[0] == 1
        assert result["stats.count_g_0"].values[0] == 1
        assert result["stats.count_t_0"].values[0] == 1
        assert result["stats.count_a_0"].values[0] == 1

    def test_compute_windows_multiple_windows(self):
        """Test that multiple windows are all present in wide format."""
        feature = BasicStatsFeature()
        fs = FeatureSet(feature, name="stats")
        record = SeqRecord(Seq("ACGTACGTACGTACGT"), id="test")

        result = fs.compute_windows(record, window_size=8, step=4)

        # Should have windows at positions 0, 4, 8
        assert "stats.length_0" in result.columns
        assert "stats.length_4" in result.columns
        assert "stats.length_8" in result.columns

        # All should have length 8
        assert result["stats.length_0"].values[0] == 8
        assert result["stats.length_4"].values[0] == 8
        assert result["stats.length_8"].values[0] == 8
