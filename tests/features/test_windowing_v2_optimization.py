"""Tests for windowing v2 union-of-windows optimization and stride semantics."""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class CallCountingCodonFeature:
    """A test feature that counts how many times compute_vector is called."""

    def __init__(self):
        """Initialize with call counter."""
        self.compute_vector_call_count = 0
        self.last_positions_arg = None

    @property
    def position_space(self) -> PositionSpace:
        """Return CODON position space."""
        return PositionSpace.CODON

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specs for each feature key."""
        return {
            "codon_idx": AggregationSpec(aggregation_fn=np.mean),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-codon indices, tracking calls and positions parameter."""
        self.compute_vector_call_count += 1
        self.last_positions_arg = kwargs.get("positions", None)

        seq = str(record.seq).upper()
        num_codons = len(seq) // 3

        # Support positions parameter for optimization
        positions = kwargs.get("positions", None)
        if positions is not None:
            # Return sparse vector for requested positions only
            codon_idx = np.array(positions, dtype=float)
        else:
            # Return full vector
            codon_idx = np.arange(num_codons, dtype=float)

        return {"codon_idx": codon_idx}


class TestUnionOfWindowsOptimization:
    """Tests for union-of-windows position computation optimization."""

    def test_compute_vector_called_once_per_feature(self):
        """Test that compute_vector is called exactly once per feature per ORF."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 30 nt = 10 codons, create 3 windows
        record = SeqRecord(Seq("ATG" * 10), id="test")

        # Call compute_orf_windows_v2 with multiple windows
        fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=12, step_nt=6)

        # Verify compute_vector was called exactly once
        assert feature.compute_vector_call_count == 1

    def test_positions_parameter_passed_to_compute_vector(self):
        """Test that union positions are passed to compute_vector."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 30 nt = 10 codons
        record = SeqRecord(Seq("ATG" * 10), id="test")

        # Call with windows
        fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=12, step_nt=6)

        # Verify positions parameter was passed
        assert feature.last_positions_arg is not None
        assert isinstance(feature.last_positions_arg, np.ndarray)

    def test_union_positions_are_sorted_unique(self):
        """Test that union positions are sorted and unique."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 30 nt = 10 codons
        record = SeqRecord(Seq("ATG" * 10), id="test")

        # Windows: 0-12 (codons 0-4), 6-18 (codons 2-6), 12-24 (codons 4-8), 18-30 (codons 6-10)
        fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=12, step_nt=6)

        # Get the positions that were passed
        positions = feature.last_positions_arg

        # Verify sorted
        assert np.all(positions[:-1] <= positions[1:])

        # Verify unique
        assert len(positions) == len(np.unique(positions))

    def test_union_includes_all_window_positions(self):
        """Test that union includes all positions used by any window."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 30 nt = 10 codons
        record = SeqRecord(Seq("ATG" * 10), id="test")

        # Windows with step=6: 0-12 (codons 0-4), 6-18 (2-6), 12-24 (4-8), 18-30 (6-10)
        # Union should be all codons 0-9
        fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=12, step_nt=6)

        positions = feature.last_positions_arg
        expected_union = set(range(10))  # All 10 codons

        assert set(positions) == expected_union


class TestSparseWindowOptimization:
    """Tests for optimization with sparse (non-overlapping) windows."""

    def test_sparse_windows_compute_subset_of_positions(self):
        """Test that sparse windows only compute needed positions."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 100 codons (300 nt), large step creates sparse windows
        record = SeqRecord(Seq("ATG" * 100), id="test")

        # Windows: 0-12 (codons 0-4), 30-42 (10-14), 60-72 (20-24), 90-102 (30-34), etc.
        fs.compute_orf_windows_v2(record, orf=(0, 300), window_nt=12, step_nt=30)

        positions = feature.last_positions_arg

        # Should have fewer positions than full sequence
        assert len(positions) < 100

        # With step=30 (10 codons) and window=12 (4 codons), we get 10 windows
        # covering positions 0-4, 10-14, 20-24, 30-34, 40-44, 50-54, 60-64, 70-74, 80-84, 90-94
        # Total: 40 positions out of 100
        assert len(positions) == 40

    def test_overlapping_windows_union_is_larger(self):
        """Test that overlapping windows produce larger union."""
        feature1 = CallCountingCodonFeature()
        feature2 = CallCountingCodonFeature()

        record = SeqRecord(Seq("ATG" * 100), id="test")

        # Sparse windows
        fs1 = FeatureSet(feature1, name="sparse")
        fs1.compute_orf_windows_v2(record, orf=(0, 300), window_nt=12, step_nt=30)
        sparse_union_size = len(feature1.last_positions_arg)

        # Dense windows
        fs2 = FeatureSet(feature2, name="dense")
        fs2.compute_orf_windows_v2(record, orf=(0, 300), window_nt=12, step_nt=3)
        dense_union_size = len(feature2.last_positions_arg)

        # Dense should have more positions in union
        assert dense_union_size > sparse_union_size


class TestStrideSemantics:
    """Tests for correct stride sampling semantics."""

    def test_step_nt_6_means_window_stride_not_codon_stride(self):
        """Test that step_nt=6 moves windows by 2 codons, doesn't subsample codons."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        # 18 nt = 6 codons
        record = SeqRecord(Seq("ATG" * 6), id="test")

        # step_nt=6 means windows are placed every 2 codons
        # Window 0: codons 0-3 (nt 0-12)
        # Window 6: codons 2-5 (nt 6-18)
        result = fs.compute_orf_windows_v2(record, orf=(0, 18), window_nt=12, step_nt=6)

        # Each window should aggregate ALL codons in its range, not subsample
        # Window 0: codons 0,1,2,3 -> mean = (0+1+2+3)/4 = 1.5
        assert abs(result["test.codon_idx_0"].iloc[0] - 1.5) < 1e-10

        # Window 6: codons 2,3,4,5 -> mean = (2+3+4+5)/4 = 3.5
        assert abs(result["test.codon_idx_6"].iloc[0] - 3.5) < 1e-10

    def test_codon_space_step_must_be_multiple_of_3(self):
        """Test that step_nt must be multiple of 3 for codon-space features."""
        feature = CallCountingCodonFeature()
        fs = FeatureSet(feature, name="test")

        record = SeqRecord(Seq("ATG" * 5), id="test")

        # step_nt=4 is not a multiple of 3
        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=9, step_nt=4)


class TestBackwardCompatibility:
    """Tests ensuring backward compatibility with features that don't use positions."""

    def test_feature_without_positions_support_still_works(self):
        """Test that features ignoring positions parameter still work correctly."""

        class LegacyCodonFeature:
            """Feature that doesn't use positions parameter."""

            @property
            def position_space(self):
                return PositionSpace.CODON

            @property
            def vector_keys(self):
                return {"idx": AggregationSpec(aggregation_fn=np.mean)}

            def compute_vector(self, record, **kwargs):
                # Ignore positions parameter, always compute full vector
                seq = str(record.seq).upper()
                num_codons = len(seq) // 3
                return {"idx": np.arange(num_codons, dtype=float)}

        feature = LegacyCodonFeature()
        fs = FeatureSet(feature, name="legacy")

        record = SeqRecord(Seq("ATG" * 10), id="test")
        result = fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=12, step_nt=6)

        # Should work correctly despite not using positions optimization
        assert "legacy.idx_0" in result.columns
        assert abs(result["legacy.idx_0"].iloc[0] - 1.5) < 1e-10
