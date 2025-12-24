"""Tests for PositionalFeature protocol and aggregation."""

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import (
    AggregationSpec,
    PositionSpace,
    geometric_mean,
)


class ToyResidueFeature:
    """A toy positional feature that computes GC content per residue."""

    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space."""
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specs for each feature key."""
        return {
            "gc": AggregationSpec(aggregation_fn=np.mean),
            "gc_geomean": AggregationSpec(aggregation_fn=geometric_mean),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-residue GC indicator (1.0 for G/C, 0.0 otherwise)."""
        seq = str(record.seq).upper()
        gc_vector = np.array([1.0 if b in "GC" else 0.0 for b in seq])
        return {
            "gc": gc_vector,
            "gc_geomean": gc_vector + 0.1,  # Add small offset to test geomean
        }


class ToyCodonFeature:
    """A toy positional feature that computes features per codon."""

    @property
    def position_space(self) -> PositionSpace:
        """Return CODON position space."""
        return PositionSpace.CODON

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specs for each feature key."""
        return {
            "start_with_a": AggregationSpec(aggregation_fn=np.mean),
            "codon_sum": AggregationSpec(aggregation_fn=np.sum),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-codon features."""
        seq = str(record.seq).upper()
        # Process sequence in codons (groups of 3)
        num_codons = len(seq) // 3
        start_with_a = np.array([1.0 if seq[i * 3] == "A" else 0.0 for i in range(num_codons)])
        codon_values = np.array([float(i + 1) for i in range(num_codons)])

        return {
            "start_with_a": start_with_a,
            "codon_sum": codon_values,
        }


class TestPositionalFeature:
    """Tests for PositionalFeature protocol implementation."""

    def test_toy_residue_feature_has_required_properties(self):
        """Test that toy residue feature has required protocol properties."""
        feature = ToyResidueFeature()

        # Check position_space property
        assert hasattr(feature, "position_space")
        assert feature.position_space == PositionSpace.RESIDUE

        # Check vector_keys property
        assert hasattr(feature, "vector_keys")
        assert isinstance(feature.vector_keys, dict)
        assert "gc" in feature.vector_keys
        assert isinstance(feature.vector_keys["gc"], AggregationSpec)

        # Check compute_vector method
        assert hasattr(feature, "compute_vector")
        assert callable(feature.compute_vector)

    def test_toy_codon_feature_has_required_properties(self):
        """Test that toy codon feature has required protocol properties."""
        feature = ToyCodonFeature()

        assert feature.position_space == PositionSpace.CODON
        assert isinstance(feature.vector_keys, dict)
        assert "start_with_a" in feature.vector_keys
        assert callable(feature.compute_vector)

    def test_compute_vector_residue_feature(self):
        """Test compute_vector returns correct per-residue values."""
        feature = ToyResidueFeature()
        record = SeqRecord(Seq("ACGTGCTA"), id="test")

        result = feature.compute_vector(record)

        # Check structure
        assert isinstance(result, dict)
        assert "gc" in result
        assert "gc_geomean" in result

        # Check values - GC at positions 2,3,4,5 out of 8
        gc_vector = result["gc"]
        assert isinstance(gc_vector, np.ndarray)
        assert len(gc_vector) == 8
        expected = np.array([0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        np.testing.assert_array_equal(gc_vector, expected)

    def test_compute_vector_codon_feature(self):
        """Test compute_vector returns correct per-codon values."""
        feature = ToyCodonFeature()
        # Sequence with 3 codons: ATG, GCA, TAA
        record = SeqRecord(Seq("ATGGCATAA"), id="test")

        result = feature.compute_vector(record)

        # Check structure
        assert isinstance(result, dict)
        assert "start_with_a" in result
        assert "codon_sum" in result

        # Check values - only first codon starts with A
        start_vector = result["start_with_a"]
        assert len(start_vector) == 3
        expected_start = np.array([1.0, 0.0, 0.0])
        np.testing.assert_array_equal(start_vector, expected_start)

        # Check codon_sum values
        codon_sum_vector = result["codon_sum"]
        assert len(codon_sum_vector) == 3
        expected_sum = np.array([1.0, 2.0, 3.0])
        np.testing.assert_array_equal(codon_sum_vector, expected_sum)


class TestAggregationHelper:
    """Tests for aggregation functions."""

    def test_mean_aggregation(self):
        """Test mean aggregation using numpy mean."""
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        spec = AggregationSpec(aggregation_fn=np.mean)

        result = spec.aggregation_fn(values)

        assert result == 3.0

    def test_sum_aggregation(self):
        """Test sum aggregation using numpy sum."""
        values = np.array([1.0, 2.0, 3.0])
        spec = AggregationSpec(aggregation_fn=np.sum)

        result = spec.aggregation_fn(values)

        assert result == 6.0

    def test_max_aggregation(self):
        """Test max aggregation using numpy max."""
        values = np.array([1.0, 5.0, 3.0, 2.0])
        spec = AggregationSpec(aggregation_fn=np.max)

        result = spec.aggregation_fn(values)

        assert result == 5.0

    def test_geometric_mean_positive_values(self):
        """Test geometric mean with positive values."""
        # Geometric mean of [1, 2, 4, 8] = (1*2*4*8)^(1/4) = 32^0.25 = 2.378...
        values = np.array([1.0, 2.0, 4.0, 8.0])

        result = geometric_mean(values)

        # Geometric mean should match exp(mean(log(values)))
        expected = np.exp(np.mean(np.log(values)))
        assert abs(result - expected) < 1e-10

    def test_geometric_mean_with_zero(self):
        """Test geometric mean returns 0.0 when any value is zero."""
        values = np.array([1.0, 2.0, 0.0, 4.0])

        result = geometric_mean(values)

        assert result == 0.0

    def test_geometric_mean_with_negative(self):
        """Test geometric mean returns 0.0 when any value is negative."""
        values = np.array([1.0, 2.0, -1.0, 4.0])

        result = geometric_mean(values)

        assert result == 0.0

    def test_geometric_mean_empty_array(self):
        """Test geometric mean returns 0.0 for empty array."""
        values = np.array([])

        result = geometric_mean(values)

        assert result == 0.0

    def test_geometric_mean_all_ones(self):
        """Test geometric mean of all ones is 1.0."""
        values = np.array([1.0, 1.0, 1.0, 1.0])

        result = geometric_mean(values)

        assert abs(result - 1.0) < 1e-10


class TestWindowAggregation:
    """Tests for aggregating per-position values into windows."""

    def test_aggregate_residue_window(self):
        """Test aggregating residue features over a window."""
        feature = ToyResidueFeature()
        record = SeqRecord(Seq("ACGTGCTA"), id="test")

        # Compute per-position values
        vectors = feature.compute_vector(record)

        # Simulate windowing: aggregate first 4 positions
        window_start = 0
        window_end = 4
        window_values = vectors["gc"][window_start:window_end]

        # Apply mean aggregation
        spec = feature.vector_keys["gc"]
        aggregated = spec.aggregation_fn(window_values)

        # Expected: [0,1,1,0] -> mean = 0.5
        assert aggregated == 0.5
        assert aggregated == 0.5

    def test_aggregate_codon_window(self):
        """Test aggregating codon features over a window."""
        feature = ToyCodonFeature()
        # 4 codons: ATG, GCA, TAA, CCC
        record = SeqRecord(Seq("ATGGCATAACCC"), id="test")

        # Compute per-codon values
        vectors = feature.compute_vector(record)

        # Aggregate first 2 codons for codon_sum
        window_start = 0
        window_end = 2
        window_values = vectors["codon_sum"][window_start:window_end]

        # Apply sum aggregation
        spec = feature.vector_keys["codon_sum"]
        aggregated = spec.aggregation_fn(window_values)

        # Expected: [1.0, 2.0] -> sum = 3.0
        assert aggregated == 3.0

    def test_aggregate_multiple_windows(self):
        """Test aggregating multiple overlapping windows."""
        feature = ToyResidueFeature()
        record = SeqRecord(Seq("GCGCATAT"), id="test")

        # Compute per-position values
        vectors = feature.compute_vector(record)
        gc_vector = vectors["gc"]

        # GC pattern: [1,1,1,1,0,0,0,0]

        # Window 1: positions 0-4 -> [1,1,1,1] -> mean = 1.0
        window1 = feature.vector_keys["gc"].aggregation_fn(gc_vector[0:4])
        assert window1 == 1.0

        # Window 2: positions 2-6 -> [1,1,0,0] -> mean = 0.5
        window2 = feature.vector_keys["gc"].aggregation_fn(gc_vector[2:6])
        assert window2 == 0.5

        # Window 3: positions 4-8 -> [0,0,0,0] -> mean = 0.0
        window3 = feature.vector_keys["gc"].aggregation_fn(gc_vector[4:8])
        assert window3 == 0.0


class TestPositionSpace:
    """Tests for PositionSpace enum."""

    def test_position_space_values(self):
        """Test PositionSpace enum has expected values."""
        assert PositionSpace.RESIDUE.value == "residue"
        assert PositionSpace.CODON.value == "codon"

    def test_position_space_comparison(self):
        """Test PositionSpace enum comparison."""
        feature1 = ToyResidueFeature()
        feature2 = ToyCodonFeature()

        assert feature1.position_space == PositionSpace.RESIDUE
        assert feature2.position_space == PositionSpace.CODON
        assert feature1.position_space != feature2.position_space
