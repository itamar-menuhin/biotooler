"""Tests for FeatureSet global computation methods (compute_global and compute_orf_global_v2)."""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class SimpleNonIncrementalFeature:
    """A simple non-incremental feature for testing compute_global."""

    def __call__(self, record: SeqRecord) -> dict[str, float]:
        """Return GC content."""
        seq = str(record.seq).upper()
        if len(seq) == 0:
            return {"gc_content": 0.0}
        gc_count = seq.count("G") + seq.count("C")
        return {"gc_content": gc_count / len(seq)}


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
            "gc": AggregationSpec(aggregation_fn=np.mean),
            "count": AggregationSpec(aggregation_fn=np.sum),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
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
            "start_a": AggregationSpec(aggregation_fn=np.mean),
            "codon_count": AggregationSpec(aggregation_fn=np.sum),
        }

    def compute_vector(
        self, record: SeqRecord, **kwargs
    ) -> dict[str, np.ndarray]:
        """Compute per-codon features."""
        seq = str(record.seq).upper()
        num_codons = len(seq) // 3
        positions = kwargs.get("positions", None)

        if positions is not None:
            # Compute only for requested codon positions
            start_a = np.array(
                [1.0 if seq[i * 3] == "A" else 0.0 for i in positions]
            )
            codon_count = np.ones(len(positions))
        else:
            # Compute for all codons (backward compatibility)
            start_a = np.array(
                [1.0 if seq[i * 3] == "A" else 0.0 for i in range(num_codons)]
            )
            codon_count = np.ones(num_codons)

        return {
            "start_a": start_a,
            "codon_count": codon_count,
        }


class TestComputeGlobal:
    """Tests for compute_global method."""

    def test_compute_global_emits_global_suffix(self):
        """Test that compute_global uses _GLOBAL suffix instead of window indices."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="gc")
        record = SeqRecord(Seq("ATGCGCATGC"), id="test")

        result = fs.compute_global(record)

        # Should have _GLOBAL suffix, not _0, _3, etc.
        assert "gc.gc_content_GLOBAL" in result.columns
        # Should NOT have any window index columns
        assert "gc.gc_content_0" not in result.columns
        assert "gc.gc_content_3" not in result.columns

    def test_compute_global_returns_single_row(self):
        """Test that compute_global returns a single-row DataFrame."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ATGCGCATGC"), id="seq1")

        result = fs.compute_global(record)

        assert result.shape[0] == 1

    def test_compute_global_metadata_columns(self):
        """Test that compute_global includes correct metadata columns."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ATGCGCATGC"), id="test_id")

        result = fs.compute_global(record)

        assert "record_id" in result.columns
        assert "region_start" in result.columns
        assert "region_end" in result.columns
        assert result["record_id"].iloc[0] == "test_id"
        assert result["region_start"].iloc[0] == 0
        assert result["region_end"].iloc[0] == 10

    def test_compute_global_with_region(self):
        """Test that compute_global respects region parameter."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="gc")
        record = SeqRecord(Seq("ATGCGCATGC"), id="test")

        result = fs.compute_global(record, region=(2, 8))

        assert result["region_start"].iloc[0] == 2
        assert result["region_end"].iloc[0] == 8
        # GC content of "GCGCAT" (4 GC out of 6) = 0.667
        assert abs(result["gc.gc_content_GLOBAL"].iloc[0] - 4 / 6) < 1e-10

    def test_compute_global_column_ordering(self):
        """Test that metadata columns come before feature columns."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ATGCGCATGC"), id="seq1")

        result = fs.compute_global(record)

        cols = list(result.columns)
        # First three should be metadata
        assert cols[0] == "record_id"
        assert cols[1] == "region_start"
        assert cols[2] == "region_end"


class TestComputeOrfGlobalV2:
    """Tests for compute_orf_global_v2 method."""

    def test_compute_orf_global_v2_emits_global_suffix(self):
        """Test that compute_orf_global_v2 uses _GLOBAL suffix."""
        fs = FeatureSet(ToyResidueFeature(), name="gc_feat")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        # Should have _GLOBAL suffix
        assert "gc_feat.gc_GLOBAL" in result.columns
        assert "gc_feat.count_GLOBAL" in result.columns
        # Should NOT have any window index columns
        assert "gc_feat.gc_0" not in result.columns
        assert "gc_feat.gc_3" not in result.columns

    def test_compute_orf_global_v2_returns_single_row(self):
        """Test that compute_orf_global_v2 returns a single-row DataFrame."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        assert result.shape[0] == 1

    def test_compute_orf_global_v2_metadata_columns(self):
        """Test that compute_orf_global_v2 includes correct metadata columns."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test_id")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns
        assert result["record_id"].iloc[0] == "test_id"
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15

    def test_compute_orf_global_v2_aggregates_full_orf(self):
        """Test that compute_orf_global_v2 aggregates over the entire ORF."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        # Sequence: ATGAAACCCGGGTTT (15 nt)
        # GC content: [0,0,1,0,0,0,1,1,1,1,1,1,0,0,0]
        # Mean = 7/15 ≈ 0.4667
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        # Check GC mean aggregation over full sequence
        expected_gc = 7 / 15
        assert abs(result["test.gc_GLOBAL"].iloc[0] - expected_gc) < 1e-10

        # Check count sum aggregation (should be full ORF length)
        assert result["test.count_GLOBAL"].iloc[0] == 15.0

    def test_compute_orf_global_v2_with_codon_feature(self):
        """Test compute_orf_global_v2 with a codon-space feature."""
        fs = FeatureSet(ToyCodonFeature(), name="codon_feat")
        # Sequence: ATGAAACCCGGGTTT (15 nt = 5 codons)
        # Codons: ATG AAA CCC GGG TTT
        # Start with A: [1, 1, 0, 0, 0] -> mean = 2/5 = 0.4
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        # Check codon aggregation
        assert abs(result["codon_feat.start_a_GLOBAL"].iloc[0] - 2 / 5) < 1e-10
        assert result["codon_feat.codon_count_GLOBAL"].iloc[0] == 5.0

    def test_compute_orf_global_v2_non_positional_raises_error(self):
        """Test that non-positional features raise ValueError."""

        def simple_feature(record):
            return {"length": len(record.seq)}

        fs = FeatureSet(simple_feature, name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        with pytest.raises(
            ValueError, match="does not implement PositionalFeature protocol"
        ):
            fs.compute_orf_global_v2(record, orf=(0, 15))

    def test_compute_orf_global_v2_column_ordering(self):
        """Test that metadata columns come before feature columns."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        result = fs.compute_orf_global_v2(record, orf=(0, 15))

        cols = list(result.columns)
        # First three should be metadata
        assert cols[0] == "record_id"
        assert cols[1] == "orf_start"
        assert cols[2] == "orf_end"


class TestGlobalColumnOrderDeterministic:
    """Tests for deterministic column ordering in global methods."""

    def test_global_column_order_is_deterministic(self):
        """Test that column order in compute_global is deterministic."""
        fs = FeatureSet(SimpleNonIncrementalFeature(), name="test")
        record = SeqRecord(Seq("ATGCGCATGC"), id="seq1")

        # Run multiple times and check order is consistent
        results = [fs.compute_global(record) for _ in range(5)]

        first_cols = list(results[0].columns)
        for result in results[1:]:
            assert list(result.columns) == first_cols

    def test_orf_global_v2_column_order_is_deterministic(self):
        """Test that column order in compute_orf_global_v2 is deterministic."""
        fs = FeatureSet(ToyResidueFeature(), name="test")
        record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="seq1")

        # Run multiple times and check order is consistent
        results = [fs.compute_orf_global_v2(record, orf=(0, 15)) for _ in range(5)]

        first_cols = list(results[0].columns)
        for result in results[1:]:
            assert list(result.columns) == first_cols

    def test_global_feature_columns_sorted(self):
        """Test that feature columns are sorted alphabetically."""

        class MultiFeature:
            """Feature that returns multiple keys in non-alphabetical order."""

            def __call__(self, record: SeqRecord) -> dict[str, float]:
                # Return in non-alphabetical order
                return {
                    "zebra": 1.0,
                    "alpha": 2.0,
                    "charlie": 3.0,
                }

        fs = FeatureSet(MultiFeature(), name="multi")
        record = SeqRecord(Seq("ATGCGCATGC"), id="seq1")

        result = fs.compute_global(record)

        # Get feature columns (skip metadata)
        feature_cols = [c for c in result.columns if c.startswith("multi.")]

        # Should be sorted alphabetically by feature key
        assert feature_cols == [
            "multi.alpha_GLOBAL",
            "multi.charlie_GLOBAL",
            "multi.zebra_GLOBAL",
        ]
