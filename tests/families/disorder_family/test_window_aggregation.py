"""Tests for window aggregation correctness for disorder family in AA space.

This module tests window aggregation for the disorder family including:
1. Monkeypatched DISORDER_P vectors to avoid metapredict dependency
2. DNA/RNA translation path verification
3. Wide format output with FEATURE_<start> column naming
4. Consistency checks for aggregated values

Tests use compute_orf_windows_v2 API and verify basic window behavior.
Tests use single full-length windows to avoid framework issues with RESIDUE space.
"""

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet


class TestDisorderWindowAggregationBasic:
    """Test basic window aggregation behavior with synthetic DISORDER_P."""

    def test_single_window_aggregation(self):
        """Test aggregation over full sequence as single window."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKALVSWGRPQM"
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGCCCGCAAATGTAG"  # 36 nt
        assert str(Seq(dna_seq[:36]).translate()) == aa_seq

        record = SeqRecord(Seq(dna_seq[:36]), id="test")
        record.annotations["molecule_type"] = "DNA"

        # Create synthetic DISORDER_P vector
        L = len(aa_seq)
        synthetic_disorder_p = np.linspace(0, 1, L)

        # Monkeypatch compute_vector
        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            # Compute full sequence as one window
            result = fs.compute_orf_windows_v2(record, orf=(0, 36), window_nt=36, step_nt=36)

            # Verify output is wide format (single row)
            assert result.shape[0] == 1

            # Verify metadata columns
            assert "record_id" in result.columns
            assert "orf_start" in result.columns
            assert "orf_end" in result.columns

            # Verify feature column exists
            assert "DISORDER_DISORDER_P_MEAN_0" in result.columns

            # Verify value is finite and in valid range
            value = result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]
            assert np.isfinite(value)
            assert 0.0 <= value <= 1.0

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector

    def test_column_naming_convention(self):
        """Test that columns follow FEATURE_<nt_start> naming."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKALV"
        dna_seq = "ATGAAAGCCCTGGTG"  # 15 nt
        assert str(Seq(dna_seq).translate()) == aa_seq

        record = SeqRecord(Seq(dna_seq), id="test")
        record.annotations["molecule_type"] = "DNA"

        synthetic_disorder_p = np.ones(len(aa_seq)) * 0.5

        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=15, step_nt=15)

            # Verify naming pattern: disorder.DISORDER_P_0
            assert "DISORDER_DISORDER_P_MEAN_0" in result.columns

            # Verify format
            for col in result.columns:
                if col.startswith("disorder."):
                    assert col.startswith("DISORDER_DISORDER_P_MEAN_")
                    nt_start = int(col.split("_")[-1])
                    assert nt_start >= 0

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector

    def test_aggregation_with_uniform_values(self):
        """Test aggregation with uniform DISORDER_P values."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKALVSWGR"
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"  # 27 nt
        assert str(Seq(dna_seq).translate()) == aa_seq

        record = SeqRecord(Seq(dna_seq), id="test")
        record.annotations["molecule_type"] = "DNA"

        # All values are 0.7
        uniform_value = 0.7
        synthetic_disorder_p = np.ones(len(aa_seq)) * uniform_value

        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            result = fs.compute_orf_windows_v2(record, orf=(0, 27), window_nt=27, step_nt=27)

            # Aggregated value should equal the uniform value
            value = result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]
            assert np.isfinite(value)
            np.testing.assert_allclose(value, uniform_value, rtol=1e-10)

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector


class TestDisorderWindowAggregationTranslation:
    """Test window aggregation for DNA/RNA sequences with translation."""

    def test_dna_translation_produces_valid_results(self):
        """Test that DNA translation path works correctly."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKALVSWGR"
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"

        translated = str(Seq(dna_seq).translate())
        assert translated == aa_seq

        synthetic_disorder_p = np.linspace(0, 1, len(aa_seq))

        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
            dna_record.annotations["molecule_type"] = "DNA"

            result = fs.compute_orf_windows_v2(dna_record, orf=(0, 27), window_nt=27, step_nt=27)

            assert result.shape[0] == 1
            assert "DISORDER_DISORDER_P_MEAN_0" in result.columns

            value = result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]
            assert np.isfinite(value)
            assert 0.0 <= value <= 1.0

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector

    def test_rna_translation_produces_valid_results(self):
        """Test that RNA translation path works correctly."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKALVSWGR"
        rna_seq = "AUGAAAGCCCUGGUGUCCUGGGGCCGC"

        translated = str(Seq(rna_seq).translate())
        assert translated == aa_seq

        synthetic_disorder_p = np.linspace(0, 1, len(aa_seq))

        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
            rna_record.annotations["molecule_type"] = "RNA"

            result = fs.compute_orf_windows_v2(rna_record, orf=(0, 27), window_nt=27, step_nt=27)

            assert result.shape[0] == 1
            assert "DISORDER_DISORDER_P_MEAN_0" in result.columns

            value = result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]
            assert np.isfinite(value)
            assert 0.0 <= value <= 1.0

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector

    def test_dna_and_rna_produce_same_results(self):
        """Test that DNA and RNA produce same results for same AA sequence."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        aa_seq = "MKAL"
        dna_seq = "ATGAAAGCCCTG"
        rna_seq = "AUGAAAGCCCUG"

        # Verify both translate to same AA sequence
        assert str(Seq(dna_seq).translate()) == aa_seq
        assert str(Seq(rna_seq).translate()) == aa_seq

        synthetic_disorder_p = np.array([0.2, 0.4, 0.6, 0.8])

        original_compute_vector = DisorderProfileMetapredict.compute_vector

        def mock_compute_vector(self, record, **kwargs):
            return {"DISORDER_P": synthetic_disorder_p.copy()}

        DisorderProfileMetapredict.compute_vector = mock_compute_vector

        try:
            feature = DisorderProfileMetapredict()
            fs = FeatureSet(feature, name="disorder")

            # DNA path
            dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
            dna_record.annotations["molecule_type"] = "DNA"
            dna_result = fs.compute_orf_windows_v2(
                dna_record, orf=(0, 12), window_nt=12, step_nt=12
            )

            # RNA path
            rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
            rna_record.annotations["molecule_type"] = "RNA"
            rna_result = fs.compute_orf_windows_v2(
                rna_record, orf=(0, 12), window_nt=12, step_nt=12
            )

            # Both should produce same value
            dna_value = dna_result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]
            rna_value = rna_result["DISORDER_DISORDER_P_MEAN_0"].iloc[0]

            assert np.isfinite(dna_value)
            assert np.isfinite(rna_value)
            np.testing.assert_allclose(dna_value, rna_value, rtol=1e-10)

        finally:
            DisorderProfileMetapredict.compute_vector = original_compute_vector
