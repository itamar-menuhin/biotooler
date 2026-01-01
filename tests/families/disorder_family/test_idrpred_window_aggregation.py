"""Window aggregation correctness tests for IDRPred (AA space + translation path).

This module tests window aggregation for the IDRPred feature including:
1. Monkeypatched IDRPRED_IDR mask vectors to avoid idrpred CLI dependency
2. DNA/RNA translation path verification
3. Wide format output with correct column naming
4. Manual aggregation verification for window correctness

Tests use compute_orf_windows_v2 API and verify window aggregation correctness.
Tests do NOT require idrpred CLI to be installed (using monkeypatching).

Note: IDRPred operates in RESIDUE (per-AA) space after translation. When windowing
DNA/RNA sequences, the windowing engine works in NT space but the feature returns
values in AA space (one value per 3 nt). Tests use full-length windows to avoid
complexity of partial window mapping.
"""

import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask
from biotooler.features import FeatureSet


class TestIDRPredWindowMeanMatchesManual:
    """Test that window-aggregated values equal manual aggregation for each window."""

    def test_idrpred_window_mean_matches_manual(self, monkeypatch):
        """Test window aggregation correctness with deterministic mask vector.

        This test monkeypatches IDRPredConsensusMask.compute_vector() to return
        a known mask vector for a specific protein length, then verifies that
        windowed aggregation matches manual calculation of mean(mask).

        Uses a single full-length window to test aggregation behavior.
        """
        # Define test sequence and expected mask
        # DNA: 36 nt -> 12 AA protein
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGCCCGCAAATG"  # 36 nt
        protein_seq = str(Seq(dna_seq).translate())
        protein_length = len(protein_seq)
        assert protein_length == 12, f"Expected 12 AA, got {protein_length}"

        # Create deterministic mask: alternating pattern for easy verification
        # [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
        known_mask = np.array([1.0 if i % 2 == 0 else 0.0 for i in range(protein_length)])
        assert np.sum(known_mask) == 6.0  # 6 ones
        expected_mean = 6.0 / 12.0  # 0.5

        # Monkeypatch compute_vector to return our known mask
        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        # Create feature and FeatureSet
        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        # Create DNA record
        record = SeqRecord(Seq(dna_seq), id="test_seq")
        record.annotations["molecule_type"] = "DNA"

        # Compute single window covering entire sequence
        result = fs.compute_orf_windows_v2(record, orf=(0, 36), window_nt=36, step_nt=36)

        # Verify output shape
        assert result.shape[0] == 1, "Expected single-row (wide format)"

        # Verify metadata columns
        assert result["record_id"].iloc[0] == "test_seq"
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 36

        # Verify computed value matches expected value
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in result.columns
        np.testing.assert_allclose(
            result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0],
            expected_mean,
            rtol=1e-10,
            err_msg="Window aggregation doesn't match manual calculation",
        )

    def test_idrpred_window_mean_with_varied_mask(self, monkeypatch):
        """Test window aggregation with a more complex mask pattern.

        This test uses a non-uniform mask to verify aggregation correctness
        with different IDR densities.
        """
        # DNA: 30 nt -> 10 AA protein
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGCCCG"  # 30 nt
        protein_seq = str(Seq(dna_seq).translate())
        protein_length = len(protein_seq)
        assert protein_length == 10

        # Create mask with specific pattern:
        # [1, 1, 1, 0, 0, 0, 1, 1, 0, 0]
        # 5 ones out of 10 -> mean = 0.5
        known_mask = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0])
        expected_mean = np.mean(known_mask)
        assert abs(expected_mean - 0.5) < 1e-10

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test_seq")
        record.annotations["molecule_type"] = "DNA"

        # Single window covering entire sequence
        result = fs.compute_orf_windows_v2(record, orf=(0, 30), window_nt=30, step_nt=30)

        # Verify computed value
        np.testing.assert_allclose(
            result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0],
            expected_mean,
            rtol=1e-10,
        )

    def test_idrpred_window_mean_all_idr(self, monkeypatch):
        """Test window aggregation with all positions in IDR (all 1s)."""
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGC"  # 24 nt -> 8 AA
        protein_length = 8

        # All positions are IDR
        known_mask = np.ones(protein_length)

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test_seq")
        record.annotations["molecule_type"] = "DNA"

        # Single window covering entire sequence
        result = fs.compute_orf_windows_v2(record, orf=(0, 24), window_nt=24, step_nt=24)

        # All IDR means should be 1.0
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in result.columns
        np.testing.assert_allclose(result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0], 1.0, rtol=1e-10)

    def test_idrpred_window_mean_no_idr(self, monkeypatch):
        """Test window aggregation with no IDR positions (all 0s)."""
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGC"  # 24 nt -> 8 AA
        protein_length = 8

        # No positions are IDR
        known_mask = np.zeros(protein_length)

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test_seq")
        record.annotations["molecule_type"] = "DNA"

        # Single window covering entire sequence
        result = fs.compute_orf_windows_v2(record, orf=(0, 24), window_nt=24, step_nt=24)

        # All IDR means should be 0.0
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in result.columns
        np.testing.assert_allclose(result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0], 0.0, rtol=1e-10)


class TestIDRPredTranslationPathWindowing:
    """Test DNA/RNA translation path produces identical results as direct computation."""

    def test_idrpred_translation_path_windowing(self, monkeypatch):
        """Test that DNA and RNA translation paths produce identical results.

        This test verifies that the translation path (DNA/RNA -> protein -> windowing)
        produces identical results for both DNA and RNA inputs encoding the same protein.
        """
        # Define protein sequence and its encoding DNA/RNA
        protein_seq = "MKALVSWGR"  # 9 AA
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"  # 27 nt encodes MKALVSWGR
        rna_seq = dna_seq.replace("T", "U")

        # Verify translations
        assert str(Seq(dna_seq).translate()) == protein_seq
        assert str(Seq(rna_seq).translate()) == protein_seq

        # Create deterministic mask for the protein
        known_mask = np.array([1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0])
        assert len(known_mask) == len(protein_seq)
        expected_mean = np.mean(known_mask)

        # Monkeypatch compute_vector
        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        # DNA path
        dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
        dna_record.annotations["molecule_type"] = "DNA"
        dna_result = fs.compute_orf_windows_v2(dna_record, orf=(0, 27), window_nt=27, step_nt=27)

        # RNA path
        rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
        rna_record.annotations["molecule_type"] = "RNA"
        rna_result = fs.compute_orf_windows_v2(rna_record, orf=(0, 27), window_nt=27, step_nt=27)

        # Both should have same columns
        assert set(dna_result.columns) == set(rna_result.columns)

        # Both should produce identical values for all feature columns
        for col in dna_result.columns:
            if col.startswith("idrpred."):
                dna_val = dna_result[col].iloc[0]
                rna_val = rna_result[col].iloc[0]
                np.testing.assert_allclose(
                    dna_val,
                    rna_val,
                    rtol=1e-10,
                    err_msg=f"DNA and RNA paths differ for column {col}",
                )
                # Also verify they match expected mean
                np.testing.assert_allclose(
                    dna_val,
                    expected_mean,
                    rtol=1e-10,
                )

    def test_rna_translation_equals_dna_translation(self, monkeypatch):
        """Test that RNA and DNA translation paths produce identical window results."""
        # Define sequences
        protein_seq = "MKALV"  # 5 AA
        dna_seq = "ATGAAAGCCCTGGTG"  # 15 nt
        rna_seq = "AUGAAAGCCCUGGUG"  # 15 nt

        # Verify translations
        assert str(Seq(dna_seq).translate()) == protein_seq
        assert str(Seq(rna_seq).translate()) == protein_seq

        # Create mask
        known_mask = np.array([1.0, 0.0, 1.0, 1.0, 0.0])
        expected_mean = np.mean(known_mask)

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        # DNA path
        dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
        dna_record.annotations["molecule_type"] = "DNA"
        dna_result = fs.compute_orf_windows_v2(dna_record, orf=(0, 15), window_nt=15, step_nt=15)

        # RNA path
        rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
        rna_record.annotations["molecule_type"] = "RNA"
        rna_result = fs.compute_orf_windows_v2(rna_record, orf=(0, 15), window_nt=15, step_nt=15)

        # Verify both produce same result
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in dna_result.columns
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in rna_result.columns

        dna_val = dna_result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0]
        rna_val = rna_result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0]

        np.testing.assert_allclose(dna_val, rna_val, rtol=1e-10)

        # Verify values match manual calculation
        np.testing.assert_allclose(dna_val, expected_mean, rtol=1e-10)
        np.testing.assert_allclose(rna_val, expected_mean, rtol=1e-10)

    def test_translation_preserves_mask_length(self, monkeypatch):
        """Test that translation path correctly preserves protein length in windowing."""
        # 24 nt DNA -> 8 AA protein
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGC"
        protein_seq = str(Seq(dna_seq).translate())
        assert len(protein_seq) == 8

        # Mask for 8 AA
        known_mask = np.array([1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0])
        expected_mean = np.mean(known_mask)  # 5/8 = 0.625

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        # DNA path with full-length window
        dna_record = SeqRecord(Seq(dna_seq), id="test")
        dna_record.annotations["molecule_type"] = "DNA"

        result = fs.compute_orf_windows_v2(dna_record, orf=(0, 24), window_nt=24, step_nt=24)

        # Should have single window at 0
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in result.columns

        # Verify calculation
        np.testing.assert_allclose(
            result["IDRPRED_IDRPRED_IDR_MEAN_0"].iloc[0], expected_mean, rtol=1e-10
        )


class TestIDRPredWindowOutputFormat:
    """Test that output follows correct wide format with proper column naming."""

    def test_output_is_wide_format(self, monkeypatch):
        """Test that output is in wide format (single row)."""
        dna_seq = "ATGAAAGCCCTGGTG"  # 15 nt -> 5 AA
        known_mask = np.array([1.0, 0.0, 1.0, 0.0, 1.0])

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test")
        record.annotations["molecule_type"] = "DNA"

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=15, step_nt=15)

        # Should be single row (wide format)
        assert result.shape[0] == 1

    def test_column_naming_convention(self, monkeypatch):
        """Test that columns follow FeatureSet naming convention: name.key_nt_start."""
        dna_seq = "ATGAAAGCCCTGGTG"  # 15 nt -> 5 AA
        known_mask = np.array([1.0, 0.0, 1.0, 0.0, 1.0])

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test")
        record.annotations["molecule_type"] = "DNA"

        # Single full-length window at nt position 0
        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=15, step_nt=15)

        # Check column naming pattern
        assert "IDRPRED_IDRPRED_IDR_MEAN_0" in result.columns

        # Verify no unexpected columns
        feature_cols = [c for c in result.columns if c.startswith("IDRPRED_")]
        assert len(feature_cols) == 1

    def test_metadata_columns_present(self, monkeypatch):
        """Test that metadata columns (record_id, orf_start, orf_end) are present."""
        dna_seq = "ATGAAAGCCCTGGTG"  # 15 nt -> 5 AA
        known_mask = np.array([1.0, 0.0, 1.0, 0.0, 1.0])

        def mock_compute_vector(self, record, **kwargs):
            return {"IDRPRED_IDR": known_mask.copy()}

        monkeypatch.setattr(IDRPredConsensusMask, "compute_vector", mock_compute_vector)

        feature = IDRPredConsensusMask()
        fs = FeatureSet(feature, name="idrpred")

        record = SeqRecord(Seq(dna_seq), id="test_record")
        record.annotations["molecule_type"] = "DNA"

        result = fs.compute_orf_windows_v2(record, orf=(0, 15), window_nt=15, step_nt=15)

        # Check metadata columns
        assert "record_id" in result.columns
        assert "orf_start" in result.columns
        assert "orf_end" in result.columns

        # Verify values
        assert result["record_id"].iloc[0] == "test_record"
        assert result["orf_start"].iloc[0] == 0
        assert result["orf_end"].iloc[0] == 15
