"""Integration tests for IDRPred feature.

These tests verify that IDRPredConsensusMask works correctly when the idrpred
CLI is available. Tests are skipped if idrpred is not installed.
"""

import shutil

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import PositionSpace

# Skip all tests if idrpred CLI is not available
pytestmark = pytest.mark.skipif(
    shutil.which("idrpred") is None,
    reason="idrpred CLI not available on PATH",
)


class TestIDRPredConsensusMaskBasic:
    """Basic tests for IDRPredConsensusMask feature."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the feature."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        assert feature is not None

    def test_position_space_is_residue(self):
        """Test that feature operates in RESIDUE position space."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        assert feature.position_space == PositionSpace.RESIDUE

    def test_vector_keys_contains_idrpred_idr(self):
        """Test that vector_keys contains IDRPRED_IDR with mean aggregation."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        vector_keys = feature.vector_keys
        assert "IDRPRED_IDR" in vector_keys
        assert vector_keys["IDRPRED_IDR"].aggregation_fn == np.mean

    def test_compute_vector_returns_dict_with_idrpred_idr(self):
        """Test that compute_vector returns a dictionary with IDRPRED_IDR key."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        assert isinstance(result, dict)
        assert "IDRPRED_IDR" in result
        assert isinstance(result["IDRPRED_IDR"], np.ndarray)

    def test_vector_length_matches_sequence_length(self):
        """Test that vector length equals protein sequence length."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()

        # Test with various protein sequence lengths
        test_sequences = [
            "MKALV",  # 5 residues
            "MKALVSWGRPQM",  # 12 residues
            "MKALVSWGRPQMTESTSEQ",  # 19 residues
        ]

        for seq in test_sequences:
            record = SeqRecord(Seq(seq), id="test")
            record.annotations["molecule_type"] = "protein"
            result = feature.compute_vector(record)

            assert len(result["IDRPRED_IDR"]) == len(seq), (
                f"Expected IDRPRED_IDR vector length {len(seq)}, "
                f"got {len(result['IDRPRED_IDR'])}"
            )

    def test_idrpred_idr_values_are_binary(self):
        """Test that IDRPRED_IDR values are binary (0 or 1)."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        idrpred_idr = result["IDRPRED_IDR"]
        # Values should be in {0, 1}
        assert np.all((idrpred_idr == 0) | (idrpred_idr == 1)), (
            "IDRPRED_IDR values should be binary (0 or 1)"
        )

    def test_idrpred_idr_dtype_is_float(self):
        """Test that IDRPRED_IDR array has float dtype."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        idrpred_idr = result["IDRPRED_IDR"]
        assert idrpred_idr.dtype in (np.float32, np.float64), (
            f"Expected float dtype, got {idrpred_idr.dtype}"
        )

    def test_empty_sequence_returns_empty_vector(self):
        """Test that empty sequence returns empty IDRPRED_IDR vector."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        assert len(result["IDRPRED_IDR"]) == 0


class TestIDRPredConsensusMaskDeterminism:
    """Test determinism of IDRPred prediction."""

    def test_same_sequence_produces_same_result(self):
        """Test that the same sequence produces identical results."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        seq = "MKALVSWGRPQMTEST"

        # Compute twice with the same sequence
        record1 = SeqRecord(Seq(seq), id="test1")
        record1.annotations["molecule_type"] = "protein"
        result1 = feature.compute_vector(record1)

        record2 = SeqRecord(Seq(seq), id="test2")
        record2.annotations["molecule_type"] = "protein"
        result2 = feature.compute_vector(record2)

        # Results should be identical
        np.testing.assert_array_equal(
            result1["IDRPRED_IDR"],
            result2["IDRPRED_IDR"],
            err_msg="Same sequence should produce identical IDR predictions",
        )

    def test_multiple_calls_on_same_feature_instance(self):
        """Test that multiple calls on the same feature instance are deterministic."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        seq = "MKALVSWGR"
        record = SeqRecord(Seq(seq), id="test")
        record.annotations["molecule_type"] = "protein"

        # Call compute_vector multiple times
        results = [feature.compute_vector(record) for _ in range(3)]

        # All results should be identical
        for i in range(1, len(results)):
            np.testing.assert_array_equal(
                results[0]["IDRPRED_IDR"],
                results[i]["IDRPRED_IDR"],
                err_msg=f"Call {i} produced different result than call 0",
            )


class TestIDRPredConsensusMaskTranslation:
    """Test translation behavior for DNA/RNA inputs."""

    def test_dna_input_is_translated(self):
        """Test that DNA input is automatically translated to protein."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()

        # Create DNA record
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"  # Encodes MKALVSWGR
        dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
        dna_record.annotations["molecule_type"] = "DNA"

        result = feature.compute_vector(dna_record)

        # Should have length matching translated protein (9 amino acids)
        expected_length = len(dna_seq) // 3
        assert len(result["IDRPRED_IDR"]) == expected_length

    def test_rna_input_is_translated(self):
        """Test that RNA input is automatically translated to protein."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()

        # Create RNA record (same sequence as DNA test but with U instead of T)
        rna_seq = "AUGAAAGCCCUGGUGUCCUGGGGCCGC"  # Encodes MKALVSWGR
        rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
        rna_record.annotations["molecule_type"] = "RNA"

        result = feature.compute_vector(rna_record)

        # Should have length matching translated protein (9 amino acids)
        expected_length = len(rna_seq) // 3
        assert len(result["IDRPRED_IDR"]) == expected_length

    def test_dna_translation_equals_direct_protein(self):
        """Test that DNA->protein path equals direct protein path."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()

        # Define protein sequence and its encoding DNA
        protein_seq = "MKALVSWGR"
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"  # Encodes MKALVSWGR

        # Compute on protein directly
        protein_record = SeqRecord(Seq(protein_seq), id="test_protein")
        protein_record.annotations["molecule_type"] = "protein"
        protein_result = feature.compute_vector(protein_record)

        # Compute on DNA (will be translated)
        dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
        dna_record.annotations["molecule_type"] = "DNA"
        dna_result = feature.compute_vector(dna_record)

        # Results should be identical
        np.testing.assert_array_equal(
            protein_result["IDRPRED_IDR"],
            dna_result["IDRPRED_IDR"],
            err_msg="DNA translation path should produce same result as direct protein",
        )

    def test_rna_translation_equals_direct_protein(self):
        """Test that RNA->protein path equals direct protein path."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()

        # Define protein sequence and its encoding RNA
        protein_seq = "MKALVSWGR"
        rna_seq = "AUGAAAGCCCUGGUGUCCUGGGGCCGC"  # Encodes MKALVSWGR

        # Compute on protein directly
        protein_record = SeqRecord(Seq(protein_seq), id="test_protein")
        protein_record.annotations["molecule_type"] = "protein"
        protein_result = feature.compute_vector(protein_record)

        # Compute on RNA (will be translated)
        rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
        rna_record.annotations["molecule_type"] = "RNA"
        rna_result = feature.compute_vector(rna_record)

        # Results should be identical
        np.testing.assert_array_equal(
            protein_result["IDRPRED_IDR"],
            rna_result["IDRPRED_IDR"],
            err_msg="RNA translation path should produce same result as direct protein",
        )

    def test_translation_with_stop_codon(self):
        """Test translation behavior with terminal stop codon."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        # Default behavior strips terminal stop
        feature = IDRPredConsensusMask(strip_terminal_stop=True)

        # DNA with terminal stop codon
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGCTAA"  # MKALVSWGR + stop
        dna_record = SeqRecord(Seq(dna_seq), id="test")
        dna_record.annotations["molecule_type"] = "DNA"

        result = feature.compute_vector(dna_record)

        # Should have 9 residues (stop stripped)
        assert len(result["IDRPRED_IDR"]) == 9


class TestIDRPredConsensusMaskExport:
    """Test that the feature is properly exported."""

    def test_feature_in_get_features(self):
        """Test that IDRPredConsensusMask is in get_features() list."""
        from biotooler.families.disorder import get_features
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        features = get_features()
        assert IDRPredConsensusMask in features

    def test_get_features_returns_list_of_classes(self):
        """Test that get_features() returns a list of classes, not instances."""
        from biotooler.families.disorder import get_features

        features = get_features()
        assert isinstance(features, list)
        assert len(features) > 0
        # Check that items are classes, not instances
        for feature_class in features:
            assert isinstance(feature_class, type)


class TestIDRPredConsensusMaskCoexistence:
    """Test that IDRPredConsensusMask can coexist with metapredict features."""

    def test_output_key_is_distinct_from_disorder_p(self):
        """Test that IDRPRED_IDR is distinct from DISORDER_P."""
        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask

        feature = IDRPredConsensusMask()
        vector_keys = feature.vector_keys

        # Should have IDRPRED_IDR, not DISORDER_P
        assert "IDRPRED_IDR" in vector_keys
        assert "DISORDER_P" not in vector_keys

    def test_can_compute_both_metapredict_and_idrpred(self):
        """Test that both metapredict and idrpred features can be used together."""
        # This test requires metapredict to be installed
        pytest.importorskip("metapredict", reason="metapredict not installed")

        from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        # Create both features
        metapredict_feature = DisorderProfileMetapredict()
        idrpred_feature = IDRPredConsensusMask()

        # Test on same record
        record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="test")
        record.annotations["molecule_type"] = "protein"

        # Compute both features
        metapredict_result = metapredict_feature.compute_vector(record)
        idrpred_result = idrpred_feature.compute_vector(record)

        # Both should return results with distinct keys
        assert "DISORDER_P" in metapredict_result
        assert "IDRPRED_IDR" in idrpred_result

        # Both should have same length (protein sequence length)
        assert len(metapredict_result["DISORDER_P"]) == len(record.seq)
        assert len(idrpred_result["IDRPRED_IDR"]) == len(record.seq)

        # Values should be different (one is continuous probability, one is binary)
        assert metapredict_result["DISORDER_P"].dtype == np.float64
        assert idrpred_result["IDRPRED_IDR"].dtype == np.float64
