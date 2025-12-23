"""Tests for metapredict backend disorder prediction feature.

Tests verify:
1. Skip cleanly if metapredict not installed
2. AA input: length matches, dtype float, min>=0, max<=1
3. DNA/RNA path equals AA direct path on translated sequence
4. Determinism for same input
"""

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features.aggregation import PositionSpace

# Try to import metapredict - skip tests if not available
pytest.importorskip("metapredict", reason="metapredict not installed")


class TestDisorderProfileMetapredictBasic:
    """Basic tests for DisorderProfileMetapredict feature."""

    def test_import_and_instantiate(self):
        """Test that we can import and instantiate the feature."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        assert feature is not None

    def test_position_space_is_residue(self):
        """Test that feature operates in RESIDUE position space."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        assert feature.position_space == PositionSpace.RESIDUE

    def test_vector_keys_contains_disorder_p(self):
        """Test that vector_keys contains DISORDER_P with mean aggregation."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        vector_keys = feature.vector_keys
        assert "DISORDER_P" in vector_keys
        assert vector_keys["DISORDER_P"].aggregation_fn == np.mean

    def test_compute_vector_returns_dict_with_disorder_p(self):
        """Test that compute_vector returns a dictionary with DISORDER_P key."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        assert isinstance(result, dict)
        assert "DISORDER_P" in result
        assert isinstance(result["DISORDER_P"], np.ndarray)

    def test_vector_length_matches_sequence_length(self):
        """Test that vector length equals protein sequence length."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()

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

            assert len(result["DISORDER_P"]) == len(seq), (
                f"Expected DISORDER_P vector length {len(seq)}, "
                f"got {len(result['DISORDER_P'])}"
            )

    def test_disorder_p_values_in_valid_range(self):
        """Test that DISORDER_P values are probabilities in [0, 1]."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        disorder_p = result["DISORDER_P"]
        assert np.all(disorder_p >= 0.0), "DISORDER_P values should be >= 0"
        assert np.all(disorder_p <= 1.0), "DISORDER_P values should be <= 1"

    def test_disorder_p_dtype_is_float(self):
        """Test that DISORDER_P array has float dtype."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        record = SeqRecord(Seq("MKALVSWGR"), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        disorder_p = result["DISORDER_P"]
        assert disorder_p.dtype in (np.float32, np.float64), (
            f"Expected float dtype, got {disorder_p.dtype}"
        )

    def test_empty_sequence_returns_empty_vector(self):
        """Test that empty sequence returns empty DISORDER_P vector."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        record = SeqRecord(Seq(""), id="test")
        record.annotations["molecule_type"] = "protein"
        result = feature.compute_vector(record)

        assert len(result["DISORDER_P"]) == 0


class TestDisorderProfileMetapredictDeterminism:
    """Test determinism of metapredict disorder prediction."""

    def test_same_sequence_produces_same_result(self):
        """Test that the same sequence produces identical results."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
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
            result1["DISORDER_P"],
            result2["DISORDER_P"],
            err_msg="Same sequence should produce identical disorder predictions",
        )

    def test_multiple_calls_on_same_feature_instance(self):
        """Test that multiple calls on the same feature instance are deterministic."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()
        seq = "MKALVSWGR"
        record = SeqRecord(Seq(seq), id="test")
        record.annotations["molecule_type"] = "protein"

        # Call compute_vector multiple times
        results = [feature.compute_vector(record) for _ in range(3)]

        # All results should be identical
        for i in range(1, len(results)):
            np.testing.assert_array_equal(
                results[0]["DISORDER_P"],
                results[i]["DISORDER_P"],
                err_msg=f"Call {i} produced different result than call 0",
            )


class TestDisorderProfileMetapredictTranslation:
    """Test translation behavior for DNA/RNA inputs."""

    def test_dna_input_is_translated(self):
        """Test that DNA input is automatically translated to protein."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()

        # Create DNA record
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGC"  # Encodes MKALVSWGR
        dna_record = SeqRecord(Seq(dna_seq), id="test_dna")
        dna_record.annotations["molecule_type"] = "DNA"

        result = feature.compute_vector(dna_record)

        # Should have length matching translated protein (9 amino acids)
        expected_length = len(dna_seq) // 3
        assert len(result["DISORDER_P"]) == expected_length

    def test_rna_input_is_translated(self):
        """Test that RNA input is automatically translated to protein."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()

        # Create RNA record (same sequence as DNA test but with U instead of T)
        rna_seq = "AUGAAAGCCCUGGUGSCCUGGGGCCGC"  # Encodes MKALVSWGR
        rna_record = SeqRecord(Seq(rna_seq), id="test_rna")
        rna_record.annotations["molecule_type"] = "RNA"

        result = feature.compute_vector(rna_record)

        # Should have length matching translated protein (9 amino acids)
        expected_length = len(rna_seq) // 3
        assert len(result["DISORDER_P"]) == expected_length

    def test_dna_translation_equals_direct_protein(self):
        """Test that DNA->protein path equals direct protein path."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()

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
            protein_result["DISORDER_P"],
            dna_result["DISORDER_P"],
            err_msg="DNA translation path should produce same result as direct protein",
        )

    def test_rna_translation_equals_direct_protein(self):
        """Test that RNA->protein path equals direct protein path."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        feature = DisorderProfileMetapredict()

        # Define protein sequence and its encoding RNA
        protein_seq = "MKALVSWGR"
        rna_seq = "AUGAAAGCCCUGGUGSCCUGGGGCCGC"  # Encodes MKALVSWGR

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
            protein_result["DISORDER_P"],
            rna_result["DISORDER_P"],
            err_msg="RNA translation path should produce same result as direct protein",
        )

    def test_translation_with_stop_codon(self):
        """Test translation behavior with terminal stop codon."""
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        # Default behavior strips terminal stop
        feature = DisorderProfileMetapredict(strip_terminal_stop=True)

        # DNA with terminal stop codon
        dna_seq = "ATGAAAGCCCTGGTGTCCTGGGGCCGCTAA"  # MKALVSWGR + stop
        dna_record = SeqRecord(Seq(dna_seq), id="test")
        dna_record.annotations["molecule_type"] = "DNA"

        result = feature.compute_vector(dna_record)

        # Should have 9 residues (stop stripped)
        assert len(result["DISORDER_P"]) == 9


class TestDisorderProfileMetapredictExport:
    """Test that the feature is properly exported."""

    def test_feature_in_get_features(self):
        """Test that DisorderProfileMetapredict is in get_features() list."""
        from biotooler.families.disorder import get_features
        from biotooler.families.disorder.metapredict_backend import (
            DisorderProfileMetapredict,
        )

        features = get_features()
        assert DisorderProfileMetapredict in features

    def test_get_features_returns_list_of_classes(self):
        """Test that get_features() returns a list of classes, not instances."""
        from biotooler.families.disorder import get_features

        features = get_features()
        assert isinstance(features, list)
        assert len(features) > 0
        # Check that items are classes, not instances
        for feature_class in features:
            assert isinstance(feature_class, type)
