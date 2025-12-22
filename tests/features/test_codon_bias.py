"""Tests for codon bias feature computation."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from codonbias.scores import (
    CodonAdaptationIndex,
    EffectiveNumberOfCodons,
    FrequencyOfOptimalCodons,
)

from biotooler.families.codon_bias import CodonBiasFeature
from biotooler.features import FeatureSet


class TestCodonBiasBasic:
    """Tests for basic codon bias feature computation."""

    def test_codon_bias_scalar_matches_codonbias_direct(self):
        """Test that CodonBiasFeature matches direct codonbias calls."""
        # Create reference sequence for CAI
        ref_seq = "ATGATGATGATGATGATGATGATGATGATGATG"

        # Create models
        cai = CodonAdaptationIndex(ref_seq)
        enc = EffectiveNumberOfCodons()

        # Create feature
        feature = CodonBiasFeature([cai, enc], names=["CAI", "ENC"])

        # Test sequence
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute with feature
        result = feature(record)

        # Compute directly with codonbias
        expected_cai = cai.get_score(test_seq)
        expected_enc = enc.get_score(test_seq)

        # Check results match
        assert "CAI" in result
        assert "ENC" in result
        assert abs(result["CAI"] - expected_cai) < 1e-6
        assert abs(result["ENC"] - expected_enc) < 1e-6

    def test_codon_bias_with_default_names(self):
        """Test that default names are derived from model class names."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        enc = EffectiveNumberOfCodons()

        feature = CodonBiasFeature([cai, enc])

        test_seq = "ATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")
        result = feature(record)

        # Check that names match model class names
        assert "CodonAdaptationIndex" in result
        assert "EffectiveNumberOfCodons" in result

    def test_names_length_must_match_models(self):
        """Test that names must match length of models."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        enc = EffectiveNumberOfCodons()

        with pytest.raises(ValueError, match="Length of names.*must match"):
            CodonBiasFeature([cai, enc], names=["CAI"])


class TestCodonBiasWindowing:
    """Tests for codon bias feature computation in windowing mode."""

    def test_baseline_window_matches_direct_slicing_calls(self):
        """Test that baseline windowing matches direct slice calls."""
        # Create models
        ref_seq = "ATGATGATGATGATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)

        # Create feature set
        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        # Test sequence with multiple windows
        test_seq = "ATGATGATGATGATGATGATGATGATGATG"  # 10 codons = 30 nt
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute windowed features
        result = fs.compute_orf_windows(record, orf=(0, 30), window_nt=9, step_nt=3)

        # Manually compute expected values using slice
        expected_0 = cai.get_score(test_seq, slice=slice(0, 9))
        expected_3 = cai.get_score(test_seq, slice=slice(3, 12))
        expected_6 = cai.get_score(test_seq, slice=slice(6, 15))

        # Check results
        assert "cb.CAI_0" in result.columns
        assert "cb.CAI_3" in result.columns
        assert "cb.CAI_6" in result.columns
        assert abs(result["cb.CAI_0"].iloc[0] - expected_0) < 1e-6
        assert abs(result["cb.CAI_3"].iloc[0] - expected_3) < 1e-6
        assert abs(result["cb.CAI_6"].iloc[0] - expected_6) < 1e-6

    def test_rolling_matches_baseline_for_models_with_weights(self):
        """Test that rolling computation matches baseline for models with weights."""
        # CAI has log_weights, so should use incremental computation
        ref_seq = "ATGATGATGATGATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)

        # Verify CAI has log_weights
        assert hasattr(cai, "log_weights")

        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATGATGATGATGATGATGATGATG"  # 10 codons
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute with windowing (uses incremental if available)
        result = fs.compute_orf_windows(record, orf=(0, 30), window_nt=9, step_nt=3)

        # Compute expected values directly
        expected_0 = cai.get_score(test_seq[0:9])
        expected_3 = cai.get_score(test_seq[3:12])
        expected_6 = cai.get_score(test_seq[6:15])

        # Check with tolerance (incremental may have small floating point differences)
        assert abs(result["cb.CAI_0"].iloc[0] - expected_0) < 1e-3
        assert abs(result["cb.CAI_3"].iloc[0] - expected_3) < 1e-3
        assert abs(result["cb.CAI_6"].iloc[0] - expected_6) < 1e-3

    def test_rolling_falls_back_to_baseline_if_weights_unavailable(self):
        """Test that models without weights fall back to baseline in rolling mode."""
        # ENC doesn't have weights or log_weights attributes
        enc = EffectiveNumberOfCodons()

        # Verify ENC doesn't have log_weights or weights
        assert not hasattr(enc, "log_weights")
        assert not hasattr(enc, "weights")

        feature = CodonBiasFeature([enc], names=["ENC"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATGATGATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # This should work even though ENC doesn't have weights
        # (it falls back to baseline get_score in emit)
        result = fs.compute_orf_windows(record, orf=(0, 30), window_nt=9, step_nt=3)

        # Verify results are computed
        assert "cb.ENC_0" in result.columns
        assert "cb.ENC_3" in result.columns
        assert "cb.ENC_6" in result.columns

        # Check values are reasonable (ENC typically ranges from 20 to 61)
        assert 20 <= result["cb.ENC_0"].iloc[0] <= 61
        assert 20 <= result["cb.ENC_3"].iloc[0] <= 61
        assert 20 <= result["cb.ENC_6"].iloc[0] <= 61


class TestCodonBiasValidation:
    """Tests for validation and error handling."""

    def test_protein_alphabet_raises_error(self):
        """Test that protein sequences raise clear error."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])

        # Create protein sequence
        protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein")
        protein_record.annotations["molecule_type"] = "protein"

        with pytest.raises(
            ValueError, match="Codon bias features apply only to DNA/RNA sequences"
        ):
            feature(protein_record)

    def test_rna_u_to_t_conversion(self):
        """Test that RNA sequences are converted from U to T."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])

        # Create RNA sequence with U
        rna_seq = "AUGUGAUGUGAUG"
        rna_record = SeqRecord(Seq(rna_seq), id="rna")

        # Compute features
        result = feature(rna_record)

        # Compare with DNA version
        dna_seq = "ATGATGATGATG"
        dna_record = SeqRecord(Seq(dna_seq), id="dna")
        dna_result = feature(dna_record)

        # Results should be the same after conversion
        assert abs(result["CAI"] - dna_result["CAI"]) < 1e-6

    def test_window_step_must_be_multiple_of_3(self):
        """Test that window step must be multiple of 3 for codon windows."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # step_nt=2 is not a multiple of 3
        with pytest.raises(ValueError, match="step_nt must be a multiple of 3"):
            fs.compute_orf_windows(record, orf=(0, 21), window_nt=9, step_nt=2)

        # step_nt=3 should work
        result = fs.compute_orf_windows(record, orf=(0, 21), window_nt=9, step_nt=3)
        assert result is not None

    def test_protein_in_windowing_raises_error(self):
        """Test that protein sequences raise error in windowing mode."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        # Create protein sequence (short is sufficient for validation test)
        protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein")
        protein_record.annotations["molecule_type"] = "protein"

        # Error is raised by FeatureSet.compute_orf_windows before reaching CodonBiasFeature
        with pytest.raises(
            ValueError, match="ORF window computation is only supported for DNA/RNA"
        ):
            fs.compute_orf_windows(protein_record, orf=(0, 9), window_nt=9, step_nt=3)


class TestCodonBiasMultipleModels:
    """Tests for using multiple codon bias models together."""

    def test_multiple_models_in_single_feature(self):
        """Test computing multiple codon bias metrics simultaneously."""
        ref_seq = "ATGATGATGATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        enc = EffectiveNumberOfCodons()
        fop = FrequencyOfOptimalCodons(ref_seq)

        feature = CodonBiasFeature([cai, enc, fop], names=["CAI", "ENC", "FOP"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATGATGATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        result = fs.compute_orf_windows(record, orf=(0, 30), window_nt=12, step_nt=6)

        # Check all models produced output
        assert "cb.CAI_0" in result.columns
        assert "cb.ENC_0" in result.columns
        assert "cb.FOP_0" in result.columns
        assert "cb.CAI_6" in result.columns
        assert "cb.ENC_6" in result.columns
        assert "cb.FOP_6" in result.columns

    def test_window_naming_uses_absolute_start_index(self):
        """Test that window suffixes use absolute nucleotide start position."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        result = fs.compute_orf_windows(record, orf=(0, 21), window_nt=9, step_nt=3)

        # Check column names use absolute start positions
        # Should have windows at positions 0, 3, 6, 9, 12
        assert "cb.CAI_0" in result.columns
        assert "cb.CAI_3" in result.columns
        assert "cb.CAI_6" in result.columns
        assert "cb.CAI_9" in result.columns
        assert "cb.CAI_12" in result.columns


class TestCodonBiasEdgeCases:
    """Tests for edge cases and special scenarios."""

    def test_single_window_computation(self):
        """Test computation with only one window."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])
        fs = FeatureSet(feature, name="cb")

        test_seq = "ATGATGATG"  # Exactly one window
        record = SeqRecord(Seq(test_seq), id="test")

        result = fs.compute_orf_windows(record, orf=(0, 9), window_nt=9, step_nt=3)

        # Should have only one window at position 0
        feature_cols = [c for c in result.columns if c.startswith("cb.CAI_")]
        assert len(feature_cols) == 1
        assert "cb.CAI_0" in result.columns

    def test_lowercase_sequence(self):
        """Test that lowercase sequences are handled correctly."""
        ref_seq = "atgatgatgatg"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])

        test_seq = "atgatgatg"
        record = SeqRecord(Seq(test_seq), id="test")

        # Should not raise an error
        result = feature(record)
        assert "CAI" in result

    def test_mixed_case_rna(self):
        """Test mixed case RNA with U/u conversion."""
        ref_seq = "ATGATGATGATG"
        cai = CodonAdaptationIndex(ref_seq)
        feature = CodonBiasFeature([cai], names=["CAI"])

        # Mixed case RNA with both uppercase and lowercase U
        rna_seq = "AUGugAUGuGAUG"
        rna_record = SeqRecord(Seq(rna_seq), id="rna")

        result = feature(rna_record)

        # Should work and produce result
        assert "CAI" in result
        assert isinstance(result["CAI"], float)


class TestCodonBiasFromReference:
    """Tests for CodonBiasFeature.from_reference factory method."""

    def test_from_reference_builds_and_computes(self):
        """Test that from_reference builds models and computes correctly."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        # Create reference set
        ref_set = ReferenceSequenceSet(
            cds={
                "gene1": "ATGATGATGATGATGATGATG",
                "gene2": "ATGATGATGATGATGATGATG",
            }
        )

        # Build feature from reference
        feature = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"])

        # Test sequence
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute features
        result = feature(record)

        # Check results
        assert "CAI" in result
        assert "ENC" in result
        assert isinstance(result["CAI"], float)
        assert isinstance(result["ENC"], float)
        assert 0 <= result["CAI"] <= 1  # CAI is typically 0-1
        assert 20 <= result["ENC"] <= 61  # ENC is typically 20-61

    def test_from_reference_matches_direct_codonbias_usage(self):
        """Test that from_reference produces same results as direct codonbias usage."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        # Create reference set
        ref_cds = {
            "gene1": "ATGATGATGATGATGATGATG",
            "gene2": "ATGATGATGATGATGATGATG",
        }
        ref_set = ReferenceSequenceSet(cds=ref_cds)

        # Build feature from reference
        feature = CodonBiasFeature.from_reference(ref_set, ["CAI"])

        # Build feature directly with codonbias
        concatenated_ref = "".join(ref_set.cds_strings())
        cai_direct = CodonAdaptationIndex(concatenated_ref)
        feature_direct = CodonBiasFeature([cai_direct], names=["CAI"])

        # Test sequence
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute with both
        result_from_ref = feature(record)
        result_direct = feature_direct(record)

        # Results should match
        assert abs(result_from_ref["CAI"] - result_direct["CAI"]) < 1e-6

    def test_from_reference_with_abbreviations(self):
        """Test resolving score abbreviations."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Use various abbreviations
        feature = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC", "FOP"])

        # Check that names are set correctly (should use abbreviations)
        assert feature.names == ["CAI", "ENC", "FOP"]

        # Check that models are instantiated
        assert len(feature.models) == 3

    def test_from_reference_with_class_names(self):
        """Test resolving full class names."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Use full class names
        feature = CodonBiasFeature.from_reference(
            ref_set, ["CodonAdaptationIndex", "EffectiveNumberOfCodons"]
        )

        # Check that names are set correctly (should use abbreviations)
        assert feature.names == ["CAI", "ENC"]

        # Check that models are instantiated
        assert len(feature.models) == 2

    def test_from_reference_with_class_objects(self):
        """Test resolving class objects."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Use class objects
        feature = CodonBiasFeature.from_reference(
            ref_set, [CodonAdaptationIndex, EffectiveNumberOfCodons]
        )

        # Check that names are set correctly
        assert feature.names == ["CAI", "ENC"]

        # Check that models are instantiated
        assert len(feature.models) == 2
        assert isinstance(feature.models[0], CodonAdaptationIndex)
        assert isinstance(feature.models[1], EffectiveNumberOfCodons)

    def test_from_reference_with_custom_names(self):
        """Test providing custom names."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Provide custom names
        feature = CodonBiasFeature.from_reference(
            ref_set, ["CAI", "ENC"], names=["CodonAdaptIndex", "EffectiveNumCodons"]
        )

        # Check custom names are used
        assert feature.names == ["CodonAdaptIndex", "EffectiveNumCodons"]

    def test_from_reference_with_score_kwargs(self):
        """Test passing kwargs to score constructors."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Pass kwargs for CAI
        feature = CodonBiasFeature.from_reference(
            ref_set,
            ["CAI"],
            score_kwargs={"CAI": {"genetic_code": 1, "k_mer": 1}},
        )

        # Check model was instantiated (no error means kwargs were accepted)
        assert len(feature.models) == 1

    def test_from_reference_empty_cds_raises_error(self):
        """Test that empty CDS raises clear error."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        # Create reference set with empty CDS
        ref_set = ReferenceSequenceSet(cds={})

        # Should raise ValueError
        with pytest.raises(
            ValueError,
            match="ReferenceSequenceSet must contain CDS sequences",
        ):
            CodonBiasFeature.from_reference(ref_set, ["CAI"])

    def test_from_reference_invalid_score_identifier(self):
        """Test that invalid score identifier raises clear error."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Use invalid identifier
        with pytest.raises(ValueError, match="Cannot resolve score identifier"):
            CodonBiasFeature.from_reference(ref_set, ["INVALID_SCORE"])

    def test_from_reference_enc_without_ref_seq(self):
        """Test that ENC works without ref_seq (doesn't require reference)."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # ENC doesn't require ref_seq
        feature = CodonBiasFeature.from_reference(ref_set, ["ENC"])

        # Should work
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")
        result = feature(record)

        assert "ENC" in result
        assert isinstance(result["ENC"], float)

    def test_from_reference_mixed_score_types(self):
        """Test mixing abbreviations, class names, and class objects."""
        from biotooler.core.reference_sequences import ReferenceSequenceSet

        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG"}
        )

        # Mix different score identifier types
        feature = CodonBiasFeature.from_reference(
            ref_set,
            ["CAI", "EffectiveNumberOfCodons", FrequencyOfOptimalCodons],
            names=["CAI", "ENC", "FOP"],
        )

        # Check all models instantiated correctly
        assert len(feature.models) == 3
        assert feature.names == ["CAI", "ENC", "FOP"]

        # Test computation works
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")
        result = feature(record)

        assert "CAI" in result
        assert "ENC" in result
        assert "FOP" in result
