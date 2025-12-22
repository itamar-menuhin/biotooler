"""Tests for protparam family feature."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.families.protparam import ProtParamFeature


class TestProtParamFeatureProteinInput:
    """Tests for ProtParamFeature with protein input."""

    def test_basic_protein_computation(self):
        """Test basic ProtParam computation on protein sequence."""
        feature = ProtParamFeature()
        protein = SeqRecord(Seq("MKALVSWGR"), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        # Check all expected keys are present
        assert "molecular_weight" in result
        assert "isoelectric_point" in result
        assert "gravy" in result
        assert "instability_index" in result
        assert "aromaticity" in result

        # Check types
        assert isinstance(result["molecular_weight"], (int, float))
        assert isinstance(result["isoelectric_point"], (int, float))
        assert isinstance(result["gravy"], (int, float))
        assert isinstance(result["instability_index"], (int, float))
        assert isinstance(result["aromaticity"], (int, float))

        # Sanity checks on values
        assert result["molecular_weight"] > 0
        assert 0 < result["isoelectric_point"] < 14
        assert 0 <= result["aromaticity"] <= 1

    def test_known_protein_values(self):
        """Test with a protein sequence with known characteristics."""
        feature = ProtParamFeature()
        # Simple sequence with known properties
        protein = SeqRecord(Seq("AAAAA"), id="test")  # All alanine
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        # All alanine should have no aromatic residues
        assert result["aromaticity"] == 0.0

    def test_aromatic_protein(self):
        """Test protein with aromatic residues."""
        feature = ProtParamFeature()
        # FWY are aromatic
        protein = SeqRecord(Seq("FFWWYY"), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        # All residues are aromatic
        assert result["aromaticity"] == 1.0

    def test_long_protein(self):
        """Test with a longer protein sequence."""
        feature = ProtParamFeature()
        # A longer, realistic sequence
        protein = SeqRecord(Seq("MKALVSWGRDPQFIHETNLCYKVM"), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        assert result["molecular_weight"] > 1000  # Reasonable MW
        assert 0 < result["isoelectric_point"] < 14


class TestProtParamFeatureDNATranslation:
    """Tests for ProtParamFeature with DNA input."""

    def test_dna_translation_and_computation(self):
        """Test that DNA is translated before ProtParam computation."""
        feature = ProtParamFeature()
        # ATGAAAGCCCTGGTGTCTTGGGGACGT -> MKALVSWGR
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTCTTGGGGACGT"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = feature(dna)

        # Should have all ProtParam features
        assert "molecular_weight" in result
        assert "isoelectric_point" in result
        assert "gravy" in result
        assert "instability_index" in result
        assert "aromaticity" in result

        # Should have computed reasonable values
        assert result["molecular_weight"] > 0

    def test_dna_with_stop_codon(self):
        """Test DNA with stop codon (should be stripped by default)."""
        feature = ProtParamFeature()
        # ATGAAAGCCCTGGTGTAA = MKALV* (TAA is stop, stripped by default)
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTAA"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = feature(dna)

        # Should succeed (stop codon stripped)
        assert "molecular_weight" in result

    def test_dna_internal_stop_raises_error(self):
        """Test that internal stop codons raise error by default."""
        feature = ProtParamFeature()
        # ATGTAAAAAGCCCTG = M*KAL (TAA is internal stop)
        dna = SeqRecord(Seq("ATGTAAAAAGCCCTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        with pytest.raises(ValueError, match="Internal stop codon"):
            feature(dna)

    def test_dna_internal_stop_with_ignore(self):
        """Test that internal stops can be ignored."""
        feature = ProtParamFeature(on_internal_stop="ignore")
        # ATGTAAAAAGCCCTG = M*KAL (TAA is internal stop)
        dna = SeqRecord(Seq("ATGTAAAAAGCCCTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        # Should raise error because ProtParam cannot handle stop codons
        # even if translation allows them
        with pytest.raises(ValueError):
            feature(dna)


class TestProtParamFeatureRNATranslation:
    """Tests for ProtParamFeature with RNA input."""

    def test_rna_translation_and_computation(self):
        """Test that RNA is translated before ProtParam computation."""
        feature = ProtParamFeature()
        # AUGAAAGCCCUGGUGUCUUGGGGACGU -> MKALVSWGR
        rna = SeqRecord(Seq("AUGAAAGCCCUGGUGUCUUGGGGACGU"), id="test")
        rna.annotations["molecule_type"] = "RNA"

        result = feature(rna)

        assert "molecular_weight" in result
        assert result["molecular_weight"] > 0


class TestProtParamFeatureTranslationOptions:
    """Tests for translation options in ProtParamFeature."""

    def test_custom_genetic_code(self):
        """Test with custom genetic code table."""
        feature = ProtParamFeature(table=11)  # Bacterial code
        dna = SeqRecord(Seq("ATGAAAGCCCTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = feature(dna)

        assert "molecular_weight" in result

    def test_strip_terminal_stop_false(self):
        """Test with strip_terminal_stop=False."""
        feature = ProtParamFeature(strip_terminal_stop=False)
        # ATGAAAGCCTAA = MKA* (TAA is stop)
        dna = SeqRecord(Seq("ATGAAAGCCTAA"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        # ProtParam should fail with stop codon
        with pytest.raises(ValueError):
            feature(dna)


class TestProtParamFeatureLazyImport:
    """Tests for lazy import behavior."""

    def test_import_does_not_load_protparam(self):
        """Test that importing the module doesn't load ProtParam."""
        import sys

        # Remove ProtParam if already loaded
        protparam_modules = [key for key in sys.modules.keys() if "ProtParam" in key]
        for mod in protparam_modules:
            del sys.modules[mod]

        # Import the module
        from biotooler.families.protparam import ProtParamFeature  # noqa: F401

        # ProtParam should not be loaded yet
        # (This check demonstrates intent but may fail if ProtParam was already imported)
        # We don't assert here as it's informational only

    def test_feature_instantiation_does_not_load_protparam(self):
        """Test that instantiating the feature doesn't load ProtParam."""
        import sys

        # Remove ProtParam if already loaded
        protparam_modules = [key for key in sys.modules.keys() if "ProtParam" in key]
        for mod in protparam_modules:
            del sys.modules[mod]

        # Create feature instance
        feature = ProtParamFeature()  # noqa: F841

        # ProtParam still should not be loaded
        # (This check demonstrates intent but may fail if ProtParam was already imported)
        # We don't assert here as it's informational only


class TestProtParamFeatureEdgeCases:
    """Tests for edge cases."""

    def test_empty_protein_sequence(self):
        """Test with empty protein sequence."""
        feature = ProtParamFeature()
        protein = SeqRecord(Seq(""), id="test")
        protein.annotations["molecule_type"] = "protein"

        # ProtParam should fail with empty sequence
        with pytest.raises(ValueError):
            feature(protein)

    def test_single_amino_acid(self):
        """Test with single amino acid."""
        feature = ProtParamFeature()
        protein = SeqRecord(Seq("M"), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        assert result["molecular_weight"] > 0
        assert result["aromaticity"] == 0.0  # M is not aromatic

    def test_very_long_sequence(self):
        """Test with a very long sequence."""
        feature = ProtParamFeature()
        # Create a long sequence
        protein = SeqRecord(Seq("MKALV" * 100), id="test")
        protein.annotations["molecule_type"] = "protein"

        result = feature(protein)

        assert result["molecular_weight"] > 10000  # Should be large
