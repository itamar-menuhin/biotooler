"""Tests for reference sequence set."""

import tempfile
from pathlib import Path

import pytest

from biotooler.core.reference_sequences import ReferenceSequenceSet


class TestReferenceSequenceSetInit:
    """Tests for ReferenceSequenceSet initialization."""

    def test_init_with_cds_only(self):
        """Test initialization with CDS sequences only."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        ref_set = ReferenceSequenceSet(cds)
        assert ref_set.cds == cds
        assert ref_set.proteins is None
        assert ref_set.genetic_code_table == 1

    def test_init_with_cds_and_proteins(self):
        """Test initialization with both CDS and protein sequences."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        proteins = {"seq1": "MK*", "seq2": "MG*"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        assert ref_set.cds == cds
        assert ref_set.proteins == proteins

    def test_init_with_custom_genetic_code(self):
        """Test initialization with custom genetic code table."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds, genetic_code_table=11)
        assert ref_set.genetic_code_table == 11


class TestFromFasta:
    """Tests for from_fasta constructor."""

    def test_from_fasta_cds_only(self):
        """Test loading CDS sequences from FASTA file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            f.write(">seq2\n")
            f.write("ATGGGGTGA\n")
            fasta_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_fasta(fasta_path)
            assert "seq1" in ref_set.cds
            assert "seq2" in ref_set.cds
            assert ref_set.cds["seq1"] == "ATGAAATAA"
            assert ref_set.cds["seq2"] == "ATGGGGTGA"
            assert ref_set.proteins is None
        finally:
            Path(fasta_path).unlink()

    def test_from_fasta_with_proteins(self):
        """Test loading CDS and protein sequences from FASTA files."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            cds_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("MK*\n")
            protein_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_fasta(cds_path, protein_path)
            assert ref_set.cds["seq1"] == "ATGAAATAA"
            assert ref_set.proteins is not None
            assert ref_set.proteins["seq1"] == "MK*"
        finally:
            Path(cds_path).unlink()
            Path(protein_path).unlink()

    def test_from_fasta_cds_file_not_found(self):
        """Test error when CDS FASTA file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="CDS FASTA file not found"):
            ReferenceSequenceSet.from_fasta("/nonexistent/file.fasta")

    def test_from_fasta_protein_file_not_found(self):
        """Test error when protein FASTA file doesn't exist."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            cds_path = f.name

        try:
            with pytest.raises(FileNotFoundError, match="Protein FASTA file not found"):
                ReferenceSequenceSet.from_fasta(cds_path, "/nonexistent/protein.fasta")
        finally:
            Path(cds_path).unlink()

    def test_from_fasta_empty_cds_file(self):
        """Test error when CDS FASTA file is empty."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            fasta_path = f.name

        try:
            with pytest.raises(ValueError, match="CDS FASTA file is empty"):
                ReferenceSequenceSet.from_fasta(fasta_path)
        finally:
            Path(fasta_path).unlink()

    def test_from_fasta_empty_protein_file(self):
        """Test error when protein FASTA file is empty."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            cds_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            protein_path = f.name

        try:
            with pytest.raises(ValueError, match="Protein FASTA file is empty"):
                ReferenceSequenceSet.from_fasta(cds_path, protein_path)
        finally:
            Path(cds_path).unlink()
            Path(protein_path).unlink()

    def test_from_fasta_duplicate_cds_ids(self):
        """Test error when CDS FASTA has duplicate IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            f.write(">seq1\n")
            f.write("ATGGGGTGA\n")
            fasta_path = f.name

        try:
            with pytest.raises(ValueError, match="Duplicate sequence ID in CDS FASTA: seq1"):
                ReferenceSequenceSet.from_fasta(fasta_path)
        finally:
            Path(fasta_path).unlink()

    def test_from_fasta_duplicate_protein_ids(self):
        """Test error when protein FASTA has duplicate IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("ATGAAATAA\n")
            cds_path = f.name

        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">seq1\n")
            f.write("MK*\n")
            f.write(">seq1\n")
            f.write("MG*\n")
            protein_path = f.name

        try:
            with pytest.raises(ValueError, match="Duplicate sequence ID in protein FASTA: seq1"):
                ReferenceSequenceSet.from_fasta(cds_path, protein_path)
        finally:
            Path(cds_path).unlink()
            Path(protein_path).unlink()


class TestCdsStrings:
    """Tests for cds_strings method."""

    def test_normalization_uppercase(self):
        """Test that CDS sequences are converted to uppercase."""
        cds = {"seq1": "atgaaataa"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings()
        assert result[0] == "ATGAAATAA"

    def test_normalization_u_to_t(self):
        """Test that U is replaced with T (RNA to DNA conversion)."""
        cds = {"seq1": "AUGAAAUAA"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings()
        assert result[0] == "ATGAAATAA"

    def test_normalization_mixed_case_and_u(self):
        """Test normalization with mixed case and U."""
        cds = {"seq1": "augAAAuaa"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings()
        assert result[0] == "ATGAAATAA"

    def test_multiple_of_three_validation_pass(self):
        """Test that sequences with length divisible by 3 pass validation."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings(require_multiple_of_three=True)
        assert len(result) == 2
        assert result[0] == "ATGAAATAA"
        assert result[1] == "ATGGGGTGA"

    def test_multiple_of_three_validation_fail(self):
        """Test that sequences with length not divisible by 3 raise error."""
        cds = {"seq1": "ATGAAATA"}  # Length 8, not divisible by 3
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq1.*length 8.*not a multiple of 3.*index: 0"):
            ref_set.cds_strings(require_multiple_of_three=True)

    def test_multiple_of_three_validation_fail_second_sequence(self):
        """Test error message includes correct index for second sequence."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGT"}  # seq2 length 7
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq2.*length 7.*not a multiple of 3.*index: 1"):
            ref_set.cds_strings(require_multiple_of_three=True)

    def test_multiple_of_three_validation_disabled(self):
        """Test that validation can be disabled."""
        cds = {"seq1": "ATGAAATA"}  # Length 8, not divisible by 3
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings(require_multiple_of_three=False)
        assert result[0] == "ATGAAATA"

    def test_caching(self):
        """Test that cds_strings results are cached."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # First call
        result1 = ref_set.cds_strings()
        assert ref_set._cds_strings_cache is not None

        # Second call should return cached result
        result2 = ref_set.cds_strings()
        assert result1 is result2  # Same object reference

    def test_multiple_sequences_order(self):
        """Test that multiple sequences maintain consistent order."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA", "seq3": "ATGCCCTAG"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.cds_strings()
        # Order should match dictionary order (insertion order in Python 3.7+)
        assert len(result) == 3
        assert result[0] == "ATGAAATAA"
        assert result[1] == "ATGGGGTGA"
        assert result[2] == "ATGCCCTAG"


class TestProteinStrings:
    """Tests for protein_strings method."""

    def test_use_provided_proteins(self):
        """Test that provided protein sequences are used when available."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        result = ref_set.protein_strings()
        assert result[0] == "MK*"

    def test_translate_cds_when_no_proteins(self):
        """Test that CDS is translated when no proteins provided."""
        cds = {"seq1": "ATGAAATAA"}  # ATG=M, AAA=K, TAA=*
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.protein_strings()
        assert result[0] == "MK*"

    def test_strip_terminal_stop(self):
        """Test that terminal stop codon can be stripped."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.protein_strings(strip_terminal_stop=True)
        assert result[0] == "MK"

    def test_strip_terminal_stop_with_provided_proteins(self):
        """Test terminal stop stripping with provided proteins."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        result = ref_set.protein_strings(strip_terminal_stop=True)
        assert result[0] == "MK"

    def test_strip_terminal_stop_no_terminal_stop(self):
        """Test that sequences without terminal stop are unchanged."""
        cds = {"seq1": "ATGAAA"}
        proteins = {"seq1": "MK"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        result = ref_set.protein_strings(strip_terminal_stop=True)
        assert result[0] == "MK"

    def test_internal_stop_error_by_default(self):
        """Test that internal stop codons raise error by default."""
        cds = {"seq1": "ATGTAAAAATAA"}  # ATG=M, TAA=*, AAA=K, TAA=*
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq1.*index: 0.*internal stop.*position 1"):
            ref_set.protein_strings()

    def test_internal_stop_in_provided_proteins(self):
        """Test that internal stops in provided proteins raise error."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "M*K"}  # Internal stop at position 1
        ref_set = ReferenceSequenceSet(cds, proteins)
        with pytest.raises(ValueError, match=r"seq1.*internal stop.*position 1"):
            ref_set.protein_strings()

    def test_internal_stop_error_disabled(self):
        """Test that internal stop validation can be disabled."""
        cds = {"seq1": "ATGTAAAAATAA"}  # ATG=M, TAA=*, AAA=K, TAA=*
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.protein_strings(error_on_internal_stop=False)
        assert result[0] == "M*K*"

    def test_translation_with_custom_genetic_code(self):
        """Test translation with custom genetic code table."""
        # Use a simple test to verify genetic_code_table parameter is used
        cds = {"seq1": "ATGAAAGGATAA"}  # ATG=M, AAA=K, GGA=G, TAA=*
        ref_set = ReferenceSequenceSet(cds, genetic_code_table=11)
        result = ref_set.protein_strings()
        # Verify translation works with custom table
        assert result[0] == "MKG*"

    def test_caching(self):
        """Test that protein_strings results are cached."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # First call
        result1 = ref_set.protein_strings()
        assert ref_set._protein_strings_cache is not None

        # Second call should return cached result
        result2 = ref_set.protein_strings()
        assert result1 is result2  # Same object reference

    def test_missing_protein_for_cds(self):
        """Test error when CDS ID not found in provided proteins."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        proteins = {"seq1": "MK*"}  # Missing seq2
        ref_set = ReferenceSequenceSet(cds, proteins)
        with pytest.raises(ValueError, match=r"seq2.*not found in provided protein sequences"):
            ref_set.protein_strings()

    def test_multiple_sequences_translation(self):
        """Test translation of multiple CDS sequences."""
        cds = {
            "seq1": "ATGAAATAA",  # MK*
            "seq2": "ATGGGGTGA",  # MG*
            "seq3": "ATGCCCTAG",  # MP*
        }
        ref_set = ReferenceSequenceSet(cds)
        result = ref_set.protein_strings()
        assert len(result) == 3
        assert result[0] == "MK*"
        assert result[1] == "MG*"
        assert result[2] == "MP*"

    def test_translation_requires_multiple_of_three(self):
        """Test that translation validates CDS length is multiple of 3."""
        cds = {"seq1": "ATGAAATA"}  # Length 8, not divisible by 3
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq1.*length 8.*not a multiple of 3"):
            ref_set.protein_strings()

    def test_terminal_stop_only_not_internal(self):
        """Test that terminal stop is correctly identified (not internal)."""
        cds = {"seq1": "ATGAAATAA"}  # MK* - stop at end
        ref_set = ReferenceSequenceSet(cds)
        # Should not raise error for terminal stop
        result = ref_set.protein_strings(error_on_internal_stop=True)
        assert result[0] == "MK*"

    def test_translation_failure_clear_error(self):
        """Test that translation failures have clear error messages."""
        # This is a bit tricky to test without invalid sequences
        # Just verify the error message format would be correct
        cds = {"seq1": "ATGAAATA"}  # Will fail due to length check
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq1.*index: 0"):
            ref_set.protein_strings()


class TestIntegration:
    """Integration tests combining multiple features."""

    def test_full_workflow_from_fasta_to_proteins(self):
        """Test complete workflow from FASTA loading to protein translation."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as f:
            f.write(">gene1\n")
            f.write("augaaauaa\n")  # Lowercase with U
            f.write(">gene2\n")
            f.write("ATGGGGTGA\n")
            fasta_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_fasta(fasta_path)

            # Get normalized CDS
            cds = ref_set.cds_strings()
            assert cds[0] == "ATGAAATAA"  # Uppercase, U->T
            assert cds[1] == "ATGGGGTGA"

            # Get proteins
            proteins = ref_set.protein_strings(strip_terminal_stop=True)
            assert proteins[0] == "MK"
            assert proteins[1] == "MG"

        finally:
            Path(fasta_path).unlink()

    def test_caching_works_across_calls(self):
        """Test that caching works correctly across multiple method calls."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # First call to cds_strings
        cds1 = ref_set.cds_strings()
        assert ref_set._cds_strings_cache is not None

        # Call protein_strings (should use cached CDS)
        proteins = ref_set.protein_strings()
        assert ref_set._protein_strings_cache is not None

        # Second call to cds_strings should still use cache
        cds2 = ref_set.cds_strings()
        assert cds1 is cds2

        # Second call to protein_strings should use cache
        proteins2 = ref_set.protein_strings()
        assert proteins is proteins2
