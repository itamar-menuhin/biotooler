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


class TestFromCsv:
    """Tests for from_csv constructor."""

    def test_from_csv_cds_only(self):
        """Test loading CDS sequences from CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence\n")
            f.write("seq1,ATGAAATAA\n")
            f.write("seq2,ATGGGGTGA\n")
            csv_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
            assert "seq1" in ref_set.cds
            assert "seq2" in ref_set.cds
            assert ref_set.cds["seq1"] == "ATGAAATAA"
            assert ref_set.cds["seq2"] == "ATGGGGTGA"
            assert ref_set.proteins is None
        finally:
            Path(csv_path).unlink()

    def test_from_csv_with_proteins(self):
        """Test loading CDS and protein sequences from CSV file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,cds,protein\n")
            f.write("seq1,ATGAAATAA,MK*\n")
            f.write("seq2,ATGGGGTGA,MG*\n")
            csv_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_csv(csv_path, "id", "cds", protein_column="protein")
            assert ref_set.cds["seq1"] == "ATGAAATAA"
            assert ref_set.proteins is not None
            assert ref_set.proteins["seq1"] == "MK*"
            assert ref_set.proteins["seq2"] == "MG*"
        finally:
            Path(csv_path).unlink()

    def test_from_csv_file_not_found(self):
        """Test error when CSV file doesn't exist."""
        with pytest.raises(FileNotFoundError, match="CSV file not found"):
            ReferenceSequenceSet.from_csv("/nonexistent/file.csv", "id", "sequence")

    def test_from_csv_empty_file(self):
        """Test error when CSV file is empty."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="Failed to read CSV file"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_missing_id_column(self):
        """Test error when ID column is missing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("sequence\n")
            f.write("ATGAAATAA\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="ID column 'id' not found"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_missing_cds_column(self):
        """Test error when CDS column is missing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id\n")
            f.write("seq1\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="CDS column 'sequence' not found"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_missing_protein_column(self):
        """Test error when specified protein column is missing."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence\n")
            f.write("seq1,ATGAAATAA\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="Protein column 'protein' not found"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence", "protein")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_duplicate_ids(self):
        """Test error when CSV has duplicate IDs."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence\n")
            f.write("seq1,ATGAAATAA\n")
            f.write("seq1,ATGGGGTGA\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="Duplicate sequence ID in CSV: seq1"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_missing_cds_value(self):
        """Test error when CDS value is missing for a row."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence\n")
            f.write("seq1,ATGAAATAA\n")
            f.write("seq2,\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="CDS sequence is missing for ID 'seq2'"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_missing_protein_value(self):
        """Test error when protein value is missing for a row."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence,protein\n")
            f.write("seq1,ATGAAATAA,MK*\n")
            f.write("seq2,ATGGGGTGA,\n")
            csv_path = f.name

        try:
            with pytest.raises(ValueError, match="Protein sequence is missing for ID 'seq2'"):
                ReferenceSequenceSet.from_csv(csv_path, "id", "sequence", "protein")
        finally:
            Path(csv_path).unlink()

    def test_from_csv_with_custom_genetic_code(self):
        """Test loading with custom genetic code table."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("id,sequence\n")
            f.write("seq1,ATGAAATAA\n")
            csv_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_csv(
                csv_path, "id", "sequence", genetic_code_table=11
            )
            assert ref_set.genetic_code_table == 11
        finally:
            Path(csv_path).unlink()

    def test_from_csv_different_column_names(self):
        """Test using different column names."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("gene_id,nucleotide_seq,amino_acid_seq\n")
            f.write("gene1,ATGAAATAA,MK*\n")
            csv_path = f.name

        try:
            ref_set = ReferenceSequenceSet.from_csv(
                csv_path, "gene_id", "nucleotide_seq", "amino_acid_seq"
            )
            assert ref_set.cds["gene1"] == "ATGAAATAA"
            assert ref_set.proteins is not None
            assert ref_set.proteins["gene1"] == "MK*"
        finally:
            Path(csv_path).unlink()


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


class TestValidate:
    """Tests for validate method."""

    def test_validate_cds_success(self):
        """Test successful CDS validation."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        ref_set = ReferenceSequenceSet(cds)
        # Should not raise
        ref_set.validate(kind="cds")

    def test_validate_cds_empty_sequence(self):
        """Test that empty CDS sequences raise error."""
        cds = {"seq1": "ATGAAATAA", "seq2": ""}
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match="CDS sequence 'seq2' is empty"):
            ref_set.validate(kind="cds")

    def test_validate_cds_whitespace_only(self):
        """Test that whitespace-only CDS sequences raise error."""
        cds = {"seq1": "ATGAAATAA", "seq2": "   "}
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match="CDS sequence 'seq2' is empty"):
            ref_set.validate(kind="cds")

    def test_validate_cds_not_multiple_of_three(self):
        """Test that CDS not multiple of 3 raises error."""
        cds = {"seq1": "ATGAAATA"}  # Length 8
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match=r"seq1.*length 8.*not a multiple of 3"):
            ref_set.validate(kind="cds")

    def test_validate_protein_success(self):
        """Test successful protein validation."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        # Should not raise
        ref_set.validate(kind="protein")

    def test_validate_protein_empty_sequence(self):
        """Test that empty protein sequences raise error."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        proteins = {"seq1": "MK*", "seq2": ""}
        ref_set = ReferenceSequenceSet(cds, proteins)
        with pytest.raises(ValueError, match="Protein sequence 'seq2' is empty"):
            ref_set.validate(kind="protein")

    def test_validate_protein_internal_stop(self):
        """Test that internal stops in proteins raise error."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "M*K"}  # Internal stop
        ref_set = ReferenceSequenceSet(cds, proteins)
        with pytest.raises(ValueError, match=r"seq1.*internal stop.*position 1"):
            ref_set.validate(kind="protein")

    def test_validate_both_success(self):
        """Test successful validation of both CDS and proteins."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        proteins = {"seq1": "MK*", "seq2": "MG*"}
        ref_set = ReferenceSequenceSet(cds, proteins)
        # Should not raise
        ref_set.validate(kind="both")

    def test_validate_both_default(self):
        """Test that default kind is 'both'."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)
        # Should not raise (validates both CDS and translated proteins)
        ref_set.validate()

    def test_validate_invalid_kind(self):
        """Test that invalid kind raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)
        with pytest.raises(ValueError, match="Invalid kind 'invalid'"):
            ref_set.validate(kind="invalid")


class TestWithProteins:
    """Tests for with_proteins method."""

    def test_with_proteins_add_new(self):
        """Test adding new proteins to CDS-only set."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        ref_set = ReferenceSequenceSet(cds)

        new_proteins = {"seq1": "MK*", "seq2": "MG*"}
        result = ref_set.with_proteins(new_proteins)

        # Should return new instance
        assert result is not ref_set
        assert result.proteins == new_proteins
        assert result.cds == cds
        assert result.genetic_code_table == ref_set.genetic_code_table

    def test_with_proteins_merge_with_existing(self):
        """Test merging new proteins with existing proteins."""
        cds = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA", "seq3": "ATGCCCTAG"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)

        new_proteins = {"seq2": "MG*", "seq3": "MP*"}
        result = ref_set.with_proteins(new_proteins)

        assert result.proteins == {"seq1": "MK*", "seq2": "MG*", "seq3": "MP*"}

    def test_with_proteins_duplicate_same_sequence(self):
        """Test that providing same sequence for existing ID is allowed."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)

        # Should not raise
        result = ref_set.with_proteins({"seq1": "MK*"})
        assert result.proteins == proteins

    def test_with_proteins_duplicate_different_sequence(self):
        """Test that different sequence for existing ID raises error."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)

        with pytest.raises(
            ValueError,
            match=r"Protein sequence 'seq1' already exists with a different sequence",
        ):
            ref_set.with_proteins({"seq1": "MKK*"})

    def test_with_proteins_empty_sequence(self):
        """Test that empty protein sequence raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        with pytest.raises(ValueError, match="Provided protein sequence 'seq1' is empty"):
            ref_set.with_proteins({"seq1": ""})

    def test_with_proteins_whitespace_only(self):
        """Test that whitespace-only protein sequence raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        with pytest.raises(ValueError, match="Provided protein sequence 'seq1' is empty"):
            ref_set.with_proteins({"seq1": "   "})

    def test_with_proteins_original_unchanged(self):
        """Test that original set is unchanged."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        result = ref_set.with_proteins({"seq1": "MK*"})

        # Original should be unchanged
        assert ref_set.proteins is None
        assert result.proteins is not None


class TestExtendCds:
    """Tests for extend_cds method."""

    def test_extend_cds_add_new(self):
        """Test adding new CDS sequences."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        new_cds = {"seq2": "ATGGGGTGA", "seq3": "ATGCCCTAG"}
        result = ref_set.extend_cds(new_cds)

        # Should return new instance
        assert result is not ref_set
        assert result.cds == {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA", "seq3": "ATGCCCTAG"}
        assert result.proteins == ref_set.proteins

    def test_extend_cds_duplicate_same_sequence(self):
        """Test that providing same sequence for existing ID is allowed."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # Should not raise (same sequence)
        result = ref_set.extend_cds({"seq1": "ATGAAATAA"})
        assert result.cds == cds

    def test_extend_cds_duplicate_same_sequence_normalized(self):
        """Test that normalized sequences are considered the same."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # Should not raise (lowercase and U->T normalization)
        result = ref_set.extend_cds({"seq1": "augaaauaa"})
        assert result.cds == {"seq1": "ATGAAATAA"}

    def test_extend_cds_duplicate_different_sequence(self):
        """Test that different sequence for existing ID raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        with pytest.raises(
            ValueError, match=r"CDS sequence 'seq1' already exists with a different sequence"
        ):
            ref_set.extend_cds({"seq1": "ATGGGGTGA"})

    def test_extend_cds_empty_sequence(self):
        """Test that empty CDS sequence raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        with pytest.raises(ValueError, match="Provided CDS sequence 'seq2' is empty"):
            ref_set.extend_cds({"seq2": ""})

    def test_extend_cds_whitespace_only(self):
        """Test that whitespace-only CDS sequence raises error."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        with pytest.raises(ValueError, match="Provided CDS sequence 'seq2' is empty"):
            ref_set.extend_cds({"seq2": "   "})

    def test_extend_cds_preserves_proteins(self):
        """Test that extending CDS preserves protein sequences."""
        cds = {"seq1": "ATGAAATAA"}
        proteins = {"seq1": "MK*"}
        ref_set = ReferenceSequenceSet(cds, proteins)

        result = ref_set.extend_cds({"seq2": "ATGGGGTGA"})

        assert result.proteins == proteins

    def test_extend_cds_cache_invalidated(self):
        """Test that cache is invalidated in new instance."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        # Populate cache
        ref_set.cds_strings()

        # Extend
        result = ref_set.extend_cds({"seq2": "ATGGGGTGA"})

        # New instance should have empty cache
        assert result._cds_strings_cache is None

    def test_extend_cds_original_unchanged(self):
        """Test that original set is unchanged."""
        cds = {"seq1": "ATGAAATAA"}
        ref_set = ReferenceSequenceSet(cds)

        result = ref_set.extend_cds({"seq2": "ATGGGGTGA"})

        # Original should be unchanged
        assert ref_set.cds == {"seq1": "ATGAAATAA"}
        assert result.cds == {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}


class TestMerge:
    """Tests for merge method."""

    def test_merge_disjoint_sets(self):
        """Test merging sets with no overlapping IDs."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq2": "ATGGGGTGA"}
        ref_set2 = ReferenceSequenceSet(cds2)

        result = ref_set1.merge(ref_set2)

        assert result.cds == {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        assert result is not ref_set1
        assert result is not ref_set2

    def test_merge_preserves_order(self):
        """Test that merge preserves order with this instance first."""
        cds1 = {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq3": "ATGCCCTAG"}
        ref_set2 = ReferenceSequenceSet(cds2)

        result = ref_set1.merge(ref_set2)

        # Order should be seq1, seq2, seq3 (this instance first)
        assert list(result.cds.keys()) == ["seq1", "seq2", "seq3"]

    def test_merge_duplicate_same_sequence(self):
        """Test that merging with same sequence for same ID is allowed."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq1": "ATGAAATAA"}
        ref_set2 = ReferenceSequenceSet(cds2)

        # Should not raise
        result = ref_set1.merge(ref_set2)
        assert result.cds == {"seq1": "ATGAAATAA"}

    def test_merge_duplicate_same_sequence_normalized(self):
        """Test that normalized sequences are considered the same."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq1": "augaaauaa"}
        ref_set2 = ReferenceSequenceSet(cds2)

        # Should not raise
        result = ref_set1.merge(ref_set2)
        assert result.cds == {"seq1": "ATGAAATAA"}

    def test_merge_duplicate_different_sequence(self):
        """Test that different sequences for same ID raise error."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq1": "ATGGGGTGA"}
        ref_set2 = ReferenceSequenceSet(cds2)

        with pytest.raises(
            ValueError, match=r"CDS sequence 'seq1' exists in both sets with different sequences"
        ):
            ref_set1.merge(ref_set2)

    def test_merge_with_proteins(self):
        """Test merging sets that both have proteins."""
        cds1 = {"seq1": "ATGAAATAA"}
        proteins1 = {"seq1": "MK*"}
        ref_set1 = ReferenceSequenceSet(cds1, proteins1)

        cds2 = {"seq2": "ATGGGGTGA"}
        proteins2 = {"seq2": "MG*"}
        ref_set2 = ReferenceSequenceSet(cds2, proteins2)

        result = ref_set1.merge(ref_set2)

        assert result.proteins == {"seq1": "MK*", "seq2": "MG*"}

    def test_merge_mixed_proteins(self):
        """Test merging when only one set has proteins."""
        cds1 = {"seq1": "ATGAAATAA"}
        proteins1 = {"seq1": "MK*"}
        ref_set1 = ReferenceSequenceSet(cds1, proteins1)

        cds2 = {"seq2": "ATGGGGTGA"}
        ref_set2 = ReferenceSequenceSet(cds2)

        result = ref_set1.merge(ref_set2)

        # Should preserve proteins from first set
        assert result.proteins == {"seq1": "MK*"}

    def test_merge_protein_conflict(self):
        """Test that conflicting proteins raise error."""
        cds1 = {"seq1": "ATGAAATAA"}
        proteins1 = {"seq1": "MK*"}
        ref_set1 = ReferenceSequenceSet(cds1, proteins1)

        cds2 = {"seq1": "ATGAAATAA"}  # Same CDS
        proteins2 = {"seq1": "MKK*"}  # Different protein
        ref_set2 = ReferenceSequenceSet(cds2, proteins2)

        with pytest.raises(
            ValueError,
            match=r"Protein sequence 'seq1' exists in both sets with different sequences",
        ):
            ref_set1.merge(ref_set2)

    def test_merge_preserves_genetic_code(self):
        """Test that merge preserves genetic code from first set."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1, genetic_code_table=11)

        cds2 = {"seq2": "ATGGGGTGA"}
        ref_set2 = ReferenceSequenceSet(cds2, genetic_code_table=1)

        result = ref_set1.merge(ref_set2)

        assert result.genetic_code_table == 11

    def test_merge_empty_cds_in_other(self):
        """Test that empty CDS in other set raises error."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq2": ""}
        ref_set2 = ReferenceSequenceSet(cds2)

        with pytest.raises(ValueError, match="CDS sequence 'seq2' in other set is empty"):
            ref_set1.merge(ref_set2)

    def test_merge_empty_protein_in_other(self):
        """Test that empty protein in other set raises error."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq2": "ATGGGGTGA"}
        proteins2 = {"seq2": ""}
        ref_set2 = ReferenceSequenceSet(cds2, proteins2)

        with pytest.raises(ValueError, match="Protein sequence 'seq2' in other set is empty"):
            ref_set1.merge(ref_set2)

    def test_merge_original_unchanged(self):
        """Test that original sets are unchanged."""
        cds1 = {"seq1": "ATGAAATAA"}
        ref_set1 = ReferenceSequenceSet(cds1)

        cds2 = {"seq2": "ATGGGGTGA"}
        ref_set2 = ReferenceSequenceSet(cds2)

        result = ref_set1.merge(ref_set2)

        # Originals should be unchanged
        assert ref_set1.cds == {"seq1": "ATGAAATAA"}
        assert ref_set2.cds == {"seq2": "ATGGGGTGA"}
        assert result.cds == {"seq1": "ATGAAATAA", "seq2": "ATGGGGTGA"}
