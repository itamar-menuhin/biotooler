"""Tests for ORF policy feature in ensure_protein_record."""

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_store import attach_orf
from biotooler.core.translation import ensure_protein_record


class TestOrfPolicyDefault:
    """Tests for default ORF policy behavior."""

    def test_default_policy_uses_frame0(self):
        """Test that default policy uses frame 0 for full sequence."""
        # Sequence: NNN ATG AAA GCC CTG GTG TAA
        # Frame 0: XMKALV*
        dna = SeqRecord(Seq("NNNATGAAAGCCCTGGTGTAA"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="default")

        # Default policy translates from frame 0
        assert result.annotations["translation_region"] == (0, 21)
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"

    def test_default_policy_is_truly_default(self):
        """Test that omitting orf_policy is same as orf_policy='default'."""
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result_default = ensure_protein_record(dna)
        result_explicit = ensure_protein_record(dna, orf_policy="default")

        assert str(result_default.seq) == str(result_explicit.seq)
        assert (
            result_default.annotations["translation_region"]
            == result_explicit.annotations["translation_region"]
        )
        assert (
            result_default.annotations["translation_region_source"]
            == result_explicit.annotations["translation_region_source"]
        )

    def test_invalid_orf_policy_raises_error(self):
        """Test that invalid orf_policy raises ValueError."""
        dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        dna.annotations["molecule_type"] = "DNA"

        with pytest.raises(ValueError, match="orf_policy must be"):
            ensure_protein_record(dna, orf_policy="invalid")


class TestOrfPolicyLongestOrf:
    """Tests for longest_orf policy."""

    def test_longest_orf_simple_case(self):
        """Test longest_orf policy selects the longest ORF."""
        # Sequence with multiple ORFs:
        # ATG AAA TAA ATG CCC TAG TAA
        # 0   3   6   9   12  15  18  21
        # ORF 1: (0, 9) = 9 bases (ATGAAATAA)
        # ORF 2: (0, 18) = 18 bases (ATGAAATAGATGCCCTAG)
        # ORF 3: (0, 21) = 21 bases (ATGAAATAGATGCCCTAGTAA)  <- longest
        # ORF 4: (9, 18) = 9 bases (ATGCCCTAG)
        # ORF 5: (9, 21) = 12 bases (ATGCCCTAGTAA)
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        # Should select the longest ORF: (0, 21)
        assert result.annotations["translation_region"] == (0, 21)
        assert result.annotations["translation_region_source"] == "longest_orf"
        # Translates to: MK*MP** (ATGAAATAGATGCCCTAGTAA has TAG and TAA both as stops)
        # With strip_terminal_stop=True: MK*MP* (only the last stop is stripped)
        assert str(result.seq) == "MK*MP*"

    def test_longest_orf_with_internal_stops_ignored(self):
        """Test longest_orf policy with internal stops when allowed."""
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        assert result.annotations["translation_region"] == (0, 21)
        assert result.annotations["translation_region_source"] == "longest_orf"

    def test_longest_orf_with_prefix(self):
        """Test longest_orf ignores non-ORF prefix sequence."""
        # Prefix NNN, then ORFs starting at position 3
        # NNN ATG AAA TAA
        # 0   3   6   9   12
        # Only ORF: (3, 12) = 9 bases
        seq = "NNNATGAAATAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        assert result.annotations["translation_region"] == (3, 12)
        assert result.annotations["translation_region_source"] == "longest_orf"
        assert str(result.seq) == "MK"

    def test_longest_orf_rna_sequence(self):
        """Test longest_orf policy works with RNA sequences."""
        # AUG AAA UAA AUG CCC UAG
        # 0   3   6   9   12  15  18
        # ORF 1: (0, 9) = 9 bases
        # ORF 2: (0, 18) = 18 bases  <- longest
        # ORF 3: (9, 18) = 9 bases
        seq = "AUGAAAUAAAUGCCCUAG"
        rna = SeqRecord(Seq(seq), id="test")
        rna.annotations["molecule_type"] = "RNA"

        result = ensure_protein_record(rna, orf_policy="longest_orf", on_internal_stop="ignore")

        assert result.annotations["translation_region"] == (0, 18)
        assert result.annotations["translation_region_source"] == "longest_orf"


class TestOrfPolicyTieBreaking:
    """Tests for deterministic tie-breaking in longest_orf policy."""

    def test_tie_breaking_earliest_start(self):
        """Test that ties are broken by earliest start position."""
        # Create two ORFs of same length at different positions:
        # ATG AAA TAA GGG AAA ATG CCC TAG
        # 0   3   6   9   12  15  18  21  24
        # ORF 1: (0, 9) = 9 bases   <- earliest start, wins tie
        # ORF 2: (15, 24) = 9 bases
        seq = "ATGAAATAAGGAAAATGCCCTAG"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        # Should select earliest start when lengths are equal
        assert result.annotations["translation_region"] == (0, 9)
        assert result.annotations["translation_region_source"] == "longest_orf"
        assert str(result.seq) == "MK"

    def test_tie_breaking_earliest_stop(self):
        """Test that ties are broken by earliest stop when starts are same."""
        # Create ORFs with same start but different stops:
        # ATG AAA TAG ATG CCC TAG TAA
        # 0   3   6   9   12  15  18  21
        # From start 0:
        # ORF 1: (0, 9) = 9 bases   <- earliest stop among same length from 0
        # ORF 2: (0, 18) = 18 bases
        # ORF 3: (0, 21) = 21 bases  <- longest overall, should win
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        # Should select the longest, which is (0, 21)
        assert result.annotations["translation_region"] == (0, 21)

    def test_tie_breaking_multiple_same_length_different_starts(self):
        """Test longest ORF selection when multiple starts create overlapping ORFs."""
        # When a sequence has multiple start codons and stops, the first start
        # can pair with multiple stops, creating ORFs of different lengths
        # ATG AAA TAA GGG AAA AAA ATG CCC TAG AAA AAA AAA ATG GGG TAA
        # From start at 0: Can reach stops at 6, 18, 30, 42
        # From start at 18: Can reach stops at 24, 42
        # From start at 36: Can reach stop at 42
        # The longest ORF is (0, 42) = 42 bases
        seq = "ATGAAATAAGGAAAAAATGCCCTAGAAAAAAAAATGGGGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        # Should select the longest ORF: (0, 42)
        assert result.annotations["translation_region"] == (0, 42)
        assert result.annotations["translation_region_source"] == "longest_orf"

    def test_tie_breaking_same_start_and_length(self):
        """Test tie-breaking when ORFs have same start and length."""
        # This is a pathological case, but we should handle it deterministically
        # ATG AAA TAG TAA (TAG and TAA both are stops)
        # 0   3   6   9
        # ORF 1: (0, 9) = 9 bases with TAG
        # ORF 2: (0, 12) = 12 bases with TAA
        # Different lengths, so (0, 12) should win
        seq = "ATGAAATAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        # Should select the longest: (0, 12)
        assert result.annotations["translation_region"] == (0, 12)


class TestOrfPolicyFallback:
    """Tests for fallback behavior when no ORF candidates found."""

    def test_no_orf_candidates_fallback_to_frame0(self):
        """Test that longest_orf falls back to frame0 when no ORFs found."""
        # Sequence with no start codons
        seq = "AAACCCGGGTTT"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        # Should fall back to frame 0
        assert result.annotations["translation_region"] == (0, 12)
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"

    def test_no_stop_codons_fallback_to_frame0(self):
        """Test fallback when start codon exists but no stop codons."""
        # Sequence with start but no stops
        seq = "ATGAAACCCGGG"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        # Should fall back to frame 0 (no valid ORF candidates)
        assert result.annotations["translation_region"] == (0, 12)
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"

    def test_empty_sequence_with_longest_orf(self):
        """Test longest_orf policy with empty sequence."""
        dna = SeqRecord(Seq(""), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        assert str(result.seq) == ""
        assert result.annotations["translation_region"] == (0, 0)
        assert result.annotations["translation_region_source"] == "full_sequence_frame0"


class TestOrfPolicyPriorityInteraction:
    """Tests for interaction between orf_policy and explicit/attached ORFs."""

    def test_explicit_orf_overrides_longest_orf_policy(self):
        """Test that explicit ORF takes priority over longest_orf policy."""
        # Sequence with multiple ORFs
        # ATG AAA TAA ATG CCC TAG TAA
        # Longest ORF would be (0, 21), but we explicitly use (9, 18)
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf=(9, 18), orf_policy="longest_orf")

        # Explicit ORF should win
        assert result.annotations["translation_region"] == (9, 18)
        assert result.annotations["translation_region_source"] == "explicit_orf"

    def test_attached_orf_overrides_longest_orf_policy(self):
        """Test that attached ORF takes priority over longest_orf policy."""
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (9, 18))

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        # Attached ORF should win over longest_orf policy
        assert result.annotations["translation_region"] == (9, 18)
        assert result.annotations["translation_region_source"] == "attached_orf"

    def test_longest_orf_with_use_orf_if_present_false(self):
        """Test longest_orf policy when use_orf_if_present=False."""
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (9, 18))  # Attach an ORF that's not the longest

        result = ensure_protein_record(
            dna, use_orf_if_present=False, orf_policy="longest_orf", on_internal_stop="ignore"
        )

        # Should use longest_orf policy and ignore attached ORF
        assert result.annotations["translation_region"] == (0, 21)
        assert result.annotations["translation_region_source"] == "longest_orf"

    def test_explicit_orf_overrides_attached_orf_with_longest_orf_policy(self):
        """Test that explicit ORF beats attached ORF even with longest_orf policy."""
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"
        dna = attach_orf(dna, (9, 18))

        # Explicit ORF (0, 9) overrides both attached (9, 18) and longest_orf policy
        result = ensure_protein_record(dna, orf=(0, 9), orf_policy="longest_orf")

        assert result.annotations["translation_region"] == (0, 9)
        assert result.annotations["translation_region_source"] == "explicit_orf"


class TestOrfPolicyProteinInput:
    """Tests for orf_policy with protein input."""

    def test_protein_ignores_orf_policy(self):
        """Test that protein sequences ignore orf_policy parameter."""
        protein = SeqRecord(Seq("MKALV"), id="test")
        protein.annotations["molecule_type"] = "protein"

        # orf_policy should be ignored for protein sequences
        result = ensure_protein_record(protein, orf_policy="longest_orf")

        assert str(result.seq) == "MKALV"
        assert result.annotations["molecule_type"] == "protein"
        # No translation annotations should be present
        assert "translation_performed" not in result.annotations


class TestOrfPolicyRealWorldScenarios:
    """Tests for real-world usage scenarios."""

    def test_longest_orf_with_upstream_sequence(self):
        """Test selecting longest ORF from sequence with 5' UTR."""
        # Simulate a gene with 5' UTR followed by coding sequence
        # NNN NNN ATG AAA GCC CTG GTG TAA NNN
        # Longest ORF: (6, 24) = 18 bases
        seq = "NNNNNNATGAAAGCCCTGGTGTAANNN"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        assert result.annotations["translation_region"] == (6, 24)
        assert result.annotations["translation_region_source"] == "longest_orf"
        assert str(result.seq) == "MKALV"

    def test_longest_orf_selects_best_among_nested_orfs(self):
        """Test that longest_orf correctly handles nested ORFs."""
        # ATG (start1) ... ATG (start2) ... TAA (stop for both)
        # Inner ORF is shorter, outer is longer
        # ATG AAA ATG CCC CCC CCC TAA
        # 0   3   6   9   12  15  18  21
        # ORF from 0: (0, 21) = 21 bases  <- longest
        # ORF from 6: (6, 21) = 15 bases
        seq = "ATGAAAATGCCCCCCCCCTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result = ensure_protein_record(dna, orf_policy="longest_orf")

        # Should select the longest (outermost) ORF
        assert result.annotations["translation_region"] == (0, 21)
        assert result.annotations["translation_region_source"] == "longest_orf"

    def test_longest_orf_deterministic_across_calls(self):
        """Test that longest_orf policy is deterministic."""
        # Same sequence should always give same result
        seq = "ATGAAATAGATGCCCTAGTAA"
        dna = SeqRecord(Seq(seq), id="test")
        dna.annotations["molecule_type"] = "DNA"

        result1 = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")
        result2 = ensure_protein_record(dna, orf_policy="longest_orf", on_internal_stop="ignore")

        assert (
            result1.annotations["translation_region"] == result2.annotations["translation_region"]
        )
        assert str(result1.seq) == str(result2.seq)
