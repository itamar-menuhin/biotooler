"""Tests for ViennaRNA ContextWindowFoldFeature."""

import sys
from unittest.mock import MagicMock, patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


@pytest.fixture
def mock_rna():
    """Create a mock RNA module for testing."""
    # Create mock RNA module
    rna_mock = MagicMock()

    # Mock fold_compound with realistic structure output
    def mock_fold_compound(seq):
        fc_mock = MagicMock()
        # Return a simple structure with some pairs
        # For testing, create a structure with pairs at positions 0-9, 1-8, etc.
        n = len(seq)
        structure = "." * n  # Start with all unpaired
        # Create a simple structure for testing
        if n >= 10:
            # Make pairs: 0-9, 1-8, 2-7, 3-6, 4-5
            structure_list = list(structure)
            pairs = min(5, n // 2)
            for i in range(pairs):
                if i < n and (n - 1 - i) < n:
                    structure_list[i] = "("
                    structure_list[n - 1 - i] = ")"
            structure = "".join(structure_list)

        # MFE based on sequence length (deterministic for testing)
        mfe_value = -len(seq) * 0.5
        fc_mock.mfe.return_value = (structure, mfe_value)
        return fc_mock

    rna_mock.fold_compound = mock_fold_compound
    return rna_mock


class TestContextWindowFoldFeature:
    """Tests for ContextWindowFoldFeature."""

    def test_import_does_not_load_rna(self):
        """Test that importing context_window_fold module doesn't import RNA."""
        # Clear any previously imported modules
        modules_to_clear = [
            module
            for module in list(sys.modules.keys())
            if module.startswith("biotooler.families.viennarna.context_window_fold")
            or module == "RNA"
            or module.startswith("RNA.")
        ]
        for module in modules_to_clear:
            del sys.modules[module]

        # Import the context_window_fold module
        from biotooler.families.viennarna import context_window_fold  # noqa: F401

        # Verify RNA is NOT loaded
        assert "RNA" not in sys.modules, (
            "RNA was imported when importing context_window_fold module, "
            "but it should only be loaded when features are instantiated/called"
        )

    def test_init_validates_parameters(self):
        """Test that __init__ validates parameters."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        # Test empty starts_nt
        with pytest.raises(ValueError, match="starts_nt must not be empty"):
            ContextWindowFoldFeature(
                starts_nt=[],
                window_size_nt=20,
                flank_left_nt=5,
                flank_right_nt=5,
            )

        # Test non-positive window_size_nt
        with pytest.raises(ValueError, match="window_size_nt must be positive"):
            ContextWindowFoldFeature(
                starts_nt=[0, 10],
                window_size_nt=0,
                flank_left_nt=5,
                flank_right_nt=5,
            )

        with pytest.raises(ValueError, match="window_size_nt must be positive"):
            ContextWindowFoldFeature(
                starts_nt=[0, 10],
                window_size_nt=-5,
                flank_left_nt=5,
                flank_right_nt=5,
            )

        # Test negative flank_left_nt
        with pytest.raises(ValueError, match="flank_left_nt must be non-negative"):
            ContextWindowFoldFeature(
                starts_nt=[0, 10],
                window_size_nt=20,
                flank_left_nt=-5,
                flank_right_nt=5,
            )

        # Test negative flank_right_nt
        with pytest.raises(ValueError, match="flank_right_nt must be non-negative"):
            ContextWindowFoldFeature(
                starts_nt=[0, 10],
                window_size_nt=20,
                flank_left_nt=5,
                flank_right_nt=-5,
            )

        # Test invalid mode
        with pytest.raises(ValueError, match="mode must be 'mfe'"):
            ContextWindowFoldFeature(
                starts_nt=[0, 10],
                window_size_nt=20,
                mode="invalid",
            )

    def test_compute_context_mfe_for_requested_windows(self, mock_rna):
        """Test that context MFE is computed only for requested window positions."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        # Patch require_viennarna to return our mock
        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[0, 10, 20],
                window_size_nt=15,
                flank_left_nt=5,
                flank_right_nt=5,
            )
            # Use 40 bp sequence so all windows fit
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")  # 40 bp

            result = feature(record)

            # Check that only requested windows are computed
            assert "CTX_MFE_0" in result
            assert "CTX_MFE_10" in result
            assert "CTX_MFE_20" in result
            assert "PAIR_OUT_FRAC_0" in result
            assert "PAIR_OUT_FRAC_10" in result
            assert "PAIR_OUT_FRAC_20" in result

            # Check that non-requested windows are not computed
            assert "CTX_MFE_5" not in result
            assert "CTX_MFE_15" not in result
            assert "PAIR_OUT_FRAC_5" not in result

    def test_dna_to_rna_conversion(self, mock_rna):
        """Test that DNA sequences (with T) are converted to RNA (with U)."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        # Track what sequences are passed to fold_compound
        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[0],
                window_size_nt=20,
                flank_left_nt=5,
                flank_right_nt=5,
            )
            # DNA sequence with T
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")

            feature(record)

            # Verify that the sequence passed to ViennaRNA has U instead of T
            assert len(sequences_folded) == 1
            assert "T" not in sequences_folded[0]
            assert "U" in sequences_folded[0]

    def test_pair_out_frac_in_range(self, mock_rna):
        """Test that PAIR_OUT_FRAC values are in [0, 1]."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[0, 10, 20],
                window_size_nt=15,
                flank_left_nt=5,
                flank_right_nt=5,
            )
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")

            result = feature(record)

            # Check that all PAIR_OUT_FRAC values are in [0, 1]
            for key in result:
                if key.startswith("PAIR_OUT_FRAC"):
                    value = result[key]
                    assert 0.0 <= value <= 1.0, f"{key} = {value} is not in [0, 1]"

    def test_boundary_flank_behavior_start_at_zero(self, mock_rna):
        """Test boundary behavior when window starts at 0 (no left flank available)."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[0],
                window_size_nt=10,
                flank_left_nt=5,  # Requested but not available
                flank_right_nt=5,
            )
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="test")  # 20 bp

            result = feature(record)

            # Should compute successfully
            assert "CTX_MFE_0" in result
            assert "PAIR_OUT_FRAC_0" in result

            # Context should be: 0 to 10+5 = 15
            # (no left flank because start is 0)
            assert len(sequences_folded) == 1
            assert len(sequences_folded[0]) == 15

    def test_boundary_flank_behavior_end_extends_beyond(self, mock_rna):
        """Test boundary behavior when right flank extends beyond sequence."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = ("." * len(seq), mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[10],
                window_size_nt=10,
                flank_left_nt=5,
                flank_right_nt=10,  # Will extend beyond sequence
            )
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="test")  # 20 bp

            result = feature(record)

            # Should compute successfully
            assert "CTX_MFE_10" in result
            assert "PAIR_OUT_FRAC_10" in result

            # Context should be: 10-5=5 to min(20, 20+10)=20
            # So 5 to 20 = 15 nucleotides
            assert len(sequences_folded) == 1
            assert len(sequences_folded[0]) == 15

    def test_window_beyond_sequence_skipped(self, mock_rna):
        """Test that windows extending beyond sequence are skipped."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = ContextWindowFoldFeature(
                starts_nt=[0, 10, 25],
                window_size_nt=20,
                flank_left_nt=5,
                flank_right_nt=5,
            )
            # Sequence is only 30 bp, so window at 25 (ends at 45) extends beyond
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTAC"), id="test")

            result = feature(record)

            # Only first two windows should be computed
            assert "CTX_MFE_0" in result
            assert "CTX_MFE_10" in result
            assert "CTX_MFE_25" not in result
            assert "PAIR_OUT_FRAC_0" in result
            assert "PAIR_OUT_FRAC_10" in result
            assert "PAIR_OUT_FRAC_25" not in result

    def test_parse_dotbracket_simple(self):
        """Test dot-bracket parsing with simple structures."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        feature = ContextWindowFoldFeature(starts_nt=[0], window_size_nt=10)

        # Simple hairpin: (((...)))
        pairs = feature._parse_dotbracket("(((...)))")
        assert pairs == {0: 8, 1: 7, 2: 6, 6: 2, 7: 1, 8: 0}

        # No pairs
        pairs = feature._parse_dotbracket(".........")
        assert pairs == {}

        # Mixed
        pairs = feature._parse_dotbracket("..((..))")
        assert pairs == {2: 7, 3: 6, 6: 3, 7: 2}

    def test_compute_pair_out_frac_all_unpaired(self):
        """Test PAIR_OUT_FRAC computation when all nucleotides are unpaired."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        feature = ContextWindowFoldFeature(starts_nt=[0], window_size_nt=10)

        # All unpaired
        frac = feature._compute_pair_out_frac("." * 20, window_start=5, window_end=15)
        assert frac == 0.0

    def test_compute_pair_out_frac_all_paired_within(self):
        """Test PAIR_OUT_FRAC when all window nts are paired within window."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        feature = ContextWindowFoldFeature(starts_nt=[0], window_size_nt=10)

        # Structure: ....((((()))))....
        # Window: positions 4-14 (10 nts)
        # All pairs are within window
        structure = "...." + "((((()))))" + "...."
        frac = feature._compute_pair_out_frac(structure, window_start=4, window_end=14)
        assert frac == 0.0

    def test_compute_pair_out_frac_some_paired_out(self):
        """Test PAIR_OUT_FRAC when some window nts are paired outside window."""
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        feature = ContextWindowFoldFeature(starts_nt=[0], window_size_nt=10)

        # Structure: ((((.......))))
        # Window: positions 4-14 (10 nts in middle)
        # First 4 window nts (positions 4-7) are paired to positions 10-13 (within window)
        # But let's construct: ((((....))))
        # Window 2-12: positions 2,3 paired to 10,11 (within); rest unpaired
        structure = "((((....))))"
        # Window positions 2-12 (10 nts: 2,3,4,5,6,7,8,9,10,11)
        # Pairs: 0-11, 1-10, 2-9, 3-8
        # So window nts: 2->9(within), 3->8(within), 4-7 unpaired,
        # 8->3(within), 9->2(within), 10->1(out), 11->0(out)
        frac = feature._compute_pair_out_frac(structure, window_start=2, window_end=12)
        # Window positions 2,3,4,5,6,7,8,9,10,11
        # 10 nts total
        # 10 pairs with 1 (outside), 11 pairs with 0 (outside)
        # So 2 out of 10 = 0.2
        assert abs(frac - 0.2) < 1e-6

    def test_compare_ctx_mfe_to_direct_viennarna_call(self):
        """Test that context MFE matches a direct ViennaRNA call on the same context slice.

        This test requires ViennaRNA to be installed.
        """
        pytest.importorskip("RNA", reason="ViennaRNA not installed")

        # Import the actual RNA module
        import RNA  # type: ignore[import]

        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        # Test sequence
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        record = SeqRecord(Seq(sequence), id="test")

        # Define windows to test
        starts_nt = [5, 15]
        window_size_nt = 10
        flank_left_nt = 3
        flank_right_nt = 4

        # Compute using ContextWindowFoldFeature
        feature = ContextWindowFoldFeature(
            starts_nt=starts_nt,
            window_size_nt=window_size_nt,
            flank_left_nt=flank_left_nt,
            flank_right_nt=flank_right_nt,
        )
        feature_result = feature(record)

        # Compare each window to direct ViennaRNA call
        for start in starts_nt:
            window_end = start + window_size_nt
            if window_end > len(sequence):
                continue

            # Compute context slice boundaries
            ctx_start = max(0, start - flank_left_nt)
            ctx_end = min(len(sequence), window_end + flank_right_nt)

            # Extract context slice and convert to RNA
            ctx_seq = sequence[ctx_start:ctx_end]
            rna_seq = ctx_seq.replace("T", "U")

            # Direct ViennaRNA call
            fc = RNA.fold_compound(rna_seq)
            structure, mfe = fc.mfe()

            # Compare CTX_MFE
            feature_mfe = feature_result[f"CTX_MFE_{start}"]
            assert abs(feature_mfe - mfe) < 1e-6, (
                f"CTX_MFE mismatch for window at position {start}. "
                f"Feature returned {feature_mfe}, direct call returned {mfe}"
            )

    def test_export_from_init(self):
        """Test that ContextWindowFoldFeature can be imported from __init__."""
        # This test ensures lazy import mechanism works
        from biotooler.families.viennarna import ContextWindowFoldFeature

        # Should be able to instantiate without loading RNA
        feature = ContextWindowFoldFeature(
            starts_nt=[0, 10],
            window_size_nt=20,
            flank_left_nt=5,
            flank_right_nt=5,
        )
        assert feature.starts_nt == [0, 10]
        assert feature.window_size_nt == 20
        assert feature.flank_left_nt == 5
        assert feature.flank_right_nt == 5

    def test_get_features_returns_context_window_fold(self):
        """Test that get_features() includes ContextWindowFoldFeature."""
        pytest.importorskip("RNA", reason="ViennaRNA not installed")

        from biotooler.families.viennarna import get_features

        features = get_features()
        assert len(features) > 0

        # Check that ContextWindowFoldFeature is in the list
        from biotooler.families.viennarna.context_window_fold import (
            ContextWindowFoldFeature,
        )

        assert ContextWindowFoldFeature in features
