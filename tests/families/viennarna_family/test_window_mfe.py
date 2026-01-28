"""Tests for ViennaRNA WindowMFEFeature."""

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

    # Mock fold_compound
    def mock_fold_compound(seq):
        fc_mock = MagicMock()
        # Return a simple MFE calculation based on sequence length
        # Use a deterministic formula for testing
        mfe_value = -len(seq) * 0.5  # Simple formula for testing
        fc_mock.mfe.return_value = (f"{'.' * len(seq)}", mfe_value)
        return fc_mock

    rna_mock.fold_compound = mock_fold_compound
    return rna_mock


class TestWindowMFEFeature:
    """Tests for WindowMFEFeature."""

    def test_import_does_not_load_rna(self):
        """Test that importing window_mfe module doesn't import RNA."""
        # Clear any previously imported modules
        modules_to_clear = [
            module
            for module in list(sys.modules.keys())
            if module.startswith("biotooler.families.viennarna.window_mfe")
            or module == "RNA"
            or module.startswith("RNA.")
        ]
        for module in modules_to_clear:
            del sys.modules[module]

        # Import the window_mfe module

        # Verify RNA is NOT loaded
        assert "RNA" not in sys.modules, (
            "RNA was imported when importing window_mfe module, "
            "but it should only be loaded when features are instantiated/called"
        )

    def test_init_validates_parameters(self):
        """Test that __init__ validates parameters."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Test empty window_starts
        with pytest.raises(ValueError, match="window_starts must not be empty"):
            WindowMFEFeature(window_starts=[], window_size=20)

        # Test non-positive window_size
        with pytest.raises(ValueError, match="window_size must be positive"):
            WindowMFEFeature(window_starts=[0, 10], window_size=0)

        with pytest.raises(ValueError, match="window_size must be positive"):
            WindowMFEFeature(window_starts=[0, 10], window_size=-5)

    def test_compute_mfe_for_requested_windows(self, mock_rna):
        """Test that MFE is computed only for requested window positions."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Patch require_viennarna to return our mock
        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = WindowMFEFeature(window_starts=[0, 10, 20], window_size=15)
            # Use 40 bp sequence so all windows fit
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"), id="test")  # 40 bp

            result = feature(record)

            # Check that only requested windows are computed
            assert "MFE_0" in result
            assert "MFE_10" in result
            assert "MFE_20" in result

            # Check that non-requested windows are not computed
            assert "MFE_5" not in result
            assert "MFE_15" not in result
            assert "MFE_25" not in result

            # Verify MFE values (based on our mock formula: -len * 0.5)
            assert result["MFE_0"] == -7.5  # 15 * 0.5
            assert result["MFE_10"] == -7.5
            assert result["MFE_20"] == -7.5

    def test_dna_to_rna_conversion(self, mock_rna):
        """Test that DNA sequences (with T) are converted to RNA (with U)."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Track what sequences are passed to fold_compound
        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = (f"{'.' * len(seq)}", mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = WindowMFEFeature(window_starts=[0], window_size=20)
            # DNA sequence with T
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGT"), id="test")

            feature(record)

            # Verify that the sequence passed to ViennaRNA has U instead of T
            assert len(sequences_folded) == 1
            assert "T" not in sequences_folded[0]
            assert "U" in sequences_folded[0]
            assert sequences_folded[0] == "ACGUACGUACGUACGUACGU"

    def test_rna_sequence_unchanged(self, mock_rna):
        """Test that RNA sequences (with U) are passed unchanged."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        sequences_folded = []

        def track_fold_compound(seq):
            sequences_folded.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = (f"{'.' * len(seq)}", mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = WindowMFEFeature(window_starts=[0], window_size=20)
            # RNA sequence with U
            record = SeqRecord(Seq("ACGUACGUACGUACGUACGU"), id="test")

            feature(record)

            # Verify that the sequence is unchanged
            assert len(sequences_folded) == 1
            assert sequences_folded[0] == "ACGUACGUACGUACGUACGU"

    def test_window_beyond_sequence_skipped(self, mock_rna):
        """Test that windows extending beyond sequence are skipped."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            feature = WindowMFEFeature(window_starts=[0, 10, 25], window_size=20)
            # Sequence is only 30 bp, so window at 25 (ends at 45) extends beyond
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTAC"), id="test")

            result = feature(record)

            # Only first two windows should be computed
            assert "MFE_0" in result
            assert "MFE_10" in result
            assert "MFE_25" not in result  # Extends beyond sequence

    def test_only_requested_starts_computed_via_monkeypatch(self, mock_rna):
        """Test that only requested starts are computed by counting fold_compound calls."""
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        fold_compound_calls = []

        def track_fold_compound(seq):
            fold_compound_calls.append(seq)
            fc_mock = MagicMock()
            mfe_value = -len(seq) * 0.5
            fc_mock.mfe.return_value = (f"{'.' * len(seq)}", mfe_value)
            return fc_mock

        mock_rna.fold_compound = track_fold_compound

        with patch(
            "biotooler.families.viennarna.integration.require_viennarna",
            return_value=mock_rna,
        ):
            # Request only 3 specific windows out of many possible
            feature = WindowMFEFeature(window_starts=[5, 15, 25], window_size=10)
            # 50 bp sequence could have many windows, but we request only 3
            record = SeqRecord(Seq("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"), id="test")

            result = feature(record)

            # Verify exactly 3 fold_compound calls (one per requested window)
            assert len(fold_compound_calls) == 3, (
                f"Expected 3 fold_compound calls but got {len(fold_compound_calls)}. "
                "This means the feature computed windows that weren't requested."
            )

            # Verify the correct windows were computed
            assert "MFE_5" in result
            assert "MFE_15" in result
            assert "MFE_25" in result

    def test_compare_mfe_to_direct_viennarna_call(self):
        """Test that each window's MFE matches a direct ViennaRNA call on the same substring.

        This test requires ViennaRNA to be installed.
        """
        pytest.importorskip("RNA", reason="ViennaRNA not installed")

        # Import the actual RNA module
        import RNA  # type: ignore[import]

        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        # Test sequence
        sequence = "ACGTACGTACGTACGTACGTACGTACGTACGT"
        record = SeqRecord(Seq(sequence), id="test")

        # Define windows to test
        window_starts = [0, 10, 15]
        window_size = 15

        # Compute using WindowMFEFeature
        feature = WindowMFEFeature(window_starts=window_starts, window_size=window_size)
        feature_result = feature(record)

        # Compare each window to direct ViennaRNA call
        for start in window_starts:
            end = start + window_size
            if end > len(sequence):
                continue

            # Extract substring and convert to RNA
            substring = sequence[start:end]
            rna_substring = substring.replace("T", "U")

            # Direct ViennaRNA call
            fc = RNA.fold_compound(rna_substring)
            structure, mfe = fc.mfe()

            # Compare
            feature_mfe = feature_result[f"MFE_{start}"]
            assert abs(feature_mfe - mfe) < 1e-6, (
                f"MFE mismatch for window at position {start}. "
                f"Feature returned {feature_mfe}, direct call returned {mfe}"
            )

    def test_export_from_init(self):
        """Test that WindowMFEFeature can be imported from __init__."""
        # This test ensures lazy import mechanism works
        from biotooler.families.viennarna import WindowMFEFeature

        # Should be able to instantiate without loading RNA
        feature = WindowMFEFeature(window_starts=[0, 10], window_size=20)
        assert feature.window_starts == [0, 10]
        assert feature.window_size == 20

    def test_get_features_returns_window_mfe(self):
        """Test that get_features() includes WindowMFEFeature."""
        pytest.importorskip("RNA", reason="ViennaRNA not installed")

        from biotooler.families.viennarna import get_features

        features = get_features()
        assert len(features) > 0

        # Check that WindowMFEFeature is in the list
        from biotooler.families.viennarna.window_mfe import WindowMFEFeature

        assert WindowMFEFeature in features
