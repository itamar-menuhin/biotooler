"""Tests for incremental ORF window computation (wide format)."""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.features import FeatureSet


class ToyIncrementalFeature:
    """A simple incremental feature that counts nucleotides.

    This feature maintains a running count of each nucleotide type (A, C, G, T)
    as windows slide across the sequence. It demonstrates the incremental interface
    where state is updated by removing outgoing bases and adding incoming bases.
    """

    def __call__(self, record: SeqRecord) -> dict[str, int]:
        """Make the feature callable for FeatureSet compatibility.

        This method provides the same output as the incremental methods but
        computes it from scratch for each window.
        """
        seq = str(record.seq).upper()
        return {
            "a_count": seq.count("A"),
            "c_count": seq.count("C"),
            "g_count": seq.count("G"),
            "t_count": seq.count("T"),
        }

    def init_state(
        self,
        record: SeqRecord,
        *,
        orf: tuple[int, int],
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> dict:
        """Initialize state for the first window."""
        # Extract the window sequence
        orf_start, orf_end = orf
        abs_start = orf_start + window_start
        abs_end = orf_start + window_end
        seq = str(record.seq)[abs_start:abs_end].upper()

        # Count nucleotides in first window
        state = {
            "record": record,
            "orf": orf,
            "a_count": seq.count("A"),
            "c_count": seq.count("C"),
            "g_count": seq.count("G"),
            "t_count": seq.count("T"),
        }
        return state

    def step_state(
        self,
        state: dict,
        *,
        out_start: int,
        out_end: int,
        in_start: int,
        in_end: int,
        **kwargs,
    ) -> None:
        """Update state for the next window by removing outgoing and adding incoming bases."""
        record = state["record"]
        orf = state["orf"]
        orf_start = orf[0]

        # Remove outgoing bases
        if out_end > out_start:
            abs_out_start = orf_start + out_start
            abs_out_end = orf_start + out_end
            out_seq = str(record.seq)[abs_out_start:abs_out_end].upper()
            state["a_count"] -= out_seq.count("A")
            state["c_count"] -= out_seq.count("C")
            state["g_count"] -= out_seq.count("G")
            state["t_count"] -= out_seq.count("T")

        # Add incoming bases
        if in_end > in_start:
            abs_in_start = orf_start + in_start
            abs_in_end = orf_start + in_end
            in_seq = str(record.seq)[abs_in_start:abs_in_end].upper()
            state["a_count"] += in_seq.count("A")
            state["c_count"] += in_seq.count("C")
            state["g_count"] += in_seq.count("G")
            state["t_count"] += in_seq.count("T")

    def emit(self, state: dict) -> dict[str, int]:
        """Emit feature values from current state."""
        return {
            "a_count": state["a_count"],
            "c_count": state["c_count"],
            "g_count": state["g_count"],
            "t_count": state["t_count"],
        }


class ToyNonIncrementalFeature:
    """Same feature as ToyIncrementalFeature but without incremental methods.

    This class only has the __call__ method, forcing the fallback path.
    """

    def __call__(self, record: SeqRecord) -> dict[str, int]:
        """Non-incremental computation."""
        seq = str(record.seq).upper()
        return {
            "a_count": seq.count("A"),
            "c_count": seq.count("C"),
            "g_count": seq.count("G"),
            "t_count": seq.count("T"),
        }


def test_incremental_path_produces_correct_output():
    """Test that incremental computation produces correct feature values."""
    feature = ToyIncrementalFeature()
    fs = FeatureSet(feature, name="test")

    # Sequence: ATGAAACCCGGGTTT
    # Windows (9nt, step=3):
    # 0-9:   ATGAAACCC -> A:4, C:3, G:1, T:1
    # 3-12:  AAACCCGGG -> A:3, C:3, G:3, T:0
    # 6-15:  CCCGGGTTT -> A:0, C:3, G:3, T:3
    record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")
    result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

    # Check window 0 (ATGAAACCC)
    assert result["test.a_count_0"].iloc[0] == 4
    assert result["test.c_count_0"].iloc[0] == 3
    assert result["test.g_count_0"].iloc[0] == 1
    assert result["test.t_count_0"].iloc[0] == 1

    # Check window 3 (AAACCCGGG)
    assert result["test.a_count_3"].iloc[0] == 3
    assert result["test.c_count_3"].iloc[0] == 3
    assert result["test.g_count_3"].iloc[0] == 3
    assert result["test.t_count_3"].iloc[0] == 0

    # Check window 6 (CCCGGGTTT)
    assert result["test.a_count_6"].iloc[0] == 0
    assert result["test.c_count_6"].iloc[0] == 3
    assert result["test.g_count_6"].iloc[0] == 3
    assert result["test.t_count_6"].iloc[0] == 3


def test_fallback_path_produces_correct_output():
    """Test that fallback (non-incremental) computation produces correct output."""
    feature = ToyNonIncrementalFeature()
    fs = FeatureSet(feature, name="test")

    record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")
    result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

    # Check window 0 (ATGAAACCC)
    assert result["test.a_count_0"].iloc[0] == 4
    assert result["test.c_count_0"].iloc[0] == 3
    assert result["test.g_count_0"].iloc[0] == 1
    assert result["test.t_count_0"].iloc[0] == 1

    # Check window 3 (AAACCCGGG)
    assert result["test.a_count_3"].iloc[0] == 3
    assert result["test.c_count_3"].iloc[0] == 3
    assert result["test.g_count_3"].iloc[0] == 3
    assert result["test.t_count_3"].iloc[0] == 0

    # Check window 6 (CCCGGGTTT)
    assert result["test.a_count_6"].iloc[0] == 0
    assert result["test.c_count_6"].iloc[0] == 3
    assert result["test.g_count_6"].iloc[0] == 3
    assert result["test.t_count_6"].iloc[0] == 3


def test_incremental_and_fallback_produce_identical_output():
    """Test that incremental and fallback paths produce exactly the same output."""
    # Create incremental version
    incremental_feature = ToyIncrementalFeature()
    fs_incremental = FeatureSet(incremental_feature, name="test")

    # Create fallback version (no incremental methods)
    fallback_feature = ToyNonIncrementalFeature()
    fs_fallback = FeatureSet(fallback_feature, name="test")

    # Use same test sequence
    record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")

    # Compute with both approaches
    result_incremental = fs_incremental.compute_orf_windows(
        record, orf=(0, 15), window_nt=9, step_nt=3
    )
    result_fallback = fs_fallback.compute_orf_windows(
        record, orf=(0, 15), window_nt=9, step_nt=3
    )

    # Assert DataFrames are exactly equal
    # Sort columns to ensure comparison works
    result_incremental = result_incremental.sort_index(axis=1)
    result_fallback = result_fallback.sort_index(axis=1)

    # Check shapes match
    assert result_incremental.shape == result_fallback.shape

    # Check columns match
    assert list(result_incremental.columns) == list(result_fallback.columns)

    # Check values match for all columns
    for col in result_incremental.columns:
        assert result_incremental[col].iloc[0] == result_fallback[col].iloc[0], (
            f"Mismatch in column {col}: "
            f"{result_incremental[col].iloc[0]} != {result_fallback[col].iloc[0]}"
        )


def test_incremental_with_multiple_windows():
    """Test incremental computation with many windows."""
    feature = ToyIncrementalFeature()
    fs = FeatureSet(feature, name="test")

    # Longer sequence with more windows
    # ATG + 18 codons = 57 nt total
    seq = "ATG" + "AAA" * 6 + "CCC" * 6 + "GGG" * 6
    record = SeqRecord(Seq(seq), id="test")

    result = fs.compute_orf_windows(record, orf=(0, 57), window_nt=12, step_nt=3)

    # Should have many windows
    # (57 - 12) // 3 + 1 = 16 windows
    feature_cols = [c for c in result.columns if c.startswith("test.a_count_")]
    assert len(feature_cols) == 16

    # Verify column ordering is deterministic
    window_starts = [int(c.rsplit("_", 1)[1]) for c in feature_cols]
    assert window_starts == sorted(window_starts)


def test_incremental_preserves_metadata():
    """Test that incremental path preserves metadata columns."""
    feature = ToyIncrementalFeature()
    fs = FeatureSet(feature, name="test")

    record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="my_seq")
    result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

    # Check metadata columns
    assert "record_id" in result.columns
    assert "orf_start" in result.columns
    assert "orf_end" in result.columns
    assert result["record_id"].iloc[0] == "my_seq"
    assert result["orf_start"].iloc[0] == 0
    assert result["orf_end"].iloc[0] == 15


def test_incremental_preserves_column_ordering():
    """Test that incremental path preserves deterministic column ordering."""
    feature = ToyIncrementalFeature()
    fs = FeatureSet(feature, name="test")

    record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")
    result = fs.compute_orf_windows(record, orf=(0, 15), window_nt=9, step_nt=3)

    # Check that metadata columns come first
    cols = list(result.columns)
    assert cols[0] == "record_id"
    assert cols[1] == "orf_start"
    assert cols[2] == "orf_end"

    # Check that feature columns are sorted by (feature_key, window_start)
    feature_cols = [c for c in cols if c.startswith("test.")]
    # Should be: a_count_0, a_count_3, a_count_6, c_count_0, c_count_3, c_count_6, ...
    assert feature_cols[0] == "test.a_count_0"
    assert feature_cols[1] == "test.a_count_3"
    assert feature_cols[2] == "test.a_count_6"
    assert feature_cols[3] == "test.c_count_0"
