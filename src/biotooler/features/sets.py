"""Feature computation sets for sequence analysis."""

from collections.abc import Callable
from typing import Any

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import OrfSpan, get_orf, select_orf_by_index
from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import FeatureOutput
from biotooler.core.windowing import iter_orf_codon_windows, iter_windows


class FeatureSet:
    """A collection of feature computations that can be applied to sequences.

    This class wraps one or more feature computation functions and provides methods
    to compute features over sequences or ORF windows.

    Args:
        features: Dictionary mapping feature names to computation functions,
                 or a single computation function
        name: Optional name for this feature set (used in column naming)

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> # Define feature functions
        >>> def compute_gc(record):
        ...     seq = str(record.seq)
        ...     gc_count = seq.count('G') + seq.count('C')
        ...     return {"gc_content": gc_count / len(seq) if seq else 0.0}
        >>> # Create feature set
        >>> fs = FeatureSet(compute_gc, name="gc")
    """

    def __init__(
        self,
        features: (
            Callable[[SeqRecord], FeatureOutput]
            | dict[str, Callable[[SeqRecord], FeatureOutput]]
        ),
        *,
        name: str = "features",
    ):
        """Initialize the FeatureSet.

        Args:
            features: Single function or dict of functions that compute features from SeqRecords
            name: Name for this feature set (used in output)
        """
        if callable(features):
            # Single function - wrap it
            self._features = {"compute": features}
        else:
            self._features = features
        self.name = name

    def compute_orf_windows(
        self,
        record: SeqRecord,
        *,
        orf: OrfSpan | None = None,
        orf_index: int | None = None,
        window_nt: int,
        step_nt: int,
        drop_partial: bool = True,
    ) -> pd.DataFrame:
        """Compute features for ORF windows and return in WIDE format.

        This method produces wide-format output: one row per record with columns suffixed
        by window start index (e.g., CAI_0, CAI_3, CAI_6 for step_nt=3).

        ORF Resolution Rules:
        - If `orf` is provided: use it directly
        - Else if `orf_index` is provided: find candidates and select that index
        - Else: try to get attached ORF from record; if missing, raise clear error

        Output Format:
        - Single row per record
        - Metadata columns first: record_id, orf_start, orf_end
        - Feature columns with suffixes: {name}.{feature_key}_{window_start}
        - Column ordering is deterministic: metadata first, then features sorted by
          (feature_key, window_start numeric)

        Args:
            record: DNA/RNA SeqRecord to analyze
            orf: Optional explicit ORF coordinates (start, end)
            orf_index: Optional index to select from ORF candidates
            window_nt: Size of each window in nucleotides
            step_nt: Step size between windows (must be multiple of 3)
            drop_partial: If True, drops partial windows at end

        Returns:
            Single-row DataFrame with wide-format features

        Raises:
            ValueError: If step_nt is not multiple of 3, or if protein sequence provided,
                       or if neither orf, orf_index, nor attached ORF is available
            KeyError: If trying to use attached ORF but none exists

        Examples:
            >>> from Bio.Seq import Seq
            >>> from Bio.SeqRecord import SeqRecord
            >>> def count_a(rec):
            ...     return {"a_count": str(rec.seq).count('A')}
            >>> fs = FeatureSet(count_a, name="count")
            >>> record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")
            >>> # Use explicit ORF
            >>> result = fs.compute_orf_windows(
            ...     record, orf=(0, 15), window_nt=9, step_nt=3
            ... )
            >>> result.shape
            (1, 6)
            >>> 'count.a_count_0' in result.columns
            True
            >>> 'count.a_count_3' in result.columns
            True
            >>> 'count.a_count_6' in result.columns
            True
        """
        # Check for protein sequences early with clear error message
        if "molecule_type" in record.annotations:
            mol_type = record.annotations["molecule_type"]
            if isinstance(mol_type, str) and mol_type.upper() == "PROTEIN":
                raise ValueError(
                    "ORF window computation is only supported for DNA/RNA sequences, "
                    "not protein sequences"
                )

        # Resolve ORF coordinates
        resolved_orf: OrfSpan
        if orf is not None:
            # Use explicit ORF
            resolved_orf = orf
        elif orf_index is not None:
            # Find candidates and select by index
            candidates = find_orf_candidates(record)
            resolved_orf = select_orf_by_index(candidates, orf_index)
        else:
            # Try to get attached ORF - distinguish missing ORF from protein error
            try:
                resolved_orf = get_orf(record)
            except KeyError as e:
                raise ValueError(
                    f"No ORF information provided for record {record.id!r}. "
                    "Please provide 'orf' parameter, 'orf_index' parameter, "
                    "or attach an ORF using attach_orf()."
                ) from e

        # Extract windows using iter_orf_codon_windows
        windows = list(
            iter_orf_codon_windows(
                record,
                orf=resolved_orf,
                window_nt=window_nt,
                step_nt=step_nt,
                drop_partial=drop_partial,
            )
        )

        # Build feature data dictionary for single row
        feature_data: dict[str, Any] = {
            "record_id": record.id,
            "orf_start": resolved_orf[0],
            "orf_end": resolved_orf[1],
        }

        if not windows:
            # No windows generated - return row with metadata only
            return pd.DataFrame([feature_data])

        # Collect all feature keys from first window to ensure we have all columns
        all_feature_keys: set[str] = set()

        # Compute features for each window and add with window_start suffix
        for _feat_name, feat_fn in self._features.items():
            # Check if feature supports incremental computation (duck-typing)
            has_incremental = (
                hasattr(feat_fn, "init_state")
                and hasattr(feat_fn, "step_state")
                and hasattr(feat_fn, "emit")
            )

            if has_incremental:
                # Use incremental path
                state = None
                prev_window_start = 0
                prev_window_end = 0

                for window in windows:
                    window_start = window.annotations["window_start"]
                    window_end = window.annotations["window_end"]

                    if state is None:
                        # First window: initialize state
                        state = feat_fn.init_state(  # type: ignore[union-attr]
                            record,
                            orf=resolved_orf,
                            window_start=window_start,
                            window_end=window_end,
                        )
                    else:
                        # Subsequent windows: update state incrementally
                        # out: bases leaving the window (from prev_start to curr_start)
                        # in: bases entering the window (from prev_end to curr_end)
                        out_start = prev_window_start
                        out_end = window_start
                        in_start = prev_window_end
                        in_end = window_end

                        feat_fn.step_state(  # type: ignore[union-attr]
                            state,
                            out_start=out_start,
                            out_end=out_end,
                            in_start=in_start,
                            in_end=in_end,
                        )

                    prev_window_start = window_start
                    prev_window_end = window_end

                    # Emit features for this window
                    features = feat_fn.emit(state)  # type: ignore[union-attr]
                    for key, value in features.items():
                        col_name = f"{self.name}.{key}_{window_start}"
                        feature_data[col_name] = value
                        all_feature_keys.add(key)
            else:
                # Use fallback path for non-incremental features
                for window in windows:
                    window_start = window.annotations["window_start"]
                    features = feat_fn(window)
                    for key, value in features.items():
                        col_name = f"{self.name}.{key}_{window_start}"
                        feature_data[col_name] = value
                        all_feature_keys.add(key)

        # Create single-row DataFrame
        df = pd.DataFrame([feature_data])

        # Ensure deterministic column ordering:
        # 1. Metadata columns first
        metadata_cols = ["record_id", "orf_start", "orf_end"]

        # 2. Feature columns sorted by (feature_key, window_start numeric)
        feature_cols = [c for c in df.columns if c not in metadata_cols]

        # Sort with error handling for malformed column names
        def sort_key(col: str) -> tuple[str, int]:
            parts = col.rsplit("_", 1)
            if len(parts) != 2:
                # No underscore found - sort by column name only
                return (col, 0)
            try:
                return (parts[0], int(parts[1]))
            except ValueError:
                # Can't parse as int - sort by full name
                return (col, 0)

        feature_cols.sort(key=sort_key)

        # Reorder columns
        ordered_cols = metadata_cols + feature_cols
        result_df = df[ordered_cols]
        assert isinstance(result_df, pd.DataFrame)

        return result_df

    def compute_windows(
        self,
        record: SeqRecord,
        *,
        window_size: int,
        step: int,
        region: tuple[int, int] | None = None,
        drop_partial: bool = True,
    ) -> pd.DataFrame:
        """Compute features for generic sliding windows and return in WIDE format.

        This method produces wide-format output: one row per record with columns suffixed
        by window start index (e.g., gc_content_0, gc_content_5, gc_content_10 for step=5).

        Unlike compute_orf_windows, this method:
        - Works with any sequence type (DNA/RNA/protein)
        - Does not require codon alignment (step can be any positive integer)
        - Uses region parameter instead of ORF coordinates
        - Window starts are relative to region start (first window is always _0)

        Output Format:
        - Single row per record
        - Metadata columns first: record_id, region_start, region_end
        - Feature columns with suffixes: {name}.{feature_key}_{window_start}
        - Column ordering is deterministic: metadata first, then features sorted by
          (feature_key, window_start numeric)

        Args:
            record: SeqRecord to analyze (DNA/RNA/protein)
            window_size: Size of each window in residues/bases
            step: Step size between windows
            region: Optional tuple (start, end) to restrict windowing to a subsequence.
                   If None, uses entire sequence (0, len(record.seq))
            drop_partial: If True, drops partial windows at end

        Returns:
            Single-row DataFrame with wide-format features

        Raises:
            ValueError: If window_size or step is not positive

        Examples:
            >>> from Bio.Seq import Seq
            >>> from Bio.SeqRecord import SeqRecord
            >>> def compute_length(rec):
            ...     return {"len": len(rec.seq)}
            >>> fs = FeatureSet(compute_length, name="basic")
            >>> record = SeqRecord(Seq("ACGTACGTACGT"), id="test")
            >>> # Use sliding windows with step=3
            >>> result = fs.compute_windows(
            ...     record, window_size=6, step=3
            ... )
            >>> result.shape
            (1, 6)
            >>> 'basic.len_0' in result.columns
            True
            >>> 'basic.len_3' in result.columns
            True
            >>> 'basic.len_6' in result.columns
            True
        """
        # Resolve region coordinates
        seq_str = get_seq_str(record)
        if region is None:
            region_start = 0
            region_end = len(seq_str)
            region_record = record
        else:
            region_start, region_end = region
            # Create a subrecord for the region
            region_seq_str = seq_str[region_start:region_end]
            region_record = SeqRecord(
                Seq(region_seq_str),
                id=record.id,
                description=record.description,
            )
            # Copy annotations but exclude cache keys to avoid using cached full sequence
            for key, value in record.annotations.items():
                if not key.startswith("_biotooler_"):
                    region_record.annotations[key] = value

        # Extract windows using iter_windows
        windows = list(
            iter_windows(
                region_record,
                window_size=window_size,
                step=step,
                drop_partial=drop_partial,
            )
        )

        # Build feature data dictionary for single row
        feature_data: dict[str, Any] = {
            "record_id": record.id,
            "region_start": region_start,
            "region_end": region_end,
        }

        if not windows:
            # No windows generated - return row with metadata only
            return pd.DataFrame([feature_data])

        # Collect all feature keys from first window to ensure we have all columns
        all_feature_keys: set[str] = set()

        # Compute features for each window and add with window_start suffix
        for _feat_name, feat_fn in self._features.items():
            # Check if feature supports incremental computation (duck-typing)
            has_incremental = (
                hasattr(feat_fn, "init_state")
                and hasattr(feat_fn, "step_state")
                and hasattr(feat_fn, "emit")
            )

            if has_incremental:
                # Use incremental path
                state = None
                prev_window_start = 0
                prev_window_end = 0

                for window in windows:
                    # For generic windows, use "start" annotation (relative to region_record)
                    window_start = window.annotations["start"]
                    window_end = window.annotations["end"]

                    if state is None:
                        # First window: initialize state
                        # Pass region parameter for consistency with ORF windowing
                        state = feat_fn.init_state(  # type: ignore[union-attr]
                            region_record,
                            window_start=window_start,
                            window_end=window_end,
                            region=(region_start, region_end),
                        )
                    else:
                        # Subsequent windows: update state incrementally
                        # out: bases leaving the window (from prev_start to curr_start)
                        # in: bases entering the window (from prev_end to curr_end)
                        out_start = prev_window_start
                        out_end = window_start
                        in_start = prev_window_end
                        in_end = window_end

                        feat_fn.step_state(  # type: ignore[union-attr]
                            state,
                            out_start=out_start,
                            out_end=out_end,
                            in_start=in_start,
                            in_end=in_end,
                        )

                    prev_window_start = window_start
                    prev_window_end = window_end

                    # Emit features for this window
                    features = feat_fn.emit(state)  # type: ignore[union-attr]
                    for key, value in features.items():
                        col_name = f"{self.name}.{key}_{window_start}"
                        feature_data[col_name] = value
                        all_feature_keys.add(key)
            else:
                # Use fallback path for non-incremental features
                for window in windows:
                    window_start = window.annotations["start"]
                    features = feat_fn(window)
                    for key, value in features.items():
                        col_name = f"{self.name}.{key}_{window_start}"
                        feature_data[col_name] = value
                        all_feature_keys.add(key)

        # Create single-row DataFrame
        df = pd.DataFrame([feature_data])

        # Ensure deterministic column ordering:
        # 1. Metadata columns first
        metadata_cols = ["record_id", "region_start", "region_end"]

        # 2. Feature columns sorted by (feature_key, window_start numeric)
        feature_cols = [c for c in df.columns if c not in metadata_cols]

        # Sort with error handling for malformed column names
        def sort_key(col: str) -> tuple[str, int]:
            parts = col.rsplit("_", 1)
            if len(parts) != 2:
                # No underscore found - sort by column name only
                return (col, 0)
            try:
                return (parts[0], int(parts[1]))
            except ValueError:
                # Can't parse as int - sort by full name
                return (col, 0)

        feature_cols.sort(key=sort_key)

        # Reorder columns
        ordered_cols = metadata_cols + feature_cols
        result_df = df[ordered_cols]
        assert isinstance(result_df, pd.DataFrame)

        return result_df
