"""Feature computation sets for sequence analysis."""

from collections.abc import Callable
from typing import Any

import pandas as pd
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import OrfSpan, get_orf, select_orf_by_index
from biotooler.core.types import FeatureOutput
from biotooler.core.windowing import iter_orf_codon_windows


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
        for window in windows:
            window_start = window.annotations["window_start"]

            # Compute all features for this window
            for _feat_name, feat_fn in self._features.items():
                features = feat_fn(window)
                for key, value in features.items():
                    # Create column name: {name}.{feature_key}_{window_start}
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
