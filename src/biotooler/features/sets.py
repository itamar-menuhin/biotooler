"""Feature computation sets for sequence analysis."""

from collections.abc import Callable
from typing import Any

import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.core.orf_store import OrfSpan, get_orf, select_orf_by_index
from biotooler.core.seq_utils import get_seq_str
from biotooler.core.types import FeatureOutput
from biotooler.core.windowing import (
    compute_window_indices,
    iter_orf_codon_windows,
    iter_windows,
)
from biotooler.features.aggregation import PositionSpace


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
        elif isinstance(features, dict):
            self._features = features
        else:
            # Assume it's a positional feature object (or other feature object)
            self._features = {"compute": features}
        self.name = name

    def _compute_features_over_windows(
        self,
        windows: list[SeqRecord],
        feature_data: dict[str, Any],
        window_start_key: str,
        init_state_kwargs: dict[str, Any],
    ) -> None:
        """Compute features over windows and populate feature_data dict.

        This is a shared helper method used by both compute_orf_windows and compute_windows
        to compute features over a list of windows, supporting both incremental and
        non-incremental feature computation.

        Args:
            windows: List of SeqRecord windows to compute features over
            feature_data: Dictionary to populate with feature values
            window_start_key: Annotation key to use for window start position
            init_state_kwargs: Additional kwargs to pass to init_state for incremental features
        """
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
                    window_start = window.annotations[window_start_key]
                    # Determine window_end key based on context
                    if "window_end" in window.annotations:
                        window_end = window.annotations["window_end"]
                    else:
                        window_end = window.annotations["end"]

                    if state is None:
                        # First window: initialize state
                        state = feat_fn.init_state(  # type: ignore[union-attr]
                            window_start=window_start,
                            window_end=window_end,
                            **init_state_kwargs,
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
                        col_name = f"{self.name.upper()}_{key}_{window_start}"
                        feature_data[col_name] = value
            else:
                # Use fallback path for non-incremental features
                for window in windows:
                    window_start = window.annotations[window_start_key]
                    features = feat_fn(window)
                    for key, value in features.items():
                        col_name = f"{self.name.upper()}_{key}_{window_start}"
                        feature_data[col_name] = value

    def _format_wide_dataframe(
        self, feature_data: dict[str, Any], metadata_cols: list[str]
    ) -> pd.DataFrame:
        """Format feature data into wide DataFrame with deterministic column ordering.

        This is a shared helper method used by both compute_orf_windows and compute_windows
        to create the final wide-format DataFrame with proper column ordering.

        Args:
            feature_data: Dictionary containing metadata and feature values
            metadata_cols: List of metadata column names to place first

        Returns:
            Single-row DataFrame with ordered columns
        """
        # Create single-row DataFrame
        df = pd.DataFrame([feature_data])

        # Ensure deterministic column ordering:
        # 1. Metadata columns first
        # 2. Feature columns sorted by (feature_key, window_start numeric)
        feature_cols = [c for c in df.columns if c not in metadata_cols]

        # Sort with error handling for malformed column names
        def sort_key(col: str) -> tuple[str, int]:
            """Sort key for feature columns.
            
            Handles both numeric suffixes (e.g., FAMILY_KEY_0) and GLOBAL suffix.
            GLOBAL is treated as a large number to sort after all numeric windows.
            """
            parts = col.rsplit("_", 1)
            if len(parts) != 2:
                # No underscore found - sort by column name only
                return (col, 0)
            
            # Check if last part is GLOBAL
            if parts[1] == "GLOBAL":
                # GLOBAL sorts after all numeric windows
                return (parts[0], float('inf'))
            
            # Try to parse as numeric window position
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

        Windowing Semantics:
        - step_nt moves the window start position forward
        - Each window aggregates ALL codons/residues in its [start, end) range
        - step_nt does NOT subsample codons/residues; it only controls window placement

        For positional features (those implementing position_space, vector_keys, and
        compute_vector), this method automatically uses v2 semantics: compute per-position
        values once across the full ORF, then aggregate by window boundaries. This ensures
        correct positional feature computation (e.g., codon bias, chimera).

        For legacy/incremental features, this method uses the original windowing approach:
        slice windows and call feature(window) or use incremental init/step/emit.

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

        # Validate step_nt (must be multiple of 3 for codon alignment)
        if step_nt % 3 != 0:
            raise ValueError(
                f"step_nt must be a multiple of 3 for codon-aligned windows, got {step_nt}"
            )

        # Build feature data dictionary for single row
        feature_data: dict[str, Any] = {
            "record_id": record.id,
            "orf_start": resolved_orf[0],
            "orf_end": resolved_orf[1],
        }

        # Separate positional features from non-positional features
        positional_features = []
        non_positional_features = []

        for feat_name, feat_fn in self._features.items():
            # Check if feature implements PositionalFeature protocol (duck-typing)
            has_positional = (
                hasattr(feat_fn, "position_space")
                and hasattr(feat_fn, "vector_keys")
                and hasattr(feat_fn, "compute_vector")
            )

            # Check if feature implements IncrementalFeature protocol (duck-typing)
            has_incremental = (
                hasattr(feat_fn, "init_state")
                and hasattr(feat_fn, "step_state")
                and hasattr(feat_fn, "emit")
            )

            # Routing logic:
            # 1. If feature has positional interface with non-empty vector_keys, check if complete
            #    - If complete (all outputs in vector_keys), prefer v2 path
            #    - If incomplete (mixed models like CodonBiasFeature with CAI+ENC), use incremental
            # 2. Else if feature has incremental interface, use legacy path
            # 3. Otherwise, use legacy path
            #
            # This ensures V2 windowing is used by default when features fully support it,
            # while allowing incremental path for mixed/hybrid features.

            if has_positional:
                try:
                    vector_keys = feat_fn.vector_keys  # type: ignore[union-attr]
                    if vector_keys:
                        # Check if mixed feature (has incremental + some models
                        # without get_vector). For features like CodonBiasFeature,
                        # check if names/outputs match vector_keys
                        is_mixed = False
                        if has_incremental and hasattr(feat_fn, 'names'):
                            # CodonBiasFeature case: check if all model names are in vector_keys
                            names = getattr(feat_fn, 'names', [])
                            if names and set(names) != set(vector_keys.keys()):
                                is_mixed = True

                        if is_mixed:
                            # Mixed feature: use incremental path to handle all models
                            non_positional_features.append((feat_name, feat_fn))
                        else:
                            # Pure positional: use V2 path
                            positional_features.append((feat_name, feat_fn))
                    elif has_incremental:
                        # Empty vector_keys but has incremental
                        # (e.g., CodonBiasFeature with ENC only)
                        non_positional_features.append((feat_name, feat_fn))
                    else:
                        # Empty vector_keys and no incremental
                        non_positional_features.append((feat_name, feat_fn))
                except Exception:
                    # If accessing vector_keys raises an exception, try incremental
                    if has_incremental:
                        non_positional_features.append((feat_name, feat_fn))
                    else:
                        non_positional_features.append((feat_name, feat_fn))
            elif has_incremental:
                # Only has incremental, no positional
                non_positional_features.append((feat_name, feat_fn))
            else:
                # Neither positional nor incremental - use legacy fallback
                non_positional_features.append((feat_name, feat_fn))

        # Process positional features using v2 semantics
        if positional_features:
            # Extract ORF sequence
            orf_start, orf_end = resolved_orf
            seq_str = get_seq_str(record)
            orf_seq_str = seq_str[orf_start:orf_end]
            orf_len = len(orf_seq_str)

            # Create ORF record for compute_vector
            orf_record = SeqRecord(
                Seq(orf_seq_str),
                id=record.id,
                description=record.description,
            )
            # Copy annotations but exclude cache keys
            for key, value in record.annotations.items():
                if not key.startswith("_biotooler_"):
                    orf_record.annotations[key] = value

            # Process each positional feature
            for _feat_name, feat_fn in positional_features:
                # Get position space and vector keys
                position_space = feat_fn.position_space  # type: ignore[union-attr]
                vector_keys = feat_fn.vector_keys  # type: ignore[union-attr]

                # Compute per-position vectors for the full ORF
                vectors = feat_fn.compute_vector(orf_record)  # type: ignore[union-attr]

                # Generate window boundaries using the shared indexing helper
                # Map PositionSpace enum to string for the helper function
                position_space_str = (
                    "codon" if position_space == PositionSpace.CODON else "residue"
                )

                # Get window boundaries in nucleotide space
                window_boundaries = compute_window_indices(
                    sequence_length=orf_len,
                    window_size=window_nt,
                    step=step_nt,
                    drop_partial=drop_partial,
                    position_space=position_space_str,
                    start_offset=0,
                )

                # Process each window
                for window_start_nt, window_end_nt in window_boundaries:
                    # Convert nucleotide positions to position space indices
                    if position_space == PositionSpace.CODON:
                        # For codon space, convert nt positions to codon indices
                        window_start_pos = window_start_nt // 3
                        window_end_pos = window_end_nt // 3
                    else:  # RESIDUE
                        # For residue space, positions are the same as nt positions
                        window_start_pos = window_start_nt
                        window_end_pos = window_end_nt

                    # Aggregate each vector key for this window
                    for key, agg_spec in vector_keys.items():
                        if key not in vectors:
                            raise ValueError(
                                f"Feature compute_vector did not return expected key '{key}'"
                            )

                        vector = vectors[key]
                        window_values = vector[window_start_pos:window_end_pos]

                        # Apply aggregation function
                        aggregated_value = agg_spec.aggregation_fn(window_values)

                        # Store with column name format: {FAMILY}_{key}_{AGG}_{window_start_nt}
                        col_name = f"{self.name.upper()}_{key}_{agg_spec.name}_{window_start_nt}"
                        feature_data[col_name] = aggregated_value

        # Process non-positional features using legacy path
        if non_positional_features:
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

            if windows:
                # Temporarily store features for legacy path
                original_features = self._features
                self._features = dict(non_positional_features)

                # Compute features over windows using shared helper
                self._compute_features_over_windows(
                    windows=windows,
                    feature_data=feature_data,
                    window_start_key="window_start",
                    init_state_kwargs={"record": record, "orf": resolved_orf},
                )

                # Restore original features
                self._features = original_features

        # Format as wide DataFrame using shared helper
        return self._format_wide_dataframe(
            feature_data=feature_data,
            metadata_cols=["record_id", "orf_start", "orf_end"],
        )

    def compute_orf_windows_v2(
        self,
        record: SeqRecord,
        *,
        orf: OrfSpan | None = None,
        orf_index: int | None = None,
        window_nt: int,
        step_nt: int,
        drop_partial: bool = True,
    ) -> pd.DataFrame:
        """Compute features for ORF windows using v2 windowing semantics (positional features).

        This method implements the new windowing approach for positional features:
        1. Compute per-position values on the full ORF sequence (via compute_vector)
        2. Aggregate per-position values into windows using specified aggregation functions

        Unlike compute_orf_windows (legacy), this method:
        - Does NOT slice the sequence into windows
        - Computes per-position values once for the entire ORF
        - Aggregates these values into windows based on the feature's position space

        This method only works with features implementing the PositionalFeature protocol,
        which requires:
        - position_space property (RESIDUE or CODON)
        - vector_keys property (dict mapping keys to AggregationSpec)
        - compute_vector method (returns dict of numpy arrays)

        ORF Resolution Rules (same as compute_orf_windows):
        - If `orf` is provided: use it directly
        - Else if `orf_index` is provided: find candidates and select that index
        - Else: try to get attached ORF from record; if missing, raise clear error

        Output Format (same as compute_orf_windows):
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
                       or if neither orf, orf_index, nor attached ORF is available,
                       or if feature does not implement PositionalFeature protocol
            KeyError: If trying to use attached ORF but none exists

        Examples:
            >>> from Bio.Seq import Seq
            >>> from Bio.SeqRecord import SeqRecord
            >>> import numpy as np
            >>> from biotooler.features.aggregation import AggregationSpec, PositionSpace
            >>> # Define a positional feature
            >>> class GCFeature:
            ...     @property
            ...     def position_space(self):
            ...         return PositionSpace.RESIDUE
            ...     @property
            ...     def vector_keys(self):
            ...         return {"gc": AggregationSpec(name="MEAN", aggregation_fn=np.mean)}
            ...     def compute_vector(self, record, **kwargs):
            ...         seq = str(record.seq).upper()
            ...         gc_vector = np.array([1.0 if b in 'GC' else 0.0 for b in seq])
            ...         return {"gc": gc_vector}
            >>> fs = FeatureSet(GCFeature(), name="gc_content")
            >>> record = SeqRecord(Seq("ATGAAACCCGGGTTT"), id="test")
            >>> result = fs.compute_orf_windows_v2(
            ...     record, orf=(0, 15), window_nt=9, step_nt=3
            ... )
            >>> result.shape
            (1, 6)
            >>> 'gc_content.gc_0' in result.columns
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

        # Resolve ORF coordinates (same logic as compute_orf_windows)
        resolved_orf: OrfSpan
        if orf is not None:
            resolved_orf = orf
        elif orf_index is not None:
            candidates = find_orf_candidates(record)
            resolved_orf = select_orf_by_index(candidates, orf_index)
        else:
            try:
                resolved_orf = get_orf(record)
            except KeyError as e:
                raise ValueError(
                    f"No ORF information provided for record {record.id!r}. "
                    "Please provide 'orf' parameter, 'orf_index' parameter, "
                    "or attach an ORF using attach_orf()."
                ) from e

        # Validate step_nt (must be multiple of 3 for codon alignment)
        if step_nt % 3 != 0:
            raise ValueError(
                f"step_nt must be a multiple of 3 for codon-aligned windows, got {step_nt}"
            )

        # Build feature data dictionary for single row
        feature_data: dict[str, Any] = {
            "record_id": record.id,
            "orf_start": resolved_orf[0],
            "orf_end": resolved_orf[1],
        }

        # Extract ORF sequence
        orf_start, orf_end = resolved_orf
        seq_str = get_seq_str(record)
        orf_seq_str = seq_str[orf_start:orf_end]
        orf_len = len(orf_seq_str)

        # Create ORF record for compute_vector
        orf_record = SeqRecord(
            Seq(orf_seq_str),
            id=record.id,
            description=record.description,
        )
        # Copy annotations but exclude cache keys
        for key, value in record.annotations.items():
            if not key.startswith("_biotooler_"):
                orf_record.annotations[key] = value

        # Process each feature
        for _feat_name, feat_fn in self._features.items():
            # Check if feature implements PositionalFeature protocol (duck-typing)
            has_positional = (
                hasattr(feat_fn, "position_space")
                and hasattr(feat_fn, "vector_keys")
                and hasattr(feat_fn, "compute_vector")
            )

            if not has_positional:
                raise ValueError(
                    "Feature does not implement PositionalFeature protocol. "
                    "compute_orf_windows_v2 requires features with position_space, "
                    "vector_keys, and compute_vector. Use compute_orf_windows for "
                    "non-positional features."
                )

            # Get position space and vector keys
            position_space = feat_fn.position_space  # type: ignore[union-attr]
            vector_keys = feat_fn.vector_keys  # type: ignore[union-attr]

            # Generate window boundaries using the shared indexing helper
            # Map PositionSpace enum to string for the helper function
            position_space_str = (
                "codon" if position_space == PositionSpace.CODON else "residue"
            )

            # Get window boundaries in nucleotide space
            window_boundaries = compute_window_indices(
                sequence_length=orf_len,
                window_size=window_nt,
                step=step_nt,
                drop_partial=drop_partial,
                position_space=position_space_str,
                start_offset=0,
            )

            # Compute union of all window positions in position space
            # This optimization allows features to compute only needed positions
            union_positions_set = set()
            window_position_ranges = []  # Store (start_nt, start_pos, end_pos) for each window

            for window_start_nt, window_end_nt in window_boundaries:
                # Convert nucleotide positions to position space indices
                if position_space == PositionSpace.CODON:
                    # For codon space, convert nt positions to codon indices
                    window_start_pos = window_start_nt // 3
                    window_end_pos = window_end_nt // 3
                else:  # RESIDUE
                    # For residue space, positions are the same as nt positions
                    window_start_pos = window_start_nt
                    window_end_pos = window_end_nt

                window_position_ranges.append((window_start_nt, window_start_pos, window_end_pos))
                union_positions_set.update(range(window_start_pos, window_end_pos))

            # Convert to sorted array for compute_vector
            union_positions = np.array(sorted(union_positions_set), dtype=np.int64)

            # Compute per-position vectors
            # Pass union_positions to allow features to optimize (compute only needed positions)
            # Features that don't support this optimization can ignore it and compute all positions
            vectors = feat_fn.compute_vector(orf_record, positions=union_positions)  # type: ignore[union-attr]

            # Check if feature returned sparse (positions-indexed) or full vector
            # If all union positions are consecutive starting from 0, it's full vector
            is_full_vector = (
                len(union_positions) > 0
                and union_positions[0] == 0
                and len(union_positions) == union_positions[-1] + 1
            )

            # Process each window and aggregate
            for window_start_nt, window_start_pos, window_end_pos in window_position_ranges:
                # Aggregate each vector key for this window
                for key, agg_spec in vector_keys.items():
                    if key not in vectors:
                        raise ValueError(
                            f"Feature compute_vector did not return expected key '{key}'"
                        )

                    vector = vectors[key]

                    if is_full_vector:
                        # Full vector: slice directly by positions
                        window_values = vector[window_start_pos:window_end_pos]
                    else:
                        # Sparse vector: use searchsorted to find indices in union_positions
                        idx_start = np.searchsorted(union_positions, window_start_pos)
                        idx_end = np.searchsorted(union_positions, window_end_pos)
                        window_values = vector[idx_start:idx_end]

                    # Apply aggregation function
                    aggregated_value = agg_spec.aggregation_fn(window_values)

                    # Store with column name format: {FAMILY}_{key}_{AGG}_{window_start_nt}
                    col_name = f"{self.name.upper()}_{key}_{agg_spec.name}_{window_start_nt}"
                    feature_data[col_name] = aggregated_value


        # Format as wide DataFrame using shared helper
        return self._format_wide_dataframe(
            feature_data=feature_data,
            metadata_cols=["record_id", "orf_start", "orf_end"],
        )

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
            # Cache keys like '_biotooler_seq_str' are sequence-specific and would be
            # invalid for the region subsequence, causing iter_windows to use wrong data
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

        # Compute features over windows using shared helper
        self._compute_features_over_windows(
            windows=windows,
            feature_data=feature_data,
            window_start_key="start",
            init_state_kwargs={
                "record": region_record,
                "region": (region_start, region_end),
            },
        )

        # Format as wide DataFrame using shared helper
        return self._format_wide_dataframe(
            feature_data=feature_data,
            metadata_cols=["record_id", "region_start", "region_end"],
        )
