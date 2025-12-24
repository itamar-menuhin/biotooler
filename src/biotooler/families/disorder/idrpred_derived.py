"""Derived IDRPred features computed from IDRPRED_IDR mask.

This module provides features that compute summary statistics from existing
IDRPRED_IDR (binary IDR membership mask) vectors without re-running IDRPred.
"""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar


class IDRPredDerivedScalars:
    """Compute scalar summary features from IDRPRED_IDR mask predictions.

    This feature computes summary statistics from per-residue IDR membership
    mask (IDRPRED_IDR) that has been previously computed by IDRPredConsensusMask.
    It does not call the IDRPred CLI itself.

    The feature computes:
    - IDRPRED_FRAC_IDR: Fraction of residues in IDRs (mean of mask values)
    - IDRPRED_LONGEST_IDR_LEN: Length of longest contiguous IDR region
    - IDRPRED_NUM_IDR_SEGMENTS: Number of distinct contiguous IDR segments

    No parameters are needed since IDRPRED_IDR is already a binary mask (0 or 1).

    Examples:
        >>> import numpy as np
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = IDRPredDerivedScalars()
        >>> record = SeqRecord(Seq("MKALVSWGR"), id="test")
        >>> # Assume IDRPRED_IDR has been computed and cached
        >>> idrpred_idr = [0, 0, 1, 1, 1, 0, 0, 0, 1]
        >>> record.annotations["IDRPRED_IDR"] = np.array(idrpred_idr)
        >>> result = feature(record)
        >>> result["IDRPRED_FRAC_IDR"]  # 4/9 residues in IDRs
        0.444...
        >>> result["IDRPRED_LONGEST_IDR_LEN"]  # Longest run is positions 2-4 (length 3)
        3
        >>> result["IDRPRED_NUM_IDR_SEGMENTS"]  # Two segments: [2-4] and [8]
        2

    Notes:
        - Requires IDRPRED_IDR to be present in record.annotations
        - IDRPRED_IDR should be a numpy array of binary mask values (0 or 1)
        - Empty sequences return 0 for all metrics
        - All-zero mask returns 0 for all metrics
        - All-one mask returns frac=1.0, longest=length, segments=1
    """

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute derived IDRPred features from IDRPRED_IDR in record annotations.

        Args:
            record: SeqRecord with IDRPRED_IDR in annotations. IDRPRED_IDR should be
                a numpy array or array-like of binary mask values (0 or 1).

        Returns:
            Dictionary mapping feature names to scalar values:
            - IDRPRED_FRAC_IDR: Fraction of residues in IDRs (mean of mask)
            - IDRPRED_LONGEST_IDR_LEN: Length of longest contiguous IDR region
            - IDRPRED_NUM_IDR_SEGMENTS: Number of distinct contiguous IDR segments

        Raises:
            ValueError: If IDRPRED_IDR is not found in record.annotations or
                if IDRPRED_IDR is not a valid array
        """
        # Get IDRPRED_IDR from annotations
        if "IDRPRED_IDR" not in record.annotations:
            raise ValueError(
                f"IDRPRED_IDR not found in annotations for record {record.id!r}. "
                "Ensure IDRPredConsensusMask or equivalent has been run first."
            )

        idrpred_idr = record.annotations["IDRPRED_IDR"]

        # Convert to numpy array if needed and validate
        try:
            idrpred_idr = np.asarray(idrpred_idr, dtype=np.float64)
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"IDRPRED_IDR in record {record.id!r} could not be converted to array: {e}"
            ) from e

        # Validate shape
        if idrpred_idr.ndim != 1:
            raise ValueError(
                f"IDRPRED_IDR in record {record.id!r} must be 1-dimensional, "
                f"got shape {idrpred_idr.shape}"
            )

        # Handle empty sequence
        if len(idrpred_idr) == 0:
            return {
                "IDRPRED_FRAC_IDR": 0.0,
                "IDRPRED_LONGEST_IDR_LEN": 0,
                "IDRPRED_NUM_IDR_SEGMENTS": 0,
            }

        # Convert to boolean mask for easier processing
        mask = idrpred_idr.astype(bool)

        # Compute IDRPRED_FRAC_IDR: fraction of residues in IDRs (mean of mask)
        frac_idr = float(np.mean(mask))

        # Compute IDRPRED_LONGEST_IDR_LEN: longest contiguous run
        longest_idr_len = self._compute_longest_run(mask)

        # Compute IDRPRED_NUM_IDR_SEGMENTS: number of contiguous segments
        num_idr_segments = self._compute_segment_count(mask)

        return {
            "IDRPRED_FRAC_IDR": frac_idr,
            "IDRPRED_LONGEST_IDR_LEN": longest_idr_len,
            "IDRPRED_NUM_IDR_SEGMENTS": num_idr_segments,
        }

    def _compute_longest_run(self, mask: np.ndarray) -> int:
        """Compute the length of the longest contiguous run of True values.

        Args:
            mask: Boolean numpy array

        Returns:
            Length of longest contiguous run of True values. Returns 0 if
            no True values are present.
        """
        if len(mask) == 0 or not np.any(mask):
            return 0

        # Find transitions: where mask changes from False to True or True to False
        # Prepend False and append False to handle edge cases
        padded = np.concatenate(([False], mask, [False]))
        # Find where transitions occur
        diff = np.diff(padded.astype(int))
        # Starts of runs: where diff == 1 (False -> True)
        starts = np.where(diff == 1)[0]
        # Ends of runs: where diff == -1 (True -> False)
        ends = np.where(diff == -1)[0]

        # Compute run lengths
        run_lengths = ends - starts

        # Return maximum run length
        return int(np.max(run_lengths))

    def _compute_segment_count(self, mask: np.ndarray) -> int:
        """Compute the number of contiguous segments of True values.

        Args:
            mask: Boolean numpy array

        Returns:
            Number of contiguous segments of True values. Returns 0 if
            no True values are present.
        """
        if len(mask) == 0 or not np.any(mask):
            return 0

        # Find transitions: where mask changes from False to True or True to False
        # Prepend False and append False to handle edge cases
        padded = np.concatenate(([False], mask, [False]))
        # Find where transitions occur
        diff = np.diff(padded.astype(int))
        # Starts of runs: where diff == 1 (False -> True)
        starts = np.where(diff == 1)[0]

        # Number of segments equals number of starts
        return int(len(starts))
