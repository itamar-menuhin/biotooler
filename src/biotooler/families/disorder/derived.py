"""Derived disorder features computed from DISORDER_P predictions.

This module provides features that compute summary statistics from existing
DISORDER_P (disorder probability) vectors without re-running predictors.
"""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.types import Scalar


class DisorderDerivedScalars:
    """Compute scalar summary features from DISORDER_P disorder predictions.

    This feature computes summary statistics from per-residue disorder
    probabilities (DISORDER_P) that have been previously computed by a
    disorder predictor like DisorderProfileMetapredict. It does not call
    the predictor itself.

    The feature computes:
    - DISORDER_FRAC: Fraction of residues with disorder probability >= threshold
    - DISORDER_LONGEST_IDR: Length of longest contiguous region with disorder >= threshold
    - DISORDER_MEAN: Mean disorder probability across all residues
    - DISORDER_P95: 95th percentile of disorder probabilities

    Args:
        threshold: Disorder probability threshold for determining disordered residues
            (default: 0.5). Used for DISORDER_FRAC and DISORDER_LONGEST_IDR.

    Examples:
        >>> import numpy as np
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = DisorderDerivedScalars(threshold=0.5)
        >>> record = SeqRecord(Seq("MKALVSWGR"), id="test")
        >>> # Assume DISORDER_P has been computed and cached
        >>> disorder_p = [0.1, 0.2, 0.6, 0.7, 0.8, 0.3, 0.2, 0.1, 0.9]
        >>> record.annotations["DISORDER_P"] = np.array(disorder_p)
        >>> result = feature(record)
        >>> result["DISORDER_FRAC"]  # 4/9 residues above 0.5
        0.444...
        >>> result["DISORDER_LONGEST_IDR"]  # Longest run is positions 2-4 (length 3)
        3
        >>> result["DISORDER_MEAN"]  # Mean of all probabilities
        0.433...

    Notes:
        - Requires DISORDER_P to be present in record.annotations
        - DISORDER_P should be a numpy array of disorder probabilities
        - Empty sequences return 0.0 for DISORDER_FRAC and DISORDER_LONGEST_IDR,
          and NaN for DISORDER_MEAN and DISORDER_P95
        - DISORDER_LONGEST_IDR counts the number of consecutive residues >= threshold
    """

    def __init__(self, *, threshold: float = 0.5):
        """Initialize DisorderDerivedScalars.

        Args:
            threshold: Disorder probability threshold for determining disordered
                residues (default: 0.5). Must be in range [0, 1].

        Raises:
            ValueError: If threshold is not in range [0, 1]
        """
        if not 0.0 <= threshold <= 1.0:
            raise ValueError(f"threshold must be in range [0, 1], got {threshold}")
        self.threshold = threshold

    def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
        """Compute derived disorder features from DISORDER_P in record annotations.

        Args:
            record: SeqRecord with DISORDER_P in annotations. DISORDER_P should be
                a numpy array or array-like of disorder probabilities.

        Returns:
            Dictionary mapping feature names to scalar values:
            - DISORDER_FRAC: Fraction of residues >= threshold
            - DISORDER_LONGEST_IDR: Length of longest contiguous region >= threshold
            - DISORDER_MEAN: Mean disorder probability
            - DISORDER_P95: 95th percentile of disorder probabilities

        Raises:
            ValueError: If DISORDER_P is not found in record.annotations or
                if DISORDER_P is not a valid array
        """
        # Get DISORDER_P from annotations
        if "DISORDER_P" not in record.annotations:
            raise ValueError(
                f"DISORDER_P not found in annotations for record {record.id!r}. "
                "Ensure DisorderProfileMetapredict or equivalent has been run first."
            )

        disorder_p = record.annotations["DISORDER_P"]

        # Convert to numpy array if needed and validate
        try:
            disorder_p = np.asarray(disorder_p, dtype=np.float64)
        except (ValueError, TypeError) as e:
            raise ValueError(
                f"DISORDER_P in record {record.id!r} could not be converted to array: {e}"
            ) from e

        # Validate shape
        if disorder_p.ndim != 1:
            raise ValueError(
                f"DISORDER_P in record {record.id!r} must be 1-dimensional, "
                f"got shape {disorder_p.shape}"
            )

        # Handle empty sequence
        if len(disorder_p) == 0:
            return {
                "DISORDER_FRAC": 0.0,
                "DISORDER_LONGEST_IDR": 0,
                "DISORDER_MEAN": float("nan"),
                "DISORDER_P95": float("nan"),
            }

        # Compute DISORDER_FRAC: fraction of residues >= threshold
        disordered_mask = disorder_p >= self.threshold
        disorder_frac = float(np.mean(disordered_mask))

        # Compute DISORDER_LONGEST_IDR: longest contiguous run >= threshold
        longest_idr = self._compute_longest_run(disordered_mask)

        # Compute DISORDER_MEAN: mean disorder probability
        disorder_mean = float(np.mean(disorder_p))

        # Compute DISORDER_P95: 95th percentile
        disorder_p95 = float(np.percentile(disorder_p, 95))

        return {
            "DISORDER_FRAC": disorder_frac,
            "DISORDER_LONGEST_IDR": longest_idr,
            "DISORDER_MEAN": disorder_mean,
            "DISORDER_P95": disorder_p95,
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
