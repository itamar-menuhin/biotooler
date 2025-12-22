"""Codon usage bias feature computation using codonbias package."""

from collections.abc import Sequence

import codonbias.scores
import codonbias.stats
import codonbias.utils
from Bio.SeqRecord import SeqRecord

from biotooler.core.seq_utils import get_seq_str


class CodonBiasFeature:
    """Feature that computes codon usage bias scores using codonbias models.

    This feature wraps external codonbias.scores.ScalarScore instances (e.g., CAI, ENC,
    FOP, RSCU, RCBS/DCBS, tAI, nTE, CPB/CPS) and computes scalar features on DNA/RNA
    regions.

    The feature supports two computation modes:
    1. Baseline: Call model.get_score(seq_str, slice=slice(start, end)) for each window
    2. Rolling (incremental): Maintain rolling codon counts and recompute using
       codonbias utilities for each step (when weights are accessible)

    For models where weights are not accessible via public attributes/methods, rolling
    mode automatically falls back to baseline per-window get_score(slice=...).

    Args:
        models: List of codonbias.scores.ScalarScore instances to compute
        names: Optional list of names for each model (defaults to model class names)

    Examples:
        >>> from codonbias.scores import CodonAdaptationIndex, EffectiveNumberOfCodons
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> # Create models
        >>> ref_seq = "ATGATGATGATGATG"
        >>> cai = CodonAdaptationIndex(ref_seq)
        >>> enc = EffectiveNumberOfCodons()
        >>> # Create feature
        >>> feature = CodonBiasFeature([cai, enc], names=["CAI", "ENC"])
        >>> # Compute on a sequence
        >>> record = SeqRecord(Seq("ATGATGATGATG"), id="test")
        >>> result = feature(record)
        >>> print(result.keys())
        dict_keys(['CAI', 'ENC'])
    """

    def __init__(
        self,
        models: Sequence[codonbias.scores.ScalarScore],
        *,
        names: Sequence[str] | None = None,
    ):
        """Initialize CodonBiasFeature.

        Args:
            models: Sequence of codonbias.scores.ScalarScore instances
            names: Optional sequence of names for each model (must match length of models)

        Raises:
            ValueError: If names is provided but length doesn't match models
        """
        self.models = list(models)
        if names is None:
            self.names = [type(model).__name__ for model in self.models]
        else:
            if len(names) != len(self.models):
                raise ValueError(
                    f"Length of names ({len(names)}) must match length of models "
                    f"({len(self.models)})"
                )
            self.names = list(names)

    def __call__(self, record: SeqRecord) -> dict[str, float]:
        """Compute codon bias scores for the entire sequence.

        Args:
            record: DNA or RNA SeqRecord to analyze

        Returns:
            Dictionary mapping feature names to scalar values

        Raises:
            ValueError: If record has protein alphabet
        """
        # Validate alphabet
        _validate_alphabet(record)

        # Get sequence and convert RNA to DNA
        seq_str = get_seq_str(record)
        seq_str = _convert_rna_to_dna(seq_str)

        # Compute scores for each model
        result = {}
        for name, model in zip(self.names, self.models):
            score = model.get_score(seq_str)
            result[name] = float(score)

        return result

    def init_state(
        self,
        record: SeqRecord,
        *,
        orf: tuple[int, int],
        window_start: int,
        window_end: int,
        **kwargs,
    ) -> dict:
        """Initialize state for the first window.

        Uses codonbias.stats.CodonCounter to initialize codon counts for the window.

        Args:
            record: The SeqRecord containing the sequence
            orf: Tuple of (start, end) coordinates for the ORF
            window_start: Start position of window relative to ORF start
            window_end: End position of window relative to ORF start
            **kwargs: Additional parameters

        Returns:
            State dictionary containing record, orf, codon counts, and per-model data
        """
        # Validate alphabet
        _validate_alphabet(record)

        # Extract window sequence
        orf_start, orf_end = orf
        abs_start = orf_start + window_start
        abs_end = orf_start + window_end
        seq_str = get_seq_str(record)[abs_start:abs_end]
        seq_str = _convert_rna_to_dna(seq_str)

        # Initialize codon counter
        counter = codonbias.stats.CodonCounter(seq_str)
        codon_counts = counter.get_codon_table()

        # Try to get weights for each model and determine if incremental is possible
        model_data = []
        for model in self.models:
            weights_info = _try_get_weights(model)
            model_data.append(weights_info)

        state = {
            "record": record,
            "orf": orf,
            "codon_counts": codon_counts,
            "model_data": model_data,
            "window_start": window_start,
            "window_end": window_end,
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
        """Update state for the next window by updating codon counts.

        Args:
            state: The state object to update
            out_start: Start position (relative to ORF) of bases leaving the window
            out_end: End position (relative to ORF) of bases leaving the window
            in_start: Start position (relative to ORF) of bases entering the window
            in_end: End position (relative to ORF) of bases entering the window
            **kwargs: Additional parameters
        """
        record = state["record"]
        orf = state["orf"]
        orf_start = orf[0]
        codon_counts = state["codon_counts"]

        # Update window position tracking
        # New window: from out_end to in_end (sliding the window forward)
        state["window_start"] = out_end
        state["window_end"] = in_end

        # Remove outgoing codons
        if out_end > out_start:
            abs_out_start = orf_start + out_start
            abs_out_end = orf_start + out_end
            out_seq = get_seq_str(record)[abs_out_start:abs_out_end]
            out_seq = _convert_rna_to_dna(out_seq)
            out_counter = codonbias.stats.CodonCounter(out_seq)
            out_counts = out_counter.get_codon_table()
            codon_counts -= out_counts

        # Add incoming codons
        if in_end > in_start:
            abs_in_start = orf_start + in_start
            abs_in_end = orf_start + in_end
            in_seq = get_seq_str(record)[abs_in_start:abs_in_end]
            in_seq = _convert_rna_to_dna(in_seq)
            in_counter = codonbias.stats.CodonCounter(in_seq)
            in_counts = in_counter.get_codon_table()
            codon_counts += in_counts

    def emit(self, state: dict) -> dict[str, float]:
        """Emit feature values from current state.

        For models with accessible weights, computes scores using codonbias.utils
        functions. For models without accessible weights, falls back to baseline
        get_score() on the current window sequence.

        Args:
            state: The state object

        Returns:
            Dictionary mapping feature names to scalar values
        """
        result = {}
        codon_counts = state["codon_counts"]
        model_data = state["model_data"]

        for name, model, weights_info in zip(self.names, self.models, model_data):
            if weights_info["has_weights"]:
                # Use incremental computation with weights
                weights = weights_info["weights"]
                weight_type = weights_info["type"]

                if weight_type == "log":
                    # Use geometric mean for log weights
                    score = codonbias.utils.geomean(weights, codon_counts)
                else:
                    # Use arithmetic mean for regular weights
                    score = codonbias.utils.mean(weights, codon_counts)

                result[name] = float(score)
            else:
                # Fall back to baseline: get window sequence and call get_score
                record = state["record"]
                orf = state["orf"]
                orf_start = orf[0]
                window_start = state["window_start"]
                window_end = state["window_end"]

                abs_start = orf_start + window_start
                abs_end = orf_start + window_end
                seq_str = get_seq_str(record)[abs_start:abs_end]
                seq_str = _convert_rna_to_dna(seq_str)
                score = model.get_score(seq_str)
                result[name] = float(score)

        return result


def _validate_alphabet(record: SeqRecord) -> None:
    """Validate that record is DNA or RNA, not protein.

    Args:
        record: SeqRecord to validate

    Raises:
        ValueError: If record has protein molecule_type annotation
    """
    if "molecule_type" in record.annotations:
        mol_type = record.annotations["molecule_type"]
        if isinstance(mol_type, str) and mol_type.upper() == "PROTEIN":
            raise ValueError(
                "Codon bias features apply only to DNA/RNA sequences, not protein sequences"
            )


def _convert_rna_to_dna(seq_str: str) -> str:
    """Convert RNA sequence (U) to DNA (T).

    Args:
        seq_str: Sequence string (may contain U/u)

    Returns:
        Sequence string with U/u replaced by T/t
    """
    return seq_str.replace("U", "T").replace("u", "t")


def _try_get_weights(
    model: codonbias.scores.ScalarScore,
) -> dict:
    """Try to extract weights from a codonbias model.

    Attempts to get weights in this order:
    1. model.log_weights (preferred) -> use with geomean
    2. model.weights -> use with mean
    3. model.get_weights() if exists -> use with mean
    4. None -> fallback to baseline mode

    Args:
        model: A codonbias ScalarScore model

    Returns:
        Dictionary with keys:
        - has_weights: bool indicating if weights were found
        - weights: pandas Series of weights (or None)
        - type: "log" or "linear" (or None)
    """
    # Try log_weights first
    if hasattr(model, "log_weights"):
        return {
            "has_weights": True,
            "weights": model.log_weights,
            "type": "log",
        }

    # Try weights
    if hasattr(model, "weights"):
        return {
            "has_weights": True,
            "weights": model.weights,
            "type": "linear",
        }

    # Try get_weights() method
    if hasattr(model, "get_weights") and callable(model.get_weights):
        try:
            weights = model.get_weights()
            # Check if it returns a valid weights structure
            if weights is not None:
                return {
                    "has_weights": True,
                    "weights": weights,
                    "type": "linear",
                }
        except Exception:
            # get_weights() may require arguments or fail for other reasons
            pass

    # No accessible weights - must use fallback
    return {
        "has_weights": False,
        "weights": None,
        "type": None,
    }
