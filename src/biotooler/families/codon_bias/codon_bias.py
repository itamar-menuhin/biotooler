"""Codon usage bias feature computation using codonbias package."""

import hashlib
import inspect
import json
from collections import OrderedDict
from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.lazy_import import lazy_import
from biotooler.core.seq_utils import get_seq_str
from biotooler.features.aggregation import AggregationSpec, PositionSpace

if TYPE_CHECKING:
    from biotooler.core.reference_sequences import ReferenceSequenceSet

# Lazy import codonbias modules
codonbias = lazy_import(  # type: ignore[misc]
    "codonbias", extra="codon_bias", purpose="computing codon usage bias features"
)


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
        max_cache_size: Maximum number of model sets to cache (default: 8)

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
        models: Sequence[Any],  # codonbias.scores.ScalarScore
        *,
        names: Sequence[str] | None = None,
        max_cache_size: int = 8,
    ):
        """Initialize CodonBiasFeature.

        Args:
            models: Sequence of codonbias.scores.ScalarScore instances
            names: Optional sequence of names for each model (must match length of models)
            max_cache_size: Maximum number of model sets to cache (default: 8)

        Raises:
            ValueError: If names is provided but length doesn't match models
        """
        # Ensure we have lists for consistent iteration
        self.models = list(models) if not isinstance(models, list) else models
        if names is None:
            self.names = [type(model).__name__ for model in self.models]
        else:
            if len(names) != len(self.models):
                raise ValueError(
                    f"Length of names ({len(names)}) must match length of models "
                    f"({len(self.models)})"
                )
            self.names = list(names) if not isinstance(names, list) else names

        # Initialize bounded cache for model building
        self._model_cache: OrderedDict[str, list[Any]] = OrderedDict()
        self._max_cache_size = max_cache_size

    @property
    def position_space(self) -> PositionSpace:
        """The position space for this feature (CODON).

        Returns:
            PositionSpace.CODON indicating features are computed per codon
        """
        return PositionSpace.CODON

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Mapping of feature keys to aggregation specifications.

        For CodonBiasFeature, this always returns an empty dict to ensure
        compute_orf_windows uses the incremental path (init_state/update/emit).
        This provides consistent column naming without aggregation suffixes.

        The positional interface (compute_vector) is still available for
        explicit use with compute_orf_windows_v2 when needed.

        Returns:
            Empty dictionary to force incremental windowing path
        """
        # Return empty dict to force incremental path in compute_orf_windows.
        # This ensures column names like CB_CAI_0 (without aggregation suffix)
        # instead of CB_CAI_GEOMEAN_0 for consistency with legacy behavior.
        return {}

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-codon feature values for the entire sequence.

        This method uses the codonbias score's get_vector() method to compute
        per-codon values across the full sequence context, without any slicing.
        This ensures that each positional value is computed with full sequence
        context information.

        Args:
            record: The SeqRecord containing the sequence
            **kwargs: Additional parameters that may include:
                - positions: Optional array of codon indices. For codon bias features,
                  this is currently ignored and all positions are computed, as the
                  codonbias models require full sequence context. The windowing code
                  will select the needed positions from the full vector.

        Returns:
            Dictionary mapping feature names to numpy arrays of per-codon values.
            Only scores that support get_vector() are included.
        """
        # Note: positions parameter is intentionally ignored for codon_bias features
        # because they require full sequence context for accurate computation

        # Validate alphabet
        _validate_alphabet(record)

        # Get sequence and convert RNA to DNA
        seq_str = get_seq_str(record)
        seq_str = _convert_rna_to_dna(seq_str)

        vectors = {}
        for name, model in zip(self.names, self.models, strict=True):
            # Check if model supports get_vector
            has_get_vector = hasattr(model, "get_vector") and callable(
                getattr(model, "get_vector", None)
            )

            if not has_get_vector:
                # Skip models without get_vector (they use legacy incremental path)
                continue

            # Use get_vector to compute per-codon values with full context
            vector = model.get_vector(seq_str)
            vectors[name] = vector

        return vectors

    @classmethod
    def from_reference(
        cls,
        reference_set: "ReferenceSequenceSet",
        scores: Sequence[str | type | Any],
        *,
        names: Sequence[str] | None = None,
        score_kwargs: dict[str, dict[str, Any]] | None = None,
        model_cache: OrderedDict[str, list[Any]] | None = None,
        max_cache_size: int = 8,
    ) -> "CodonBiasFeature":
        """Create CodonBiasFeature from a ReferenceSequenceSet.

        This factory method builds codon bias score models from a shared reference
        sequence set, automatically handling reference sequence requirements.

        Models are cached to avoid expensive rebuilds when the same reference set
        and parameters are used repeatedly. The cache can be shared across multiple
        from_reference calls by passing the same OrderedDict instance.

        Args:
            reference_set: ReferenceSequenceSet containing CDS sequences
            scores: Sequence of score identifiers, which can be:
                - Abbreviations like "CAI", "ENC", "FOP"
                - Class names like "CodonAdaptationIndex"
                - Class objects like CodonAdaptationIndex
                - Already instantiated score objects
            names: Optional custom names for the scores (defaults to abbreviations or class names)
            score_kwargs: Optional dict mapping score identifiers to their constructor kwargs
                Example: {"CAI": {"genetic_code": 11}, "ENC": {"bg_correction": True}}
            model_cache: Optional cache dictionary for reusing models (feature-instance scoped).
                If None, a new cache is created for this call.
            max_cache_size: Maximum number of model sets to cache (default: 8)

        Returns:
            CodonBiasFeature instance with instantiated score models

        Raises:
            ValueError: If reference_set has no CDS sequences when required by a score
            ValueError: If a score identifier cannot be resolved
            TypeError: If a score constructor fails

        Examples:
            >>> from biotooler.core.reference_sequences import ReferenceSequenceSet
            >>> ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATG", "gene2": "ATGATGATG"})
            >>> # Using abbreviations
            >>> feature = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"])
            >>> # Using class objects with kwargs
            >>> from codonbias.scores import CodonAdaptationIndex
            >>> feature = CodonBiasFeature.from_reference(
            ...     ref_set,
            ...     [CodonAdaptationIndex, "ENC"],
            ...     names=["CAI", "ENC"],
            ...     score_kwargs={"CAI": {"genetic_code": 11}}
            ... )
            >>> # Reusing cache across multiple calls
            >>> from collections import OrderedDict
            >>> cache = OrderedDict()
            >>> feature1 = CodonBiasFeature.from_reference(ref_set, ["CAI"], model_cache=cache)
            >>> feature2 = CodonBiasFeature.from_reference(ref_set, ["CAI"], model_cache=cache)
            >>> # feature2 reuses models built for feature1
        """
        # Check if reference_set has CDS sequences
        if not reference_set.cds:
            raise ValueError(
                "ReferenceSequenceSet must contain CDS sequences to build codon bias models"
            )

        # Initialize or use provided cache
        if model_cache is None:
            model_cache = OrderedDict()

        # Generate cache key
        cache_key = _generate_cache_key(reference_set, scores, score_kwargs)

        # Check if models are already cached
        if cache_key in model_cache:
            # Move to end for LRU behavior
            model_cache.move_to_end(cache_key)
            models = model_cache[cache_key]
            # Resolve names for cached models
            resolved_names = [type(model).__name__ for model in models]
        else:
            # Build models since not in cache
            # Get concatenated reference sequences for models that need ref_seq
            ref_seq_strings = reference_set.cds_strings()
            concatenated_ref_seq = "".join(ref_seq_strings)

            # Resolve score identifiers and instantiate models
            models = []
            resolved_names = []

            for score_id in scores:
                # Resolve the score class
                score_class = _resolve_score_identifier(score_id)

                # Check if it's already an instance
                if not inspect.isclass(score_class):
                    # Already instantiated
                    models.append(score_class)
                    resolved_names.append(type(score_class).__name__)
                    continue

                # Get kwargs for this score if provided
                kwargs = {}
                if score_kwargs:
                    # Try to match by various identifiers in order of specificity:
                    # 1. Original score_id (as provided by user)
                    # 2. Class name (e.g., "CodonAdaptationIndex")
                    # 3. Abbreviation (e.g., "CAI")
                    for key in [
                        score_id,
                        score_class.__name__,
                        _get_score_abbreviation(score_class.__name__),
                    ]:
                        if isinstance(key, str) and key in score_kwargs:
                            kwargs = score_kwargs[key].copy()
                            break

                # Check if the score class accepts ref_seq parameter
                sig = inspect.signature(score_class.__init__)
                if "ref_seq" in sig.parameters:
                    # Pass the concatenated reference sequence
                    kwargs["ref_seq"] = concatenated_ref_seq

                # Instantiate the score
                try:
                    model = score_class(**kwargs)
                    models.append(model)
                    resolved_names.append(score_class.__name__)
                except TypeError as e:
                    raise TypeError(f"Failed to instantiate {score_class.__name__}: {e}") from e

            # Cache the built models with FIFO eviction
            # Evict oldest if cache would exceed max size
            while len(model_cache) >= max_cache_size:
                model_cache.popitem(last=False)  # Remove oldest (FIFO)
            model_cache[cache_key] = models

        # Use provided names or resolved names
        if names is None:
            # Use abbreviations where possible, otherwise class names
            names = [
                _get_score_abbreviation(name) if _get_score_abbreviation(name) else name
                for name in resolved_names
            ]

        # Create instance with models and cache
        instance = cls(models, names=names, max_cache_size=max_cache_size)
        instance._model_cache = model_cache
        return instance

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
        for name, model in zip(self.names, self.models, strict=True):
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
        counter = codonbias.stats.CodonCounter(seq_str)  # type: ignore[attr-defined]
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
            out_counter = codonbias.stats.CodonCounter(out_seq)  # type: ignore[attr-defined]
            out_counts = out_counter.get_codon_table()
            codon_counts -= out_counts

        # Add incoming codons
        if in_end > in_start:
            abs_in_start = orf_start + in_start
            abs_in_end = orf_start + in_end
            in_seq = get_seq_str(record)[abs_in_start:abs_in_end]
            in_seq = _convert_rna_to_dna(in_seq)
            in_counter = codonbias.stats.CodonCounter(in_seq)  # type: ignore[attr-defined]
            in_counts = in_counter.get_codon_table()
            codon_counts += in_counts

    def emit(self, state: dict) -> dict[str, float]:
        """Emit feature values from current state.

        Attempts incremental computation using codonbias.utils functions when weights
        are available. Falls back to baseline get_score() if incremental fails or
        weights are unavailable.

        Args:
            state: The state object

        Returns:
            Dictionary mapping feature names to scalar values
        """
        result = {}
        codon_counts = state["codon_counts"]
        model_data = state["model_data"]

        # Get window position for baseline fallback
        record = state["record"]
        orf = state["orf"]
        orf_start = orf[0]
        window_start = state["window_start"]
        window_end = state["window_end"]
        abs_start = orf_start + window_start
        abs_end = orf_start + window_end

        for name, model, weights_info in zip(self.names, self.models, model_data, strict=True):
            # Try incremental computation first if weights are available
            score = None
            if weights_info["has_weights"]:
                try:
                    weights = weights_info["weights"]
                    weight_type = weights_info["type"]

                    if weight_type == "log":
                        # Use geometric mean for log weights
                        score = codonbias.utils.geomean(  # type: ignore[attr-defined]
                            weights, codon_counts
                        )
                    else:
                        # Use arithmetic mean for regular weights
                        score = codonbias.utils.mean(  # type: ignore[attr-defined]
                            weights, codon_counts
                        )
                except Exception:
                    # Incremental computation failed - will use baseline fallback
                    score = None

            # Fall back to baseline if incremental failed or unavailable
            if score is None:
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
    model: Any,  # codonbias.scores.ScalarScore
) -> dict:
    """Try to extract weights from a codonbias model for incremental computation.

    This function attempts to get weights for optimized incremental computation,
    but incremental mode is purely optional. If weights cannot be obtained, the
    feature will automatically fall back to baseline per-window get_score() calls.

    Attempts to get weights in this order:
    1. model.log_weights (preferred) -> use with codonbias.utils.geomean
    2. model.weights -> use with codonbias.utils.mean
    3. model.get_weights() if callable -> use with codonbias.utils.mean

    Args:
        model: A codonbias ScalarScore model

    Returns:
        Dictionary with keys:
        - has_weights: bool indicating if weights were obtained
        - weights: pandas Series of weights (or None)
        - type: "log" or "linear" (or None)
    """
    # Try log_weights first (wrapped in try-except for safety)
    try:
        if hasattr(model, "log_weights"):
            log_weights = getattr(model, "log_weights", None)
            if log_weights is not None:
                return {
                    "has_weights": True,
                    "weights": log_weights,
                    "type": "log",
                }
    except Exception:
        pass

    # Try weights
    try:
        if hasattr(model, "weights"):
            weights = getattr(model, "weights", None)
            if weights is not None:
                return {
                    "has_weights": True,
                    "weights": weights,
                    "type": "linear",
                }
    except Exception:
        pass

    # Try get_weights() method
    try:
        if hasattr(model, "get_weights"):
            get_weights_method = getattr(model, "get_weights", None)
            if get_weights_method is not None and callable(get_weights_method):
                try:
                    weights = get_weights_method()
                    # Check if it returns a valid weights structure
                    if weights is not None:
                        return {
                            "has_weights": True,
                            "weights": weights,
                            "type": "linear",
                        }
                except (TypeError, AttributeError, ValueError):
                    pass
    except Exception:
        pass

    # No accessible weights - will use baseline fallback
    return {
        "has_weights": False,
        "weights": None,
        "type": None,
    }


def _resolve_score_identifier(score_id: str | type | Any) -> type | Any:
    """Resolve a score identifier to a score class or instance.

    Args:
        score_id: Score identifier (abbreviation, class name, class object, or instance)

    Returns:
        Score class or instance

    Raises:
        ValueError: If the identifier cannot be resolved
    """
    # If it's already a class or instance, return it
    if inspect.isclass(score_id) or not isinstance(score_id, str):
        return score_id

    # Define abbreviation mapping
    # Note: Some scores have alternative abbreviations commonly used in literature:
    # - DCBS (Distance from CBS) is an alternative name for RCBS
    # - CPS (Codon Pair Score) is an alternative name for CPB
    abbreviation_map = {
        "CAI": "CodonAdaptationIndex",
        "ENC": "EffectiveNumberOfCodons",
        "FOP": "FrequencyOfOptimalCodons",
        "RSCU": "RelativeSynonymousCodonUsage",
        "RCBS": "RelativeCodonBiasScore",
        "DCBS": "RelativeCodonBiasScore",  # Alternative name for RCBS
        "tAI": "TrnaAdaptationIndex",
        "nTE": "NormalizedTranslationalEfficiency",
        "CPB": "CodonPairBias",
        "CPS": "CodonPairBias",  # Alternative name for CPB
    }

    # Try to resolve as abbreviation
    class_name = abbreviation_map.get(score_id, score_id)

    # Try to get the class from codonbias.scores
    try:
        score_class = getattr(codonbias.scores, class_name)  # type: ignore[attr-defined]
        return score_class
    except AttributeError:
        raise ValueError(
            f"Cannot resolve score identifier '{score_id}'. "
            f"Expected one of: {', '.join(abbreviation_map.keys())}, "
            f"a codonbias.scores class name, or a class object."
        ) from None


def _get_score_abbreviation(class_name: str) -> str:
    """Get the abbreviation for a score class name.

    Args:
        class_name: Score class name

    Returns:
        Abbreviation if known, otherwise empty string
    """
    name_to_abbrev = {
        "CodonAdaptationIndex": "CAI",
        "EffectiveNumberOfCodons": "ENC",
        "FrequencyOfOptimalCodons": "FOP",
        "RelativeSynonymousCodonUsage": "RSCU",
        "RelativeCodonBiasScore": "RCBS",
        "TrnaAdaptationIndex": "tAI",
        "NormalizedTranslationalEfficiency": "nTE",
        "CodonPairBias": "CPB",
    }
    return name_to_abbrev.get(class_name, "")


def _generate_cache_key(
    reference_set: "ReferenceSequenceSet",
    scores: Sequence[str | type | Any],
    score_kwargs: dict[str, dict[str, Any]] | None,
) -> str:
    """Generate a deterministic cache key for model building.

    The cache key is based on:
    - reference_set.name (if present)
    - Content hash of cds_strings()
    - Score specifications
    - Score kwargs

    Args:
        reference_set: ReferenceSequenceSet containing CDS sequences
        scores: Sequence of score identifiers
        score_kwargs: Optional dict mapping score identifiers to their constructor kwargs

    Returns:
        Deterministic cache key string
    """
    # Start with reference set name if available
    key_parts = []
    if reference_set.name:
        key_parts.append(f"name:{reference_set.name}")

    # Hash the concatenated CDS strings for deterministic content-based key
    cds_strings = reference_set.cds_strings()
    concatenated_cds = "".join(cds_strings)
    cds_hash = hashlib.sha256(concatenated_cds.encode()).hexdigest()[:32]
    key_parts.append(f"cds:{cds_hash}")

    # Add score specifications (convert to strings for consistent hashing)
    score_strs = []
    for score_id in scores:
        if inspect.isclass(score_id):
            score_strs.append(score_id.__name__)
        elif isinstance(score_id, str):
            score_strs.append(score_id)
        else:
            # For instances, use class name + instance id to ensure uniqueness
            score_strs.append(f"{type(score_id).__name__}@{id(score_id)}")
    key_parts.append(f"scores:{','.join(score_strs)}")

    # Add score kwargs if provided
    if score_kwargs:
        # Use json.dumps for deterministic serialization
        kwargs_str = json.dumps(score_kwargs, sort_keys=True)
        kwargs_hash = hashlib.sha256(kwargs_str.encode()).hexdigest()[:32]
        key_parts.append(f"kwargs:{kwargs_hash}")

    return "|".join(key_parts)
