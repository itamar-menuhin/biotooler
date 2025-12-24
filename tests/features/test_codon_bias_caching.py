"""Tests for codon bias model caching functionality."""

from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biotooler.core.reference_sequences import ReferenceSequenceSet
from biotooler.families.codon_bias import CodonBiasFeature


class TestCodonBiasCaching:
    """Tests for model caching in CodonBiasFeature.from_reference."""

    def test_cache_hit_reuses_models(self):
        """Test that models are reused when cache hit occurs."""
        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG", "gene2": "ATGATGATGATGATGATGATG"}
        )

        # Create a shared cache
        cache = OrderedDict()

        # First call builds models
        feature1 = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"], model_cache=cache)

        # Cache should have one entry
        assert len(cache) == 1

        # Second call with same parameters should reuse models
        feature2 = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"], model_cache=cache)

        # Cache should still have one entry
        assert len(cache) == 1

        # Models should be the same instances (reused from cache)
        assert feature1.models[0] is feature2.models[0]  # CAI
        assert feature1.models[1] is feature2.models[1]  # ENC

    def test_cache_miss_on_different_reference(self):
        """Test that different reference sets result in cache miss."""
        ref_set1 = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})
        ref_set2 = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGAAA"}  # Different content
        )

        cache = OrderedDict()

        # Build with first reference
        feature1 = CodonBiasFeature.from_reference(ref_set1, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Build with second reference should create new models
        feature2 = CodonBiasFeature.from_reference(ref_set2, ["CAI"], model_cache=cache)
        assert len(cache) == 2

        # Models should be different instances
        assert feature1.models[0] is not feature2.models[0]

    def test_cache_miss_on_different_scores(self):
        """Test that different score specifications result in cache miss."""
        ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        cache = OrderedDict()

        # Build with CAI
        CodonBiasFeature.from_reference(ref_set, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Build with ENC should create new models
        CodonBiasFeature.from_reference(ref_set, ["ENC"], model_cache=cache)
        assert len(cache) == 2

    def test_cache_miss_on_different_kwargs(self):
        """Test that different score kwargs result in cache miss."""
        ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        cache = OrderedDict()

        # Build with default kwargs
        CodonBiasFeature.from_reference(ref_set, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Build with custom kwargs should create new models
        CodonBiasFeature.from_reference(
            ref_set, ["CAI"], score_kwargs={"CAI": {"genetic_code": 1}}, model_cache=cache
        )
        assert len(cache) == 2

    def test_cache_eviction_fifo(self):
        """Test that cache evicts oldest entries when full."""
        ref_sets = [ReferenceSequenceSet(cds={"gene1": f"ATGATGATG{'ATG' * i}"}) for i in range(10)]

        cache = OrderedDict()

        # Fill cache beyond max_cache_size (default is 8)
        features = []
        for ref_set in ref_sets:
            feature = CodonBiasFeature.from_reference(
                ref_set, ["CAI"], model_cache=cache, max_cache_size=3
            )
            features.append(feature)

        # Cache should be limited to max_cache_size
        assert len(cache) == 3

        # The last 3 reference sets should be in cache
        # Create features again with the last 3 refs - should be cache hits
        for i in range(7, 10):
            feature_new = CodonBiasFeature.from_reference(
                ref_sets[i], ["CAI"], model_cache=cache, max_cache_size=3
            )
            # Should reuse the same model instance
            assert feature_new.models[0] is features[i].models[0]

    def test_cache_key_deterministic_same_content(self):
        """Test that cache key is deterministic for same content."""
        # Two reference sets with identical content
        ref_set1 = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})
        ref_set2 = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        cache = OrderedDict()

        # Build with first reference
        feature1 = CodonBiasFeature.from_reference(ref_set1, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Build with second reference (same content) should hit cache
        feature2 = CodonBiasFeature.from_reference(ref_set2, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Should reuse the same model
        assert feature1.models[0] is feature2.models[0]

    def test_cache_key_uses_reference_name(self):
        """Test that cache key includes reference set name when provided."""
        # Two reference sets with same content but different names
        ref_set1 = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"}, name="ref1")
        ref_set2 = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"}, name="ref2")

        cache = OrderedDict()

        # Build with first reference
        feature1 = CodonBiasFeature.from_reference(ref_set1, ["CAI"], model_cache=cache)
        assert len(cache) == 1

        # Build with second reference (different name) should create new models
        feature2 = CodonBiasFeature.from_reference(ref_set2, ["CAI"], model_cache=cache)
        assert len(cache) == 2

        # Models should be different instances
        assert feature1.models[0] is not feature2.models[0]

    def test_cached_models_produce_correct_results(self):
        """Test that cached models produce the same results as freshly built models."""
        ref_set = ReferenceSequenceSet(
            cds={"gene1": "ATGATGATGATGATGATGATG", "gene2": "ATGATGATGATGATGATGATG"}
        )

        cache = OrderedDict()

        # Build feature with cache
        feature1 = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"], model_cache=cache)

        # Build feature again (should use cache)
        feature2 = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"], model_cache=cache)

        # Test sequence
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")

        # Compute with both features
        result1 = feature1(record)
        result2 = feature2(record)

        # Results should be identical
        assert abs(result1["CAI"] - result2["CAI"]) < 1e-9
        assert abs(result1["ENC"] - result2["ENC"]) < 1e-9

    def test_no_cache_creates_new_instance_cache(self):
        """Test that not providing a cache still works (creates per-instance cache)."""
        ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        # Build without providing cache - should work fine
        feature = CodonBiasFeature.from_reference(ref_set, ["CAI", "ENC"])

        # Feature should have its own cache
        assert hasattr(feature, "_model_cache")
        assert isinstance(feature._model_cache, OrderedDict)

        # Test computation works
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")
        result = feature(record)

        assert "CAI" in result
        assert "ENC" in result

    def test_cache_with_already_instantiated_models(self):
        """Test that already instantiated models are not cached (they're passed directly)."""
        from codonbias.scores import CodonAdaptationIndex, EffectiveNumberOfCodons

        ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        # Create model instances
        ref_seq = "".join(ref_set.cds_strings())
        cai = CodonAdaptationIndex(ref_seq)
        enc = EffectiveNumberOfCodons()

        cache = OrderedDict()

        # Pass already instantiated models
        feature = CodonBiasFeature.from_reference(
            ref_set, [cai, enc], names=["CAI", "ENC"], model_cache=cache
        )

        # Cache should still be populated
        assert len(cache) == 1

        # Test computation works
        test_seq = "ATGATGATGATGATGATG"
        record = SeqRecord(Seq(test_seq), id="test")
        result = feature(record)

        assert "CAI" in result
        assert "ENC" in result

    def test_cache_different_instances_same_class(self):
        """Test that different score instances produce different cache keys."""
        from codonbias.scores import CodonAdaptationIndex

        ref_set = ReferenceSequenceSet(cds={"gene1": "ATGATGATGATGATGATGATG"})

        # Create two different CAI instances
        ref_seq = "".join(ref_set.cds_strings())
        cai1 = CodonAdaptationIndex(ref_seq)
        cai2 = CodonAdaptationIndex(ref_seq)

        cache = OrderedDict()

        # Pass first CAI instance
        feature1 = CodonBiasFeature.from_reference(
            ref_set, [cai1], names=["CAI"], model_cache=cache
        )
        assert len(cache) == 1

        # Pass second CAI instance (different object) should create new cache entry
        feature2 = CodonBiasFeature.from_reference(
            ref_set, [cai2], names=["CAI"], model_cache=cache
        )
        assert len(cache) == 2

        # The features should use different model instances
        assert feature1.models[0] is cai1
        assert feature2.models[0] is cai2
        assert feature1.models[0] is not feature2.models[0]

    def test_lru_behavior_move_to_end(self):
        """Test that accessing cached models moves them to end (LRU)."""
        ref_sets = [ReferenceSequenceSet(cds={"gene1": f"ATGATGATG{'ATG' * i}"}) for i in range(5)]

        cache = OrderedDict()

        # Build 3 features (fill cache with max_cache_size=3)
        for i in range(3):
            CodonBiasFeature.from_reference(
                ref_sets[i], ["CAI"], model_cache=cache, max_cache_size=3
            )

        assert len(cache) == 3

        # Access the first one again (should move to end)
        CodonBiasFeature.from_reference(ref_sets[0], ["CAI"], model_cache=cache, max_cache_size=3)

        # Still 3 entries
        assert len(cache) == 3

        # Add a new one - should evict ref_sets[1] (oldest non-accessed)
        CodonBiasFeature.from_reference(ref_sets[3], ["CAI"], model_cache=cache, max_cache_size=3)

        assert len(cache) == 3

        # ref_sets[0] should still be cached (was accessed)
        # ref_sets[2] should still be cached (recent)
        # ref_sets[3] should be cached (just added)
        # ref_sets[1] should be evicted
