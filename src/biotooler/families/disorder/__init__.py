"""Disorder feature family for biotooler.

Protein disorder prediction using metapredict backend with optional future backends.
"""

__all__ = ["FAMILY_META", "get_features"]

# Family metadata for discovery and documentation
FAMILY_META = {
    "name": "disorder",
    "extra": "disorder",  # pip install biotooler[disorder]
    "owner": "@itamar-menuhin",  # GitHub username or team for CODEOWNERS
    "summary": "Protein disorder prediction using metapredict",
}


def get_features() -> list[type]:
    """Get list of feature classes provided by this family.

    This function lazily imports the feature module to keep the family
    import lightweight. The metapredict and idrpred dependencies are only
    loaded when the respective features are used.

    Returns:
        List of Feature classes (not instances) provided by this family.
        Users can instantiate these classes as needed.

    Raises:
        ImportError: If metapredict is not installed (when using
            DisorderProfileMetapredict) or if idrpred CLI is not available
            (when using IDRPredConsensusMask).
    """
    # Lazy import feature classes to avoid loading dependencies at family import
    from biotooler.families.disorder.derived import DisorderDerivedScalars
    from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask
    from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict

    return [DisorderProfileMetapredict, DisorderDerivedScalars, IDRPredConsensusMask]
