"""Type definitions for biotooler."""

# Scalar types for feature outputs
Scalar = int | float | str | bool

# Feature output is a dictionary mapping feature names to scalar values
FeatureOutput = dict[str, Scalar]
