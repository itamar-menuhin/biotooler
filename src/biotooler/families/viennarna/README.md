# ViennaRNA Feature Family

## What this family provides

The ViennaRNA family will provide RNA secondary structure prediction and analysis features using the ViennaRNA package. Features are planned for future releases.

**Status**: Infrastructure only - no features implemented yet.

## Intuition

RNA secondary structure plays a crucial role in RNA function, stability, and interactions. The ViennaRNA package provides state-of-the-art algorithms for:
- Minimum free energy (MFE) structure prediction
- Partition function calculations
- Base pairing probabilities
- Suboptimal structure ensembles
- RNA-RNA interaction prediction

These analyses are fundamental for:
- Understanding RNA folding and stability
- Designing synthetic RNA constructs (siRNA, aptamers, riboswitches)
- Predicting regulatory RNA structures (riboswitches, UTRs)
- Analyzing evolutionary conservation of RNA structures

## Mathematical formulation

Features in this family will use thermodynamic models and algorithms from ViennaRNA:

### Minimum Free Energy (MFE) Structure
The MFE structure minimizes the Gibbs free energy:
```
ΔG = ΔH - TΔS
```
Where:
- `ΔG` = Gibbs free energy change
- `ΔH` = Enthalpy (stacking, hydrogen bonding)
- `T` = Temperature (typically 37°C)
- `ΔS` = Entropy (conformational freedom)

ViennaRNA uses nearest-neighbor thermodynamic parameters to calculate folding energy.

### Partition Function
The partition function Z sums over all possible structures:
```
Z = Σ exp(-ΔG_i / RT)
```
This enables calculation of base pairing probabilities and ensemble properties.

**Note**: Specific formulations will be documented when features are implemented.

## Features and output schema

**Status**: No features implemented yet.

Future features may include:
- MFE structure and energy
- Base pairing probabilities
- Ensemble diversity metrics
- Structural accessibility
- RNA-RNA interaction predictions

Output schemas will be documented when features are added.

## References

- **ViennaRNA Package**: Lorenz, R., et al. (2011) "ViennaRNA Package 2.0" *Algorithms for Molecular Biology* 6:26. https://doi.org/10.1186/1748-7188-6-26
- **ViennaRNA Web Services**: Gruber, A.R., et al. (2008) "The Vienna RNA websuite" *Nucleic Acids Research* 36:W70-W74. https://doi.org/10.1093/nar/gkn188
- **Turner Energy Parameters**: Turner, D.H. and Mathews, D.H. (2010) "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure" *Nucleic Acids Research* 38:D280-D282. https://doi.org/10.1093/nar/gkp892

## Upstream library links

- **ViennaRNA Package**: https://www.tbi.univie.ac.at/RNA/
- **ViennaRNA Documentation**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/index.html
- **Python Bindings**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/api_python.html
- **GitHub Repository**: https://github.com/ViennaRNA/ViennaRNA

## Examples

**Status**: No features implemented yet.

Example usage will be provided when features are added. The general pattern will be:

```python
from biotooler.families.viennarna import get_features

# This will work once features are implemented
# features = get_features()
# feature = features[0]()
# result = feature(rna_sequence)
```

## Edge cases and validation

**Status**: To be documented when features are implemented.

Expected considerations:
- RNA sequence validation (A, C, G, U only)
- Handling of ambiguous nucleotides
- Minimum sequence length requirements
- Temperature parameter validation
- Structure constraint handling

## Windowing correctness

### Position space

This family will operate at the **RNA nucleotide level** (RESIDUE position space). Features will be computed on RNA sequences with position-specific information where applicable.

### Vector computation

**Status**: To be determined based on specific features implemented.

When features are added, computation will use ViennaRNA's Python API. The approach will depend on the specific feature:
- **Global structure features** (e.g., MFE): Compute on full sequence via upstream library
- **Local features** (e.g., accessibility): May support positional computation via upstream library APIs

All computations will properly interface with ViennaRNA's C library through Python bindings to ensure correctness.

### Aggregation strategy

**Status**: To be determined based on specific features implemented.

Aggregation strategies will be chosen based on the biological meaning of each feature:
- **Structure-level features**: Likely no aggregation (single value per structure)
- **Position-level features**: May use appropriate aggregation (mean, max, etc.)

Documentation will be updated when features are implemented.

### Testing approach

Windowing correctness tests will be added alongside feature implementations:
1. **Upstream library correctness**: Validate results match direct ViennaRNA API calls
2. **Position space correctness**: Verify proper coordinate mapping
3. **Full-context computation**: Ensure features use appropriate sequence context
4. **Edge case handling**: Test boundary conditions and special cases

## Maintenance notes

### Dependencies
- **ViennaRNA Python bindings** (optional): Install via `pip install "biotooler[viennarna]"` or `pip install ViennaRNA`
- Package uses lazy imports to avoid loading ViennaRNA at module import time

### Implementation status
- ✅ Family infrastructure (registry, imports, tests)
- ✅ Lazy import integration with helpful error messages
- ⏳ Feature implementations (planned for future releases)

### Installation
```bash
pip install "biotooler[viennarna]"
```

### Design decisions
1. **Lazy imports**: ViennaRNA is a heavy C/C++ library, so it's only loaded when features are accessed
2. **Infrastructure first**: Family registered and integrated before feature implementation
3. **Marked as heavy**: Registry indicates this family has significant dependencies

### Future work
- Implement core structure prediction features
- Add base pairing probability features
- Support RNA-RNA interaction predictions
- Add ensemble diversity metrics
