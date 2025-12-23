# ViennaRNA Feature Family

## What this family provides

The ViennaRNA family provides RNA secondary structure prediction and analysis features using the ViennaRNA package. This family implements two types of windowing computations:

1. **Window MFE (minimum free energy)**: Substring-defined features that compute MFE for specific window substrings
2. **Accessibility (PU - Unpaired Probability)**: Vector-based features that compute per-nucleotide unpaired probabilities across the full sequence

**Implemented Features**:
- `ViennaRNAAccessibility`: Per-nucleotide unpaired probability vector (full-sequence computation)
- `WindowMFEFeature`: MFE for requested window substrings (substring-defined computation)
- `ContextWindowFoldFeature`: MFE and pairing metrics for windows with flanking context
- `WindowMFEStartRegion`: Convenience wrapper for MFE across a region
- `AccessibilityStartRegion`: Convenience wrapper for aggregated PU across a region

## Intuition

RNA secondary structure is the pattern of base pairing that forms when an RNA molecule folds upon itself. This structure is critical for RNA function, stability, and interactions. Understanding RNA folding helps in:

- **Predicting function**: Secondary structure determines which regions are accessible for binding
- **Designing therapeutics**: siRNA and antisense oligonucleotides require specific accessibility patterns
- **Understanding regulation**: Riboswitches and UTRs control gene expression through structural changes
- **Analyzing stability**: MFE indicates thermodynamic stability of RNA molecules

### Key Concepts

**Minimum Free Energy (MFE)**: The most thermodynamically stable secondary structure an RNA sequence can form. Lower (more negative) MFE values indicate more stable structures. MFE is computed for a specific sequence substring and represents a single scalar value.

**Unpaired Probability (PU / Accessibility)**: The probability that each nucleotide is unpaired (accessible) in the ensemble of all possible structures. Values range from 0.0 (always paired) to 1.0 (always unpaired). High accessibility indicates regions available for interactions. This is computed as a vector across the full sequence, with one probability per nucleotide.

**Base Pairing**: Watson-Crick pairs (G-C, A-U) and wobble pairs (G-U) that form the structure. ViennaRNA uses nearest-neighbor thermodynamic parameters to evaluate pairing stability.

## Mathematical formulation

### Minimum Free Energy (MFE)

The MFE structure minimizes the Gibbs free energy according to thermodynamic principles:

```
ΔG = ΔH - TΔS
```

Where:
- `ΔG` = Gibbs free energy change (kcal/mol) - negative values indicate stable structures
- `ΔH` = Enthalpy contribution from base stacking and hydrogen bonding
- `T` = Temperature in Kelvin (default: 310.15 K = 37°C)
- `ΔS` = Entropy contribution from conformational freedom

ViennaRNA computes MFE using the **nearest-neighbor model**, which sums contributions from adjacent base pairs and structural motifs (loops, bulges, hairpins). The algorithm uses dynamic programming (Zuker algorithm) to find the optimal structure.

**Upstream Algorithm**: ViennaRNA's `RNA.fold_compound(sequence).mfe()` returns both the dot-bracket structure notation and the MFE value. See [ViennaRNA folding documentation](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__mfe__fold.html).

### Partition Function and Unpaired Probability

The partition function Z sums Boltzmann-weighted contributions over all possible structures:

```
Z = Σ_i exp(-ΔG_i / RT)
```

Where:
- `ΔG_i` = Free energy of structure i
- `R` = Gas constant (0.001987 kcal/(mol·K))
- `T` = Temperature in Kelvin

From the partition function, ViennaRNA computes **base pairing probabilities** P(i,j) for each nucleotide pair (i, j). The **unpaired probability** (accessibility) for position i is:

```
PU(i) = 1 - Σ_j P(i,j)
```

This represents the probability that nucleotide i is NOT paired to any other position across the ensemble of all possible structures.

**Upstream Algorithm**: ViennaRNA's partition function is computed via `RNA.fold_compound(sequence).pf()`, and pairing probabilities are extracted using `plist_from_probs(cutoff)`. See [ViennaRNA partition function documentation](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__pf__fold.html).

## Features and output schema

### 1. ViennaRNAAccessibility

**Type**: PositionalFeature (vector-based)  
**Position Space**: RESIDUE (per-nucleotide)  
**Computation**: Full-sequence partition function → per-nucleotide PU vector

Computes unpaired probability (accessibility) for every nucleotide in the sequence. The feature performs full-sequence computation to maintain correct structural context.

**Output Schema**:
- `PU`: numpy array of float64, length = sequence length
  - Values in [0.0, 1.0] representing unpaired probability at each position
  - Aggregated using `np.mean` when windowing is applied

**Upstream API**: `RNA.fold_compound(seq).pf()` followed by `plist_from_probs(0.0)` to extract all base pairing probabilities.

**Use Case**: Identify accessible regions for primer binding, protein interactions, or regulatory element positioning.

### 2. WindowMFEFeature

**Type**: Wide feature (substring-defined)  
**Position Space**: RESIDUE (nucleotide positions define windows)  
**Computation**: Per-window substring folding (independent computations)

Computes MFE for specific requested window substrings. Each window is folded independently without context from surrounding sequence.

**Parameters**:
- `window_starts`: List of 0-based start positions for windows
- `window_size`: Size of each window in nucleotides

**Output Schema**:
- `MFE_<start>`: float (kcal/mol) for each requested window start
  - Example: `MFE_0`, `MFE_10`, `MFE_20`
  - Negative values indicate stable structures (more negative = more stable)

**Upstream API**: `RNA.fold_compound(window_substring).mfe()` called once per requested window.

**Use Case**: Compare stability of specific regions without flanking context influence.

### 3. ContextWindowFoldFeature

**Type**: Wide feature (substring-defined with context)  
**Position Space**: RESIDUE (nucleotide positions define windows)  
**Computation**: Per-window folding with flanking sequence context

Computes MFE for windows including flanking regions, then reports both context MFE and window-specific metrics. This captures how flanking sequence affects window structure.

**Parameters**:
- `starts_nt`: List of 0-based start positions for windows
- `window_size_nt`: Size of each window in nucleotides
- `flank_left_nt`: Number of nucleotides to include as left flank (default: 0)
- `flank_right_nt`: Number of nucleotides to include as right flank (default: 0)
- `mode`: Currently only "mfe" is supported

**Output Schema**:
- `CTX_MFE_<start>`: float (kcal/mol) - MFE of window + flanks
- `PAIR_OUT_FRAC_<start>`: float in [0.0, 1.0] - fraction of window nucleotides paired to positions outside the window

**Upstream API**: `RNA.fold_compound(context_substring).mfe()` where context_substring = left_flank + window + right_flank.

**Use Case**: Understand how surrounding sequence affects window stability and identify long-range interactions.

### 4. WindowMFEStartRegion

**Type**: Convenience wrapper around WindowMFEFeature  
**Computation**: Generates window starts for a region, delegates to WindowMFEFeature

**Parameters**:
- `region_start`: Start of region (0-based, inclusive)
- `region_end`: End of region (0-based, exclusive)
- `window_size`: Size of windows in nucleotides
- `step`: Step size between window starts

**Output Schema**: Same as WindowMFEFeature - `MFE_<start>` for each window in the region.

**Use Case**: Scan a specific region (e.g., 5' UTR, start codon region) for stability patterns.

### 5. AccessibilityStartRegion

**Type**: Convenience wrapper around ViennaRNAAccessibility  
**Computation**: Computes full PU vector once, then aggregates into windows

**Parameters**:
- `region_start`: Start of region (0-based, inclusive)
- `region_end`: End of region (0-based, exclusive)
- `window_size`: Size of windows in nucleotides
- `step`: Step size between window starts

**Output Schema**:
- `PU_<start>`: float in [0.0, 1.0] - mean unpaired probability for window starting at position <start>

**Use Case**: Quickly assess average accessibility across multiple windows in a region.

## References

- **ViennaRNA Package 2.0**: Lorenz, R., et al. (2011) "ViennaRNA Package 2.0" *Algorithms for Molecular Biology* 6:26. https://doi.org/10.1186/1748-7188-6-26
  - Primary reference for the ViennaRNA software suite
- **ViennaRNA Web Services**: Gruber, A.R., et al. (2008) "The Vienna RNA websuite" *Nucleic Acids Research* 36:W70-W74. https://doi.org/10.1093/nar/gkn188
  - Web interface and additional tools
- **Turner Energy Parameters**: Turner, D.H. and Mathews, D.H. (2010) "NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure" *Nucleic Acids Research* 38:D280-D282. https://doi.org/10.1093/nar/gkp892
  - Thermodynamic parameter database used by ViennaRNA
- **McCaskill Partition Function**: McCaskill, J.S. (1990) "The equilibrium partition function and base pair binding probabilities for RNA secondary structure" *Biopolymers* 29:1105-1119. https://doi.org/10.1002/bip.360290621
  - Algorithm for computing base pairing probabilities

## Upstream library links

- **ViennaRNA Package Homepage**: https://www.tbi.univie.ac.at/RNA/
- **ViennaRNA Documentation**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/index.html
- **Python API Documentation**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/api_python.html
- **MFE Folding Functions**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__mfe__fold.html
- **Partition Function API**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__pf__fold.html
- **GitHub Repository**: https://github.com/ViennaRNA/ViennaRNA
- **Installation Guide**: https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html

## Examples

### Example 1: ViennaRNAAccessibility (Vector-Based Feature)

Compute per-nucleotide unpaired probabilities for the entire sequence:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.viennarna import ViennaRNAAccessibility

# Create feature (no parameters needed)
accessibility = ViennaRNAAccessibility()

# Create RNA/DNA sequence record
record = SeqRecord(Seq("ACGUACGUACGUACGU"), id="example")

# Compute full PU vector (one value per nucleotide)
result = accessibility.compute_vector(record)
pu_vector = result["PU"]  # numpy array, length = 16

print(f"PU values: {pu_vector}")
print(f"Mean accessibility: {pu_vector.mean():.3f}")
print(f"Position 5 PU: {pu_vector[5]:.3f}")
```

**Key Points**:
- DNA sequences (with T) are automatically converted to RNA (with U)
- Full sequence is folded to maintain correct structural context
- Returns one PU value per nucleotide (vector-based computation)
- Values range from 0.0 (always paired) to 1.0 (always unpaired)

### Example 2: WindowMFEFeature (Substring-Defined Feature)

Compute MFE for specific window substrings (no flanking context):

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.viennarna import WindowMFEFeature

# Specify exactly which windows to compute
feature = WindowMFEFeature(
    window_starts=[0, 10, 20],  # Start positions (0-based)
    window_size=15              # Each window is 15 nt
)

# Create sequence record (must be at least 35 nt for window at position 20)
record = SeqRecord(Seq("ACGUACGUACGUACGUACGUACGUACGUACGUACGU"), id="example")

# Compute MFE for each requested window substring
result = feature(record)

print(f"MFE at position 0: {result['MFE_0']:.2f} kcal/mol")
print(f"MFE at position 10: {result['MFE_10']:.2f} kcal/mol")
print(f"MFE at position 20: {result['MFE_20']:.2f} kcal/mol")
# 'MFE_5' is NOT in result (not requested)
```

**Key Points**:
- Each window is folded independently as a substring (no surrounding context)
- Only requested window_starts are computed (efficient for sparse windows)
- Each window produces one scalar MFE value
- Windows extending beyond sequence are automatically skipped

### Example 3: ContextWindowFoldFeature (Substring with Context)

Compute MFE including flanking regions to capture long-range interactions:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.viennarna import ContextWindowFoldFeature

# Include flanking sequence in folding
feature = ContextWindowFoldFeature(
    starts_nt=[10],           # Window starts at position 10
    window_size_nt=20,        # Window is 20 nt
    flank_left_nt=5,          # Include 5 nt upstream
    flank_right_nt=5,         # Include 5 nt downstream
    mode="mfe"
)

record = SeqRecord(Seq("ACGU" * 15), id="example")  # 60 nt sequence

# Fold window [10:30] with flanks [5:35]
result = feature(record)

print(f"Context MFE: {result['CTX_MFE_10']:.2f} kcal/mol")
print(f"Fraction paired outside window: {result['PAIR_OUT_FRAC_10']:.2f}")
```

**Key Points**:
- Folds window + left_flank + right_flank as single sequence
- `CTX_MFE_<start>`: MFE of entire context slice
- `PAIR_OUT_FRAC_<start>`: Fraction of window nucleotides paired to flanking regions
- Captures how flanking sequence affects window structure

### Example 4: WindowMFEStartRegion (Convenience Wrapper)

Scan a region with automatically generated window starts:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.viennarna import WindowMFEStartRegion

# Scan region [0, 30) with 10 nt windows, step 5 nt
feature = WindowMFEStartRegion(
    region_start=0,
    region_end=30,
    window_size=10,
    step=5  # Windows at 0, 5, 10, 15, 20, 25
)

record = SeqRecord(Seq("ACGUACGU" * 5), id="example")  # 40 nt

result = feature(record)

for key in sorted(result.keys()):
    print(f"{key}: {result[key]:.2f} kcal/mol")
# Output: MFE_0, MFE_5, MFE_10, MFE_15, MFE_20, MFE_25
```

**Key Points**:
- Generates window_starts automatically: range(region_start, region_end, step)
- Delegates to WindowMFEFeature for computation
- Useful for scanning specific regions (e.g., 5' UTR, start codon region)

### Example 5: AccessibilityStartRegion (Convenience Wrapper)

Compute aggregated accessibility across windows in a region:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.viennarna import AccessibilityStartRegion

# Compute mean PU for windows in region [0, 30)
feature = AccessibilityStartRegion(
    region_start=0,
    region_end=30,
    window_size=10,
    step=5
)

record = SeqRecord(Seq("ACGUACGU" * 5), id="example")

result = feature(record)

for key in sorted(result.keys()):
    print(f"{key}: {result[key]:.3f}")
# Output: PU_0, PU_5, PU_10, PU_15, PU_20, PU_25
```

**Key Points**:
- Computes full PU vector once (efficient, no recomputation)
- Aggregates PU values into windows using mean
- Each window gets one scalar mean PU value

### Example 6: Using with Windowing Engine

Integrate with biotooler's windowing system:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.features import FeatureSet
from biotooler.families.viennarna import ViennaRNAAccessibility

# Create feature set with ViennaRNAAccessibility
feature = ViennaRNAAccessibility()
fs = FeatureSet(feature, name="viennarna")

record = SeqRecord(Seq("ACGUACGU" * 4), id="example")  # 32 nt

# Compute with sliding windows
result_df = fs.compute_orf_windows_v2(
    record,
    orf=(0, 32),        # Full sequence
    window_nt=12,       # 12 nt windows
    step_nt=6           # 6 nt step → windows at 0, 6, 12, 18
)

# Result is DataFrame with columns: viennarna.PU_0, viennarna.PU_6, viennarna.PU_12, viennarna.PU_18
print(result_df)
```

**Key Points**:
- FeatureSet handles window generation and aggregation
- ViennaRNAAccessibility returns vectors → automatically aggregated using mean
- WindowMFEFeature would require explicit window_starts (wide feature)

## Edge cases and validation

### Input Validation

**Sequence Requirements**:
- Accepts both DNA (A, C, G, T) and RNA (A, C, G, U) sequences
- DNA sequences are automatically normalized to RNA (T → U) before folding
- Uppercase and lowercase sequences are both supported (normalized to uppercase)
- Non-standard nucleotides (N, Y, R, etc.) are not supported by ViennaRNA's default energy model and will cause the underlying library to raise errors or produce undefined results

**Window Bounds**:
- Windows extending beyond sequence end are automatically skipped
- Empty window_starts lists raise `ValueError` during initialization
- Negative window sizes or steps raise `ValueError` during initialization

**Temperature and Energy Model**:
- Default temperature: 37°C (310.15 K) - matches typical physiological conditions
- Uses ViennaRNA's default Turner energy parameters (2004 version)
- No explicit temperature parameter exposed in current features

### Edge Cases

**Empty Sequences**:
- `ViennaRNAAccessibility.compute_vector()` on empty sequence returns `{"PU": np.array([], dtype=np.float64)}`
- WindowMFEFeature skips all windows if sequence is too short

**Short Sequences**:
- Sequences shorter than ~4 nucleotides may have limited structural options
- ViennaRNA still computes valid MFE and PU values, but structures are trivial

**Highly Structured vs. Unstructured**:
- Poly-A or poly-U sequences have high PU values (minimal pairing)
- GC-rich complementary sequences have low PU values (strong pairing)
- MFE becomes more negative with increasing structure stability

**Partial Windows**:
- WindowMFEFeature: Windows extending beyond sequence are skipped (not computed)
- AccessibilityStartRegion: Partial windows are clipped to sequence length and still computed
- When using windowing engine, `drop_partial=True` (default) excludes partial windows

**Numerical Precision**:
- PU values should be in [0.0, 1.0], but floating-point errors may cause slight deviations (< 1e-6)
- MFE values are in kcal/mol with precision to ~0.01 kcal/mol

### Common Pitfalls

**Context Matters**: WindowMFEFeature folds substrings independently, which may miss long-range interactions. Use ContextWindowFoldFeature when flanking context is important.

**Full-Sequence Computation**: ViennaRNAAccessibility always computes the full PU vector to maintain structural correctness, even if only specific positions are needed. This is intentional for accuracy.

**DNA vs RNA**: Input sequences with T are accepted but converted to U. Results are identical for equivalent DNA/RNA sequences (e.g., "ATG" and "AUG" give same results).

**Window Start Semantics**: 
- WindowMFEFeature with `window_starts=[0, 10]` and `window_size=15` computes:
  - Window 1: sequence[0:15]
  - Window 2: sequence[10:25]
- Overlapping windows are allowed and computed independently

## Windowing correctness

This section describes how ViennaRNA features handle windowing and ensures computational correctness.

### Position space

All ViennaRNA features operate at the **RESIDUE position space** (per-nucleotide level). Position indices are 0-based and refer to nucleotide positions in the sequence.

**Example**:
- Sequence: `ACGUACGU` (8 nucleotides)
- Position 0: A, Position 1: C, Position 2: G, etc.
- A window [0:4] contains nucleotides ACGU (positions 0, 1, 2, 3)

### Vector computation

ViennaRNA features use **two distinct windowing paradigms** depending on the feature type:

**1. Vector-Based Features (ViennaRNAAccessibility)**

Computation Strategy: Full-sequence, vector-based

The `ViennaRNAAccessibility` feature computes unpaired probabilities as a **vector across the entire sequence**:

1. **Full-sequence partition function**: Uses `RNA.fold_compound(full_seq).pf()` to compute base pairing probabilities considering all possible long-range interactions
2. **Per-nucleotide PU values**: Derives unpaired probability for each position: `PU(i) = 1 - Σ_j P(i,j)`
3. **Vector output**: Returns numpy array with length equal to sequence length
4. **Window aggregation**: When windowing is applied (e.g., via FeatureSet), PU vectors are aggregated using `np.mean` over each window

**Why full-sequence computation?** RNA secondary structure involves long-range base pairing (e.g., a nucleotide at position 10 can pair with position 50). Computing PU on a substring would ignore these long-range interactions and produce incorrect accessibility values.

**Upstream API calls**:
```python
# biotooler implementation (simplified)
fc = RNA.fold_compound(full_sequence)
fc.pf()  # Compute partition function
pairs = fc.plist_from_probs(0.0)  # Get all base pairing probabilities
# Derive PU from pairs as described in accessibility.py
```

See ViennaRNA [partition function documentation](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__pf__fold.html).

**2. Substring-Based Features (WindowMFEFeature, ContextWindowFoldFeature)**

Computation Strategy: Per-window substring folding

Window MFE features compute MFE for **specific substrings independently** using upstream ViennaRNA library APIs:

1. **Substring extraction**: Extract `seq[window_start:window_start+window_size]` for each requested window
2. **Independent folding**: Call `RNA.fold_compound(substring).mfe()` on the substring alone
3. **Scalar output**: Each window produces one MFE value (kcal/mol)
4. **No aggregation**: MFE is already a scalar per window (not a vector)

**Why substring computation?** MFE represents the single most stable structure for a given sequence. For localized stability analysis, computing MFE on a substring shows the intrinsic folding potential of that region without being dominated by stronger long-range interactions elsewhere in the sequence.

**Context variant**: `ContextWindowFoldFeature` includes flanking sequence in the fold:
```python
# Window at position 10, size 20, flanks 5 nt each
context = seq[5:35]  # 5 upstream + 20 window + 5 downstream
fc = RNA.fold_compound(context)
structure, mfe = fc.mfe()
# Then analyze which window nucleotides pair with flanks
```

**Upstream API calls**:
```python
# biotooler implementation (simplified)
window_seq = full_sequence[start:start+window_size]
rna_seq = window_seq.replace("T", "U")  # DNA to RNA
fc = RNA.fold_compound(rna_seq)
structure, mfe = fc.mfe()  # Returns (dot-bracket, energy)
```

See ViennaRNA [MFE folding documentation](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__mfe__fold.html).

**Key Difference**:
- **Accessibility (PU)**: Vector-based → full sequence computation → window aggregation
- **MFE**: Substring-based → per-window computation → no aggregation needed

### Aggregation strategy

#### Vector-Based Features (ViennaRNAAccessibility)

**Aggregation Function**: `np.mean` (arithmetic mean)

When PU vectors are aggregated into windows:
1. Compute full PU vector: `[PU_0, PU_1, ..., PU_n]`
2. For window [start:end], extract: `PU_window = PU_vector[start:end]`
3. Aggregate: `mean_PU = np.mean(PU_window)`

**Why mean?** Average unpaired probability represents the overall accessibility of a window region. A window with mean PU = 0.8 indicates generally high accessibility.

**Example**:
```python
# PU vector: [0.9, 0.8, 0.7, 0.6, 0.5, 0.4]
# Window [1:4] contains PU values: [0.8, 0.7, 0.6]
# Aggregated value: mean([0.8, 0.7, 0.6]) = 0.7
```

#### Substring-Based Features (WindowMFE*)

**Aggregation**: None (already scalar per window)

MFE features produce one scalar value per window. No aggregation is needed or performed. Each window's MFE is independent.

**Output**:
- WindowMFEFeature with `window_starts=[0, 10, 20]` produces three scalar values: `MFE_0`, `MFE_10`, `MFE_20`

### Testing approach

Windowing correctness is validated through comprehensive tests:

#### 1. Upstream Library Correctness

**Vector features** (`test_accessibility.py`):
- Compare full-vector results to manual aggregation
- Verify window engine v2 aggregations match `np.mean` of vector slices
- Test with various window sizes and steps

**Substring features** (`test_window_mfe.py`):
- Compare each window's MFE to direct ViennaRNA call on same substring
- Verify: `feature(record)["MFE_X"]` == `RNA.fold_compound(seq[X:X+size]).mfe()[1]`
- Test ensures no hidden transformations or context bleed

#### 2. Position Space Correctness

- Verify vector lengths match sequence lengths
- Confirm window start/end positions map correctly to nucleotides
- Test edge cases: sequence boundaries, partial windows, empty sequences

#### 3. Full-Context Computation

**ViennaRNAAccessibility**:
- Test that `positions` argument doesn't affect computation (full vector always computed)
- Verify long-range pairing is captured (compare with known hairpin structures)

**WindowMFEFeature**:
- Test that each window is folded independently (no context bleed between windows)
- Verify windows extending beyond sequence are skipped

#### 4. Edge Case Handling

- Empty sequences return empty vectors or skip all windows
- Short sequences (< window size) handled gracefully
- DNA/RNA normalization tested (T → U conversion)
- Determinism verified (same input → same output)

#### 5. Windowing Engine Integration

Tests verify that biotooler's windowing engine correctly:
- Calls `compute_vector()` once per sequence (not per window)
- Aggregates vector slices using specified aggregation function
- Handles window boundaries and partial windows correctly
- Produces wide-format output with correct column naming

## Maintenance notes

### Dependencies

**Required**:
- Python ≥ 3.9
- `biopython` - For SeqRecord and Seq objects
- `numpy` - For vector computations and aggregation

**Optional**:
- `ViennaRNA` Python package (required to use this family)
  - Install via: `pip install ViennaRNA` or `pip install "biotooler[viennarna]"`
  - ViennaRNA is a large C/C++ library with Python bindings
  - Not included in default biotooler installation to keep dependencies lightweight

### Installation

**Install biotooler with ViennaRNA support**:
```bash
pip install "biotooler[viennarna]"
```

**Install ViennaRNA separately**:
```bash
pip install ViennaRNA
```

**Check installation**:
```python
import RNA  # Should not raise ImportError
print(RNA.__version__)
```

If ViennaRNA is not installed, importing features will succeed but calling them will raise `ImportError` with installation instructions.

### Implementation status

- ✅ Family infrastructure (registry, imports, lazy loading)
- ✅ `ViennaRNAAccessibility` - Per-nucleotide unpaired probability (vector-based)
- ✅ `WindowMFEFeature` - Window substring MFE (substring-based)
- ✅ `ContextWindowFoldFeature` - Window MFE with flanking context
- ✅ `WindowMFEStartRegion` - Convenience wrapper for region scanning (MFE)
- ✅ `AccessibilityStartRegion` - Convenience wrapper for region scanning (PU)
- ✅ Comprehensive test coverage with ViennaRNA integration tests
- ✅ Windowing engine v2 integration

### Design decisions

**1. Lazy imports**: ViennaRNA is only loaded when features are accessed, not when the module is imported
   - Rationale: ViennaRNA is a heavy C/C++ library; lazy loading keeps biotooler import fast
   - Implementation: `require_viennarna()` helper checks and imports RNA module on first use

**2. Full-sequence computation for accessibility**: Always compute PU for entire sequence
   - Rationale: RNA structure is global; substring computation would miss long-range interactions
   - Trade-off: Slightly slower, but guarantees correct structural context

**3. Substring computation for MFE**: Fold each window independently without flanks
   - Rationale: Localized stability analysis; compare intrinsic folding potential
   - Alternative: Use `ContextWindowFoldFeature` when flanking context matters

**4. DNA to RNA normalization**: Automatically convert T → U
   - Rationale: ViennaRNA expects RNA (ACGU), but users often work with DNA sequences
   - Implementation: `.replace("T", "U")` before folding

**5. Wide feature format for MFE**: Return `{"MFE_0": -5.2, "MFE_10": -3.8, ...}`
   - Rationale: Each window is independent; wide format avoids recomputation
   - Contrast: Vector features return arrays that are aggregated by windowing engine

**6. Mean aggregation for PU**: Use `np.mean` instead of max or min
   - Rationale: Average accessibility represents overall window behavior
   - Alternative: Users can access full PU vector and apply custom aggregation

### Performance considerations

**ViennaRNAAccessibility**:
- Time complexity: O(n³) for sequence of length n (partition function algorithm)
- Space complexity: O(n²) for base pairing probability matrix
- Computed once per sequence, regardless of number of windows
- Typical timing: ~10ms for 100 nt, ~100ms for 500 nt, ~1s for 1000 nt

**WindowMFEFeature**:
- Time complexity: O(k × w³) for k windows of size w
- Space complexity: O(w²) per window
- Each window folded independently (can be parallelized in principle)
- Typical timing: ~5ms per 50 nt window

**ContextWindowFoldFeature**:
- Time complexity: O(k × c³) for k windows with context size c = w + flanks
- Larger context = slower computation but better structural accuracy

### Future work

**Potential enhancements**:
- Support for constraints (force/prohibit specific base pairs)
- Temperature parameter control (currently uses ViennaRNA default 37°C)
- Suboptimal structure ensemble analysis
- Base pairing probability matrix export
- RNA-RNA interaction prediction features
- Cofolding for multiple sequences

**Optimization opportunities**:
- Cache partition function results for repeated queries on same sequence
- Parallel computation for multiple windows (thread pool)
- Sparse PU computation for large sequences (though this sacrifices accuracy)

### Troubleshooting

**Import Error: "No module named 'RNA'"**
- Solution: Install ViennaRNA with `pip install ViennaRNA` or `pip install "biotooler[viennarna]"`

**Unexpected MFE values (very positive or NaN)**
- Check: Sequence contains only valid nucleotides (ACGU or ACGT)
- Check: Window size is reasonable (at least 4-5 nucleotides for meaningful structure)

**PU values outside [0, 1]**
- Rare floating-point errors may cause slight deviations (< 1e-6)
- If values are far outside range, check for ViennaRNA version compatibility

**Slow computation on long sequences**
- ViennaRNAAccessibility: O(n³) complexity; consider shorter sequences or targeted regions
- Use `AccessibilityStartRegion` to compute only specific regions efficiently

### Version compatibility

- Tested with ViennaRNA 2.5.x and 2.6.x
- Python bindings API is stable across recent ViennaRNA versions
- biotooler >= 0.1.0 required for windowing engine v2 support
