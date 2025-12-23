# Chimera Feature Family

## What this family provides

The chimera family provides feature computation for gene expression prediction and sequence optimization using the Chimera algorithms from the pyChimera package. These algorithms can predict the expression level of a gene in an unsupervised manner based solely on the coding sequence and the host genome. The algorithms measure sequence adaptation by comparing target genes against reference gene sets to quantify similarity patterns using maximal common substring analysis.

## Intuition

Gene expression levels are influenced not only by regulatory elements but also by the coding sequence itself. The Chimera algorithms exploit "hidden information" embedded in the redundancy of the genetic code—information beyond what is encoded in the amino acid sequence. This information can reflect evolutionary optimization for translation efficiency, mRNA stability, codon usage patterns, and other expression-related factors.

The core insight is that highly expressed genes in a host organism tend to share common sequence patterns. By measuring how well a target gene's sequence matches patterns found in a reference set of host genes (e.g., highly expressed native genes), we can:

1. **Predict expression levels**: Genes with higher similarity to reference gene patterns are predicted to express better in the host
2. **Optimize sequences**: Design synonymous codon variants that maximize similarity to reference patterns while maintaining the same protein sequence
3. **Detect adaptation**: Identify whether genes are well-adapted to their host organism's translational machinery

The Chimera algorithms work by finding maximal common substrings between target and reference sequences at every position in the target sequence, either globally (cARS/cMap) or with position-specific constraints that respect gene structure (PScARS/PScMap). The key mathematical operation is computing the longest substring starting at each position that appears in the reference set—these lengths are then averaged to produce an adaptation score.

## Mathematical formulation

### Alphabets and Units

The pyChimera library supports computation in three different alphabets, each with different granularity:

1. **Nucleotide (nt)**: Individual bases A, C, G, T (or U for RNA)
   - Window sizes and positions are measured in nucleotides
   - Most fine-grained representation

2. **Codon**: Triplets of nucleotides (e.g., ATG, GCC, TAG)
   - Window sizes and positions are measured in codons (1 codon = 3 nucleotides)
   - Standard representation for coding sequences
   - Used by default in `ChimeraFeature` implementation

3. **Amino Acid (aa)**: Single amino acids (e.g., M, A, L)
   - Used for protein sequence optimization (cMap/PScMap)
   - Not currently supported for cARS scoring in this implementation

**In this biotooler implementation**, the `ChimeraFeature` class operates at the **codon level**. This means:
- Input sequences are DNA/RNA (not protein) that are converted to codon representation via `nt2codon()`
- Window sizes in biotooler windowing APIs are specified in **codon units**
- Position space is `PositionSpace.CODON`
- Each "position" in the cARS vector corresponds to one codon (3 nucleotides)

### ChimeraARS (cARS) - Average Repetitive Substring

The cARS score measures the average length of maximal common substrings between a target gene and a reference set at every position:

```
cARS(target) = (1/L) * Σ(i=1 to L) max_length(target[i:], reference_set)
```

Where:
- L is the length of the target sequence
- max_length finds the longest substring starting at position i that appears in the reference set
- Higher scores indicate better adaptation to the reference set

**Position-Specific cARS (PScARS)** extends this by constraining matches to similar positions within genes (e.g., matching beginning regions with beginning regions), capturing position-dependent regulatory signals.

### ChimeraMap (cMap) - Sequence Optimization

Given a target amino acid sequence, cMap generates an optimized nucleotide sequence by selecting codons that maximize coverage by substrings found in the reference set:

```
optimized_seq = argmax Σ coverage(substring) * weight(substring)
```

The algorithm constructs the sequence from minimal sequence blocks that appear in reference genes, using a greedy approach with suffix arrays for efficient substring matching.

**Position-Specific cMap (PScMap)** adds position constraints, and **Multi-sequence cMap (MScMap)** generates multiple diverse optimized variants for multi-copy systems.

### pyChimera API Functions

This family uses the following core functions from the pyChimera library:

1. **`nt2codon(sequences: list[str]) -> list[list[str]]`**
   - Converts DNA/RNA nucleotide sequences to codon representation
   - Input: List of nucleotide strings (DNA with T or RNA with U)
   - Output: List of codon lists, where each codon is a 3-character string
   - Handles sequence normalization (RNA→DNA conversion with U→T)

2. **`build_suffix_array(codon_sequences, pos_spec=False)`**
   - Builds a suffix array data structure from reference sequences
   - Input: Codon sequences (output of `nt2codon()`)
   - Parameter `pos_spec`: If True, includes position-specific information for PScARS
   - Output: Suffix array object used for efficient substring matching
   - This is a preprocessing step that can be cached per reference set

3. **`calc_cARS(target_codon, suffix_array, win_params=None, max_len=40, max_pos=0.5, return_vec=False)`**
   - Computes cARS (or PScARS if `win_params` provided) score or vector
   - Parameters:
     - `target_codon`: Single codon sequence (one element from `nt2codon()` output)
     - `suffix_array`: Prebuilt suffix array from reference sequences
     - `win_params`: Dict with window parameters for PScARS (keys: 'size', 'center', 'by_start', 'by_stop')
     - `max_len`: Maximum substring length to consider (for homolog filtering, in codon units)
     - `max_pos`: Maximum position difference for position-specific matching (fraction of sequence length)
     - `return_vec`: If False, returns scalar score (mean). If True, returns per-position vector
   - Output: 
     - If `return_vec=False`: Single float (average maximal substring length across all positions)
     - If `return_vec=True`: List/array of floats (maximal substring length at each codon position)

The key design: **`return_vec=True` enables correct windowing** by computing per-position values on the full sequence, which are then aggregated into windows by biotooler's windowing machinery.

## Features and output schema

### Input

- **Reference sequences**: Set of host genes (DNA/RNA sequences) used as the adaptation template
  - Typically highly expressed genes from the host organism
  - Converted to codon representation via `nt2codon()` during initialization
  - Used to build suffix array (cached after first use)
  
- **Target sequences**: Query genes to analyze (DNA/RNA for cARS, amino acid for cMap)
  - For cARS/PScARS: DNA or RNA sequences (automatically normalized to DNA)
  - Converted to codon representation via `nt2codon()` before scoring
  
- **Parameters**:
  - `algorithm`: Which algorithm to use ("cARS", "PScARS", "cMap", "PScMap")
    - Currently only "cARS" and "PScARS" are implemented
  - `max_len`: Maximum substring length to consider (default: 40 codons)
    - Used for homolog filtering to avoid spurious matches
  - `max_pos`: Maximum position difference for position-specific algorithms (default: 0.5)
    - Fraction of sequence length; controls how far from aligned positions matches are allowed
  - `win_params`: Window parameters for position-specific algorithms (dict with keys: 'size', 'center', 'by_start', 'by_stop')
    - Only used for PScARS/PScMap variants
    - Controls position-specific matching windows within genes

### Output

**cARS/PScARS scalar output** (when called directly):
- Dictionary mapping feature names to scalar scores
- `"cARS_score"` or `"PScARS_score"`: Adaptation score (higher = better adapted)
- Scores represent average maximal common substring length across all codon positions
- Typical range: ~10-50 for bacterial genes, depends on reference set and sequence length
- Units: average substring length in codons

**cARS/PScARS vector output** (when using PositionalFeature protocol):
- Dictionary mapping feature names to numpy arrays
- `"cARS_score"` or `"PScARS_score"`: Per-position maximal substring lengths
- Array length equals number of codons in the sequence
- Each element represents the longest substring starting at that codon position that appears in the reference set
- Units: substring length in codons at each codon position
- Aggregated to windows using mean (matching scalar definition)

**cMap/PScMap output** (not yet implemented):
- Optimized nucleotide sequences
- Block composition information
- Optimization metrics

### Modes

The feature can operate on:
- **Single sequences or batches**: Uses multiprocessing for batch operations
- **Different alphabets**: Nucleotide, codon (default), or amino acid
  - This implementation uses codon alphabet for cARS/PScARS
- **Global or position-specific matching**:
  - cARS: Global matching (any position in target can match any position in reference)
  - PScARS: Position-specific matching (positions constrained by `win_params` and `max_pos`)
- **Scalar or vector computation**:
  - Direct call (`feature(record)`): Returns scalar score (mean of per-position values)
  - PositionalFeature protocol (`feature.compute_vector(record)`): Returns per-position vector for windowing

## References

- Zur, H., & Tuller, T. (2015). Exploiting hidden information interleaved in the redundancy of the genetic code without prior knowledge. *Bioinformatics*, 31(9), 1398-1404. https://doi.org/10.1093/bioinformatics/btu797

- Diament, A., Schreiber, M., & Tuller, T. (2019). ChimeraUGEM: unsupervised gene expression modeling in any given organism. *Bioinformatics*, 35(13), 2243-2250. https://doi.org/10.1093/bioinformatics/btz080

- Burghardt, L., et al. (2025). Multi-sequence Chimera Map for optimizing gene expression in multi-copy systems. (In preparation)

## Upstream library links

- pyChimera GitHub: https://github.com/CompSynthBio/pyChimera
- ChimeraUGEM website: https://www.cs.tau.ac.il/~tamirtul/ChimeraUGEM/
- Zenodo DOI: https://doi.org/10.5281/zenodo.17577902

## Examples

### Basic Usage: cARS Prediction

```python
from biotooler.families.chimera import ChimeraFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Create reference set (e.g., highly expressed E. coli genes)
reference_genes = [
    "ATGAAACGCATTAGCACCACCATTACC...",  # Gene 1
    "ATGGCAAGCGTTATTAAAGGTGTTACC...",  # Gene 2
    # ... more reference genes
]

# Target gene to analyze
target_gene = "ATGACCGTTAAAGGTGTTACCATCGAA..."

# Create feature and compute cARS score
feature = ChimeraFeature(reference_genes, algorithm="cARS")
record = SeqRecord(Seq(target_gene), id="my_gene")
result = feature(record)
print(result)  # {'cARS_score': 35.2}
```

### Advanced: Position-Specific cARS

```python
# Use position-specific algorithm with window parameters
win_params = {'size': 40, 'center': 0, 'by_start': True, 'by_stop': True}
feature = ChimeraFeature(
    reference_genes,
    algorithm="PScARS",
    win_params=win_params,
    max_len=40,
    max_pos=0.5
)
result = feature(record)
print(result)  # {'PScARS_score': 38.7}
```

## Edge cases and validation

### Windowing Support

**Chimera features support windowing only via the PositionalFeature protocol**, not via direct slicing.

The cARS score is mathematically defined as:
```
cARS(target) = (1/L) * Σ(i=1 to L) max_length(target[i:], reference_set)
```

This formula averages the maximal common substring length **starting at every position** in the target sequence. **Slicing the target sequence before computing cARS would be mathematically incorrect** because:

1. The maximal common substrings are computed relative to the full target sequence context
2. Slicing changes which substrings are maximal at each position
3. This produces results that don't reflect the true cARS definition

**Correct windowing approach**: `ChimeraFeature` implements the `PositionalFeature` protocol to ensure windowing correctness:

- **`position_space`**: Returns `PositionSpace.CODON` (operates on codon-level sequences)
- **`vector_keys`**: Returns aggregation specifications using mean (matching the scalar cARS definition)
- **`compute_vector`**: Computes per-codon maximal common substring lengths across the **full sequence** using pyChimera's `calc_cARS(..., return_vec=True)` API

This design ensures that:
1. Per-position values are computed on the full sequence (preserving full context)
2. Windows select subsets of positions from the full-sequence vector
3. Selected positions are aggregated using mean (matching scalar cARS)
4. Windowing is correct-by-construction with no possibility of silent mathematical errors

**Example with windowing:**

```python
from biotooler.families.chimera import ChimeraFeature
from biotooler.features.sets import FeatureSet

# Create chimera feature
feature = ChimeraFeature(reference_genes, algorithm="cARS")

# Use with FeatureSet for windowed computation
feature_set = FeatureSet(features={"chimera": feature})

# Compute on ORF windows (uses PositionalFeature protocol internally)
results = feature_set.compute_orf_windows_v2(
    record, 
    orf_candidates, 
    window_size=50  # 50 codons per window
)
# Each window gets correct cARS score computed from per-position vectors
# windowing machinery automatically:
# 1. Calls compute_vector() once on full sequence
# 2. Slices vector to window boundaries
# 3. Aggregates using mean
```

### Validated behavior

- Lazy import behavior verified in tests
- Module structure follows biotooler family conventions
- Integration with lazy_import helper confirmed
- PositionalFeature protocol implementation verified in tests
- Per-position vector computation with `return_vec=True` tested
- Vector aggregation matches scalar computation (both use mean)

### Known limitations

- Only cARS and PScARS algorithms currently implemented (cMap/PScMap are future work)
- Requires pyChimera installation: `pip install git+https://github.com/CompSynthBio/pyChimera` or `pip install "biotooler[chimera]"`
- Suffix array preprocessing required for reference sequences (cached after first build)
- Optimal for coding sequences; less meaningful for non-coding regions
- Reference set quality and composition directly impact prediction accuracy
- Operates at codon level—requires sequences with lengths divisible by 3 for proper codon conversion
- **Windowing requires PositionalFeature protocol**—no fallback to direct slicing (by design, to prevent incorrect results)

## Windowing correctness

### Position space

This family computes features at the **codon position level**. The `position_space` property returns `PositionSpace.CODON`, indicating that each position corresponds to one codon (3 nucleotides).

**Implementation details:**
- Input DNA/RNA sequences are converted to codon representation using pyChimera's `nt2codon()` function
- A sequence of length 300 nucleotides becomes 100 codon positions
- Window sizes in biotooler's windowing APIs are interpreted as codon units when this family is used
- The Chimera algorithms internally work with codon-level sequences to capture synonymous codon usage patterns

### Vector computation

This family **does use the PositionalFeature protocol** to support correct windowing semantics.

**Why vector computation is necessary for Chimera:**

The cARS score is defined as an average over maximal common substring lengths at every position:
```
cARS(target) = (1/L) * Σ(i=1 to L) max_length(target[i:], reference_set)
```

The critical issue: `max_length(target[i:], reference_set)` depends on the **full suffix** starting at position i. If we slice the target sequence before computing cARS, we change the suffixes and get mathematically incorrect results.

**Solution using pyChimera API:**

The `compute_vector` method calls pyChimera's `calc_cARS` function with `return_vec=True`:

```python
# Build suffix array from reference sequences (once, cached)
ref_cod = nt2codon(self._normalized_refs)
self._suffix_array = build_suffix_array(ref_cod, pos_spec=("PS" in self.algorithm))

# Convert target to codons
target_cod = nt2codon([target_seq])

# Compute per-position vector (NOT scalar)
cars_vec = calc_cARS(
    target_cod[0],
    self._suffix_array,
    win_params=self.win_params if "PS" in self.algorithm else None,
    max_len=self.max_len,
    max_pos=self.max_pos,
    return_vec=True,  # Returns array of per-position values
)
```

This returns a numpy array where `cars_vec[i]` is the maximal common substring length starting at codon position i, computed in the context of the full target sequence.

**Upstream library API:**
- Function: `chimera.calc_cARS()`
- Parameter: `return_vec=True` (instead of default `return_vec=False`)
- Returns: Array of floats (per-position maximal substring lengths) instead of single float (mean)
- Documentation: Available in pyChimera GitHub repository (https://github.com/CompSynthBio/pyChimera)

### Aggregation strategy

The per-position vector values are aggregated into window-level scores using **mean** aggregation, matching the mathematical definition of the scalar cARS score.

**Aggregation specification:**

```python
@property
def vector_keys(self) -> dict[str, AggregationSpec]:
    feature_name = f"{self.algorithm}_score"
    return {
        feature_name: AggregationSpec(aggregation_fn=np.mean),
    }
```

**Windowing semantic (union-of-windows positions in full-context, then aggregate):**

1. **Full-context computation**: `compute_vector()` is called **once** on the entire target sequence, producing a vector of length L (number of codons)
2. **Union-of-windows positions**: For a set of windows W₁, W₂, ..., Wₙ covering positions in the sequence, each window Wᵢ = [start_i, end_i) identifies a subset of positions
3. **Per-window aggregation**: For each window Wᵢ:
   - Extract the vector slice: `vec_window = cars_vec[start_i:end_i]`
   - Aggregate using mean: `score_window = np.mean(vec_window)`
4. **Result**: Each window gets a score that represents the average maximal common substring length for positions in that window, computed with full sequence context

**Why this is correct:**

- Every position's value reflects the maximal substring starting at that position **in the full sequence**
- Windows select subsets of these full-context position values
- Aggregation (mean) matches the scalar cARS definition: mean of per-position maximal lengths
- No recomputation needed for overlapping windows (efficiency)
- Mathematically equivalent to computing cARS on each window **if** the reference set and sequence context were preserved (which the vector approach ensures)

**Example:**

```python
# Sequence: 300 nt = 100 codons
# Full-context vector: [35.2, 28.1, 42.7, ..., 31.5]  # length 100

# Window 1: codons 0-49
# Vector slice: [35.2, 28.1, 42.7, ..., (50 values)]
# Aggregated score: mean([35.2, 28.1, 42.7, ...]) = 33.4

# Window 2: codons 25-74 (overlaps with Window 1)
# Vector slice: [(from position 25), ..., (50 values)]
# Aggregated score: mean([...]) = 36.8
# Note: No recomputation needed, just different slice of same vector
```

### Testing approach

Windowing correctness for this family is validated through:

1. **Full sequence context preservation**: 
   - Verify that `compute_vector()` receives the complete input sequence
   - Confirm that pyChimera's `calc_cARS` is called with full target sequence (not sliced)
   - Check that suffix array is built from complete reference sequences

2. **Vector vs scalar consistency**:
   - Compare scalar cARS (direct call with `return_vec=False`) against mean of vector (call with `return_vec=True`)
   - Both should produce identical results for full-sequence computation
   - Test formula: `scalar_score ≈ np.mean(vector_values)`

3. **Window aggregation correctness**:
   - Verify that window slicing produces correct sub-vectors
   - Confirm mean aggregation matches expected values
   - Test with overlapping windows to ensure no recomputation artifacts

4. **Reference set consistency**:
   - Ensure suffix array is built once and reused across all computations
   - Verify that reference sequences remain constant for all window computations
   - Check that `max_len`, `max_pos`, and `win_params` are correctly passed to pyChimera

5. **Algorithm parameter validation**:
   - Test both cARS (global) and PScARS (position-specific) algorithms
   - Verify that PScARS uses `win_params` correctly
   - Confirm that position-specific matching constraints are applied

6. **Edge cases**:
   - Empty sequences (should return empty vectors)
   - Short sequences (< max_len codons)
   - Sequences not divisible by 3 (codon conversion handling)
   - Windows at sequence boundaries

Tests validate that the PositionalFeature protocol implementation ensures mathematically correct windowing that matches the cARS definition while avoiding silent errors from naive slicing approaches.

## Maintenance notes

### Dependencies

- **pyChimera**: Python implementation of Chimera algorithms (optional, installed via `pip install "biotooler[chimera]"`)
- **BioPython**: For SeqRecord handling (core biotooler dependency)
- **NumPy**: Required by pyChimera for suffix array operations

### Future enhancements

- Implement cARS score computation
- Implement PScARS with position-specific scoring
- Add cMap/PScMap optimization features
- Support for caching suffix arrays
- Batch processing optimization
- Support for different sequence alphabets (nt, codon, aa)
- Integration with codon usage tables
- Homolog filtering options

### Testing strategy

- Lazy import tests ensure pyChimera is not loaded on module import
- Feature instantiation tests verify proper error messages without pyChimera
- PositionalFeature protocol tests verify correct implementation of position_space, vector_keys, and compute_vector
- Integration tests will validate cARS scores against known examples
- Performance tests for suffix array construction and querying
- Windowing tests will verify correct aggregation of per-position vectors
