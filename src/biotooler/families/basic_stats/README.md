# Basic Statistics Feature Family

## What this family provides

The basic_stats family provides fundamental sequence composition analysis for DNA, RNA, and protein sequences. It computes nucleotide/amino acid counts, frequencies, sequence length, and GC content (for DNA/RNA). The family supports both baseline per-window computation and optimized incremental computation using sliding window count vectors for efficient analysis of large sequences.

## Intuition

Understanding the basic composition of biological sequences is a foundational step in sequence analysis. The distribution of nucleotides or amino acids provides insights into:

1. **Sequence identity**: DNA vs RNA (T vs U) vs protein sequences
2. **GC content**: A key indicator of genome structure, gene expression, and evolutionary adaptation
   - High GC content is associated with thermophilic organisms, stable genomic regions, and coding sequences
   - Low GC content is common in AT-rich intergenic regions and some viral genomes
3. **Amino acid composition**: Protein properties (hydrophobicity, charge, size distribution)
4. **Sequence bias**: Over/under-representation of specific nucleotides or amino acids

Basic statistics serve as:
- **Quality control metrics** for sequence data
- **Input features** for machine learning models
- **Filtering criteria** for sequence selection
- **Baseline comparisons** for more complex analyses

This family provides these statistics efficiently, supporting both small-scale single-window analysis and large-scale sliding window studies where incremental computation significantly reduces computational overhead.

## Mathematical formulation

### Core Statistics

For a sequence S of length L with alphabet Σ (e.g., {A,C,G,T} for DNA):

**Count of character c:**
```
count(c) = |{i : S[i] = c}|
```

**Fraction of character c:**
```
fraction(c) = count(c) / L     if L > 0
            = 0                 if L = 0
```

**Sequence length:**
```
length = L = |S|
```

### GC Content (DNA/RNA only)

For DNA sequences with alphabet {A, C, G, T} or RNA sequences with alphabet {A, C, G, U}:

```
gc_fraction = (count(G) + count(C)) / L     if L > 0
            = 0                              if L = 0
```

**Biological interpretation:**
- gc_fraction ∈ [0, 1]
- gc_fraction = 0.5 indicates equal representation of GC vs AT/AU
- gc_fraction > 0.5 indicates GC-rich regions (common in coding sequences, promoters)
- gc_fraction < 0.5 indicates AT/AU-rich regions (common in intergenic regions)

### Incremental Computation

For sliding windows with step size s, instead of recomputing full counts for each window, we maintain a count vector **c** ∈ ℕ^|Σ| and update incrementally:

**Initialize** for first window [w_start, w_end):
```
c[k] = count of character Σ[k] in S[w_start:w_end]
```

**Update** for next window (remove S[out_start:out_end], add S[in_start:in_end]):
```
for each character σ in S[out_start:out_end]:
    c[index(σ)] -= 1

for each character σ in S[in_start:in_end]:
    c[index(σ)] += 1
```

This reduces time complexity from O(W × N) for baseline computation (W = window size, N = number of windows) to O(L) for incremental computation (L = sequence length), achieving speedups proportional to W/s when s << W.

## Features and output schema

### Input
- No configuration required - feature automatically detects sequence type

### Output
Dictionary mapping feature names to scalar values:

**For DNA sequences (alphabet: A, C, G, T, optionally N):**
```python
{
    "length": int,           # Sequence/window length
    "count_a": int,          # Count of A
    "count_c": int,          # Count of C
    "count_g": int,          # Count of G
    "count_t": int,          # Count of T
    "count_n": int,          # Count of N (if present)
    "fraction_a": float,     # Fraction of A (count_a / length)
    "fraction_c": float,     # Fraction of C
    "fraction_g": float,     # Fraction of G
    "fraction_t": float,     # Fraction of T
    "fraction_n": float,     # Fraction of N (if present)
    "gc_fraction": float     # (count_g + count_c) / length
}
```

**For RNA sequences (alphabet: A, C, G, U, optionally N):**
```python
{
    "length": int,
    "count_a": int,
    "count_c": int,
    "count_g": int,
    "count_u": int,          # Count of U (not T)
    "count_n": int,          # Count of N (if present)
    "fraction_a": float,
    "fraction_c": float,
    "fraction_g": float,
    "fraction_u": float,
    "fraction_n": float,     # Fraction of N (if present)
    "gc_fraction": float     # (count_g + count_c) / length
}
```

**For protein sequences (alphabet: 20 standard amino acids, optionally X, *):**
```python
{
    "length": int,
    "count_a": int,          # Alanine
    "count_c": int,          # Cysteine
    # ... (all 20 standard amino acids)
    "count_x": int,          # Unknown amino acid (if present)
    "count_*": int,          # Stop codon (if present)
    "fraction_a": float,
    "fraction_c": float,
    # ... (all amino acids)
    # Note: NO gc_fraction for protein sequences
}
```

### Modes
1. **Baseline mode**: Calls `_compute_stats()` for each window independently
2. **Incremental mode**: Uses `init_state/step_state/emit` with numpy count vectors for sliding windows

Both modes produce identical results (exact match for counts, within float tolerance for fractions).

## References

### GC content and genome biology
- Bernardi G (2000) Isochores and the evolutionary genomics of vertebrates. Gene 241(1):3-17. https://doi.org/10.1016/S0378-1119(99)00485-0
- Vinogradov AE (2003) DNA helix: the importance of being GC-rich. Nucleic Acids Research 31(7):1838-1844. https://doi.org/10.1093/nar/gkg296

### Nucleotide composition and gene expression
- Kudla G, Murray AW, Tollervey D, Plotkin JB (2009) Coding-sequence determinants of gene expression in Escherichia coli. Science 324(5924):255-258. https://doi.org/10.1126/science.1170160

### Amino acid composition analysis
- Nakashima H, Nishikawa K, Ooi T (1986) The folding type of a protein is relevant to the amino acid composition. Journal of Biochemistry 99(1):153-162. https://doi.org/10.1093/oxfordjournals.jbchem.a135454

### Sliding window algorithms
- Welch LR (1984) Hidden Markov models and the Baum-Welch algorithm. IEEE Information Theory Society Newsletter 34(4):1-13. (General sliding window concepts)

## Upstream library links

N/A - This family uses only standard Python libraries (numpy) and Biopython for sequence handling. No external domain-specific libraries are required.

## Examples

### Basic usage with DNA sequence

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.basic_stats import BasicStatsFeature

# Create feature
feature = BasicStatsFeature()

# Analyze DNA sequence
dna_record = SeqRecord(Seq("ACGTACGTGGGCCC"), id="dna1")
result = feature(dna_record)

print(f"Length: {result['length']}")           # 14
print(f"A count: {result['count_a']}")         # 2
print(f"GC fraction: {result['gc_fraction']}")  # 0.714... (10/14)
```

### RNA sequence detection

```python
# RNA sequence (contains U instead of T)
rna_record = SeqRecord(Seq("ACGUACGUGGGCCC"), id="rna1")
result = feature(rna_record)

print(f"U count: {result['count_u']}")         # 2
print(f"GC fraction: {result['gc_fraction']}")  # 0.714...
# Note: 'count_t' is not in the result for RNA
```

### Protein sequence analysis

```python
# Protein sequence
protein_record = SeqRecord(Seq("MKALVSWGRDE"), id="protein1")
result = feature(protein_record)

print(f"Length: {result['length']}")           # 11
print(f"Alanine count: {result['count_a']}")   # 1
print(f"Glycine fraction: {result['fraction_g']}")  # 0.0909... (1/11)
# Note: 'gc_fraction' is not in the result for protein
```

### Incremental computation with FeatureSet

```python
from biotooler.features import FeatureSet

# Create feature set
fs = FeatureSet(BasicStatsFeature(), name="stats")

# Compute over sliding windows (incremental)
dna_record = SeqRecord(Seq("ACGTACGTACGTACGT"), id="seq1")
result_df = fs.compute_windows(
    dna_record,
    window_size=8,
    step=4
)

# Wide format output with window-specific columns
print(result_df.columns)
# ['record_id', 'region_start', 'region_end',
#  'stats.length_0', 'stats.count_a_0', ..., 'stats.gc_fraction_0',
#  'stats.length_4', 'stats.count_a_4', ..., 'stats.gc_fraction_4',
#  'stats.length_8', ...]
```

### Comparing baseline vs incremental

```python
# Both methods produce identical results
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.basic_stats import BasicStatsFeature
from biotooler.features import FeatureSet
from biotooler.core.windowing import iter_windows

feature = BasicStatsFeature()
record = SeqRecord(Seq("ACGTACGTACGTACGT"), id="test")

# Baseline: compute each window independently
baseline_results = []
for window in iter_windows(record, window_size=8, step=4):
    baseline_results.append(feature(window))

# Incremental: use FeatureSet.compute_windows
fs = FeatureSet(feature, name="stats")
incremental_df = fs.compute_windows(record, window_size=8, step=4)

# Results match exactly
print(baseline_results[0]["count_a"])  # 2
print(incremental_df["stats.count_a_0"].values[0])  # 2
```

## Edge cases and validation

### Validated behavior

1. **Empty sequences**: Returns length=0, all counts=0, all fractions=0.0
2. **Single-character sequences**: Correctly computes fractions as 1.0 for that character
3. **Ambiguous nucleotides**: 
   - DNA with N: Counts N separately, includes N in length but NOT in gc_fraction
   - RNA with N: Same behavior as DNA
4. **Protein with X or ***: Counts separately, includes in length and fractions
5. **Sequence type detection**:
   - Presence of U → RNA
   - Presence of T (without U) → DNA
   - Presence of protein-specific amino acids → protein
   - Ambiguous cases (only ACGN) → defaults to DNA
6. **Incremental vs baseline**: Both methods produce identical counts and fractions within floating-point tolerance (< 1e-10 relative error)
7. **Alphabet consistency**: Alphabet is detected from full sequence, not individual windows, ensuring consistent feature keys across all windows

### Known limitations

1. **Modified bases**: Does not handle modified nucleotides (e.g., methylated cytosine) or non-standard amino acids beyond X and *
2. **Mixed case**: Input sequences are converted to uppercase by `get_seq_str`; lowercase letters in original sequence are not distinguished
3. **Ambiguous IUPAC codes**: Only N is explicitly handled for DNA/RNA; other ambiguity codes (R, Y, S, W, K, M, etc.) are treated as unknown characters
4. **Sequence type ambiguity**: Sequences containing only ACGN without T or U are assumed to be DNA; there is no way to force RNA interpretation without including U
5. **Codon-aware windowing**: This feature does not enforce codon boundaries; use with `compute_orf_windows` if codon alignment is required

### Performance characteristics

- **Baseline mode**: O(W) per window, O(W × N) total for N windows
- **Incremental mode**: O(W) initialization + O(s) per step, O(W + s × N) total
- **Speedup**: When s << W (e.g., step=1, window=100), incremental mode is ~W/s times faster
- **Memory**: O(|Σ|) for count vector (typically 4-5 for DNA/RNA, 20-22 for protein)

## Windowing correctness

### Position space

This family does not use the PositionalFeature protocol. Instead, it uses **incremental computation** for efficient sliding window analysis. Features are computed directly on window-level sequences rather than per-position values.

### Vector computation

Not applicable - this family uses the IncrementalFeature protocol with `init_state`, `step_state`, and `emit` methods rather than the PositionalFeature protocol with `compute_vector`.

The incremental approach:
1. `init_state`: Initializes a count vector for the first window by iterating through the window sequence
2. `step_state`: Updates the count vector by subtracting outgoing characters and adding incoming characters
3. `emit`: Computes feature values (counts and fractions) from the current count vector state

This approach ensures correctness by maintaining exact character counts and computing fractions on demand.

### Aggregation strategy

Not applicable - features are computed directly on each window rather than aggregating per-position values. Each window produces:
- Exact character counts (integers)
- Character fractions (computed as count/length)
- GC fraction for DNA/RNA sequences

### Testing approach

Windowing correctness is validated through:
1. **Baseline vs incremental comparison**: Tests verify that incremental computation produces identical results to baseline (per-window) computation
2. **Edge case handling**: Empty windows, single-character windows, and full-sequence windows
3. **Alphabet consistency**: Alphabet detected from full sequence ensures consistent feature keys across all windows
4. **Numerical precision**: Floating-point comparisons use tolerance < 1e-10 for fraction calculations

See `tests/families/test_basic_stats.py` for complete test coverage validating both modes produce identical results.

## Maintenance notes

### Dependencies
- `numpy`: For efficient count vector operations in incremental mode
- `biopython`: For SeqRecord handling
- Internal `biotooler.core.seq_utils`: For sequence string extraction with caching
- Internal `biotooler.core.types`: For Scalar type definition

### Design decisions
1. **Alphabet auto-detection**: Simplifies API by not requiring explicit sequence type parameter, but may cause issues with ambiguous sequences
2. **Lowercase feature names**: Follows convention of using lowercase for consistency (count_a, not count_A)
3. **Separate count and fraction**: Provides both absolute and relative measures for flexibility in downstream analysis
4. **GC-only for DNA/RNA**: Protein sequences do not have a gc_fraction to avoid confusion
5. **N handling**: N is counted separately but excluded from gc_fraction to match standard biology practice

### Future enhancements
- Support for custom alphabets (e.g., reduced amino acid alphabets)
- IUPAC ambiguity code handling with fractional counting
- Configurable alphabet detection (force DNA/RNA/protein interpretation)
- Support for dinucleotide/dipeptide frequencies
- Batch computation across multiple sequences
- Vectorized incremental computation for multiple features simultaneously
