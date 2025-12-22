# Chimera Feature Family

## What this family provides

The chimera family provides feature computation for gene expression prediction and sequence optimization using the Chimera algorithms. These algorithms can predict the expression level of a gene in an unsupervised manner based solely on the coding sequence and the host genome. The algorithms measure sequence adaptation by comparing target genes against reference gene sets to quantify similarity patterns.

## Intuition

Gene expression levels are influenced not only by regulatory elements but also by the coding sequence itself. The Chimera algorithms exploit "hidden information" embedded in the redundancy of the genetic code - information beyond what is encoded in the amino acid sequence. This information can reflect evolutionary optimization for translation efficiency, mRNA stability, and other expression-related factors.

The core insight is that highly expressed genes in a host organism tend to share common sequence patterns. By measuring how well a target gene's sequence matches patterns found in a reference set of host genes, we can:

1. **Predict expression levels**: Genes with higher similarity to reference gene patterns are predicted to express better
2. **Optimize sequences**: Design synonymous variants that maximize similarity to reference patterns while maintaining the same protein
3. **Detect adaptation**: Identify whether genes are well-adapted to their host organism

The Chimera algorithms work by finding maximal common substrings between target and reference sequences, either globally (cARS/cMap) or with position-specific constraints (PScARS/PScMap).

## Mathematical formulation

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

### Implementation approach

The algorithms use suffix arrays for efficient substring matching:
1. Build suffix array from reference sequences (one-time preprocessing)
2. Query target sequences against the suffix array
3. For cARS: calculate average maximal match lengths
4. For cMap: greedily select optimal codon choices

## Features and output schema

### Input

- **Reference sequences**: Set of host genes (DNA sequences) used as the adaptation template
- **Target sequences**: Query genes to analyze or optimize (DNA for cARS, amino acid for cMap)
- **Parameters**:
  - `max_len`: Maximum substring length to consider (homolog filtering)
  - `max_pos`: Maximum position difference (for position-specific variants)
  - `win_params`: Window parameters for position-specific algorithms

### Output

**cARS/PScARS output**:
- Dictionary mapping feature names to scalar scores
- `"cARS_score"` or `"PScARS_score"`: Adaptation score (higher = better adapted)
- Scores typically range from ~10-50 for typical genes

**cMap/PScMap output** (future):
- Optimized nucleotide sequences
- Block composition information
- Optimization metrics

### Modes

The feature can operate on:
- Single sequences or batches (uses multiprocessing)
- Different alphabets: nucleotide, codon, or amino acid
- Global or position-specific matching

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

### Validated

- Lazy import behavior verified in tests
- Module structure follows biotooler family conventions
- Integration with lazy_import helper confirmed

### Known limitations

- Feature computation not yet implemented (stub only)
- Requires pyChimera installation: `pip install git+https://github.com/CompSynthBio/pyChimera`
- Suffix array preprocessing required for reference sequences (can be cached)
- Optimal for coding sequences; less meaningful for non-coding regions
- Reference set quality impacts prediction accuracy

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
- Integration tests will validate cARS scores against known examples
- Performance tests for suffix array construction and querying
