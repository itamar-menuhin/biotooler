# Codon Bias Feature Family

## What this family provides

The codon bias family provides feature computation for analyzing codon usage bias in DNA and RNA sequences. It wraps external codonbias package models (CAI, ENC, FOP, RSCU, RCBS/DCBS, tAI, nTE, CPB/CPS) to compute scalar features on biological sequences. The family supports both baseline per-window computation and optimized rolling/incremental computation for efficient analysis of large sequences.

## Intuition

Codon usage bias refers to the non-uniform usage of synonymous codons in coding sequences. The genetic code is degenerate—most amino acids can be encoded by multiple codons (synonymous codons). However, organisms do not use these synonymous codons with equal frequency. This bias arises from multiple evolutionary pressures including:

1. **Translation efficiency**: Preferred codons often correspond to abundant tRNAs, enabling faster translation
2. **Translation accuracy**: Optimal codons can reduce amino acid misincorporation errors
3. **mRNA stability**: Codon choice affects secondary structure and degradation rates
4. **Gene expression regulation**: Rare codons can serve as regulatory checkpoints
5. **Protein folding**: Translation speed variations allow proper co-translational folding

Different organisms, tissues, and genes exhibit distinct codon preferences. Highly expressed genes (e.g., ribosomal proteins, glycolysis enzymes) typically show strong bias toward "optimal" codons that match abundant tRNAs. Understanding codon bias enables:
- **Prediction of gene expression levels** from sequence alone
- **Optimization of heterologous protein expression** (e.g., humanizing sequences for E. coli)
- **Identification of horizontally transferred genes** (foreign codon usage patterns)
- **Inference of translation dynamics** and evolutionary constraints

The metrics in this family quantify different aspects of codon bias, from overall bias strength (ENC) to adaptation to specific expression contexts (CAI, tAI).

## Mathematical formulation

This family implements multiple codon bias metrics from the codonbias package. Each metric captures different biological properties of codon usage patterns.

### CAI (Codon Adaptation Index)

**Purpose**: Measures how well a gene's codon usage matches that of highly expressed reference genes.

**Formula**: CAI is the geometric mean of relative adaptiveness (w) values for each codon in the sequence:

```
CAI = exp((1/L) * Σ ln(w_i))
```

where:
- L = number of codons in the sequence
- w_i = relative adaptiveness of codon i
- w is calculated from a reference set of highly expressed genes

**Relative adaptiveness** for a codon encoding amino acid X:
```
w(codon) = f(codon) / f_max(X)
```
where f(codon) is the frequency of that codon in the reference set, and f_max(X) is the frequency of the most common codon for amino acid X.

**Range**: 0 to 1, where 1 indicates perfect adaptation to reference codon usage.

**Interpretation**:
- CAI > 0.8: High expression potential, well-adapted codon usage
- CAI 0.5-0.8: Moderate expression potential
- CAI < 0.5: Poor adaptation, likely low expression

**Reference**: Sharp & Li (1987) Nucleic Acids Research 15:1281-1295

### ENC (Effective Number of Codons)

**Purpose**: Quantifies the overall level of codon bias independent of any reference set. Measures how many codons are effectively used in a sequence.

**Formula**: ENC estimates the effective number of codons used by calculating the average homozygosity for each amino acid family:

```
ENC = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆
```

where F̄ₖ is the average homozygosity for amino acid families with k synonymous codons:

```
F̄ₖ = Σ_families [(Σ_codons (n_i/n)²) - 1/(n-1)] / (k-1)
```

where n_i is the count of codon i and n is the total count for that amino acid family.

**Range**: 20 to 61
- 20: Extreme bias (only one codon per amino acid)
- 61: No bias (all synonymous codons used equally)

**Interpretation**:
- ENC < 35: Strong codon bias (typical for highly expressed genes)
- ENC 35-50: Moderate bias
- ENC > 50: Weak bias (random codon usage)

**Advantage**: Does not require a reference set; measures intrinsic bias.

**Reference**: Wright (1990) Gene 87:23-29

### FOP (Frequency of Optimal Codons)

**Purpose**: Measures the proportion of codons in a sequence that are "optimal" based on a reference set.

**Formula**:
```
FOP = N_optimal / N_total
```

where:
- N_optimal = count of optimal codons (those used most frequently in reference set)
- N_total = total number of codons (excluding Met, Trp, stop codons)

**Range**: 0 to 1

**Interpretation**:
- FOP > 0.7: High use of optimal codons, likely highly expressed
- FOP 0.4-0.7: Moderate optimality
- FOP < 0.4: Low use of optimal codons

**Comparison to CAI**: FOP is simpler (binary optimal/non-optimal) while CAI uses weighted relative adaptiveness.

**Reference**: Ikemura (1981) Journal of Molecular Biology 151:389-409

### RSCU (Relative Synonymous Codon Usage)

**Purpose**: Measures the observed frequency of a codon relative to its expected frequency under uniform usage.

**Formula**: For each codon encoding amino acid X:
```
RSCU(codon) = (observed_frequency(codon) × n_synonyms(X)) / Σ_synonyms observed_frequency
```

where n_synonyms(X) is the number of synonymous codons for amino acid X.

**Range**: 
- RSCU < 1: Codon used less than expected (under-represented)
- RSCU = 1: Codon used at expected frequency (no bias)
- RSCU > 1: Codon used more than expected (over-represented)

**Interpretation**: RSCU values reveal which specific codons are preferred or avoided. A value of 6 for a 6-fold degenerate amino acid means that codon is used exclusively.

**Application**: Identifying codon preferences for individual codons rather than genome-wide bias.

**Reference**: Sharp et al. (1986) Nucleic Acids Research 14:5125-5143

### tAI (tRNA Adaptation Index)

**Purpose**: Measures adaptation to the tRNA pool by accounting for tRNA gene copy numbers and wobble base pairing rules.

**Formula**:
```
tAI = exp((1/L) * Σ ln(W_i))
```

where:
- L = number of codons
- W_i = absolute adaptiveness of codon i, calculated from tRNA gene copy numbers

**Wobble pairing weights** (from Crick's wobble hypothesis):
- Watson-Crick pairs: weight = 1.0
- G:U wobble: weight = 0.5-0.8
- I:U, I:C, I:A wobble: weight varies

**Range**: 0 to 1, similar to CAI but based on tRNA availability rather than codon frequency.

**Advantage**: Biologically more direct—directly relates to translation machinery rather than inferring from highly expressed genes.

**Reference**: dos Reis et al. (2004) Nucleic Acids Research 32:5036-5044

### nTE (normalized Translation Efficiency)

**Purpose**: Estimates translation efficiency by combining tRNA abundance, wobble pairing, and other factors affecting translation speed.

**Formula**: nTE incorporates:
- tRNA availability (gene copy numbers)
- Codon-anticodon affinity (including wobble penalties)
- Position-dependent effects
- Normalization across sequences

The exact formula is implementation-specific but generally:
```
nTE = (TE_observed - TE_min) / (TE_max - TE_min)
```

where TE is a weighted sum of codon translation times.

**Range**: 0 to 1, with higher values indicating more efficient translation.

**Applications**: 
- Predicting protein abundance from mRNA levels
- Engineering sequences for optimal expression

**Reference**: Carbone et al. (2003) BMC Bioinformatics 4:10

### RCBS/DCBS (Relative/Distance from Codon Bias Strength)

**Purpose**: Measures the strength of codon bias as deviation from uniform usage.

**RCBS Formula**:
```
RCBS = (Σ_amino_acids Σ_codons |RSCU - 1|) / (2 × N_codons)
```

**Range**: 0 to 1
- 0: No bias (all synonymous codons equally frequent)
- 1: Maximum bias (one codon per amino acid)

**DCBS**: Distance-based measure comparing observed codon frequencies to expected under no bias.

**Interpretation**: Unlike ENC, these metrics are linear and easier to compare across sequences of different lengths.

### CPB/CPS (Codon Pair Bias/Score)

**Purpose**: Measures bias in the usage of adjacent codon pairs, which can affect translation efficiency and mRNA stability beyond individual codon effects.

**Concept**: Some codon pairs are over- or under-represented relative to expectation from individual codon frequencies. This can result from:
- mRNA secondary structure preferences
- tRNA recycling efficiency
- Evolutionary optimization of translation

**Formula** (simplified):
```
CPB(codon_i, codon_j) = ln((observed_pair_freq) / (expected_pair_freq))
CPS = average CPB across all adjacent codon pairs
```

**Range**: CPS is typically centered near 0
- CPS > 0: Over-representation of common codon pairs
- CPS < 0: Under-representation (may indicate engineered or atypical sequences)

**Application**: Codon pair deoptimization for vaccine design (maintaining amino acids but reducing translation efficiency).

**Reference**: Gutman & Hatfield (1989) PNAS 86:3699-3703

### Implementation Notes

The `CodonBiasFeature` class implements the `PositionalFeature` protocol for modern windowing (v2), which provides efficient and accurate codon bias computation using **full-context vectors with window aggregation**:

#### Windowing v2 (PositionalFeature Protocol)

For scores that support `get_vector()` (CAI, FOP, RSCU, RCBS, CPB):

1. **Full-context vector computation**: Each score's `get_vector(sequence)` is called **once** on the entire sequence to compute per-codon values with full context
2. **Window aggregation**: Per-codon values are aggregated into windows using the appropriate aggregation function:
   - **CAI, tAI**: Geometric mean (as per mathematical definition: `exp(mean(log(weights)))`)
   - **FOP, RSCU, RCBS, CPB**: Arithmetic mean
3. **No sequence slicing**: The sequence is **never sliced**—each codon's value is computed with awareness of the full sequence context

This approach ensures that:
- Each positional value has full sequence context (e.g., for RSCU calculations that depend on codon frequencies)
- Window aggregation is mathematically correct (e.g., geometric mean for CAI)
- Computation is efficient (single vector computation, then efficient array slicing for windows)

**Example**:
```python
# For a sequence with 10 codons:
# 1. Compute full-context vector: cai.get_vector(seq) -> [w0, w1, w2, ..., w9]
# 2. For window covering codons 2-5: geometric_mean([w2, w3, w4, w5])
# 3. No slicing of the original sequence—values already computed with full context
```

#### Legacy Incremental Mode (IncrementalFeature)

For scores without `get_vector()` (like ENC):

1. **Incremental codon counting**: Maintains rolling codon counts as windows slide
2. **Fallback computation**: Uses `get_score()` directly when incremental computation is not possible
3. **Baseline slicing**: As a last resort, calls `get_score(seq, slice=slice(start, end))` for each window

Models like ENC that require holistic computation (considering all codon families) use this legacy mode.

## Features and output schema

### Input
- `models`: Sequence of `codonbias.scores.ScalarScore` instances
- `names`: Optional sequence of custom names for each model (defaults to class names)

### Output
Dictionary mapping feature names to scalar float values:
```python
{
    "CAI": 0.85,
    "ENC": 45.2,
    "FOP": 0.62
}
```

### Windowing Modes

1. **Windowing v2** (`compute_orf_windows_v2`): Uses `PositionalFeature` protocol
   - Full-context vector computation via `get_vector()`
   - Window aggregation with proper aggregation functions
   - Works for: CAI, FOP, RSCU, RCBS, CPB
   
2. **Legacy windowing** (`compute_orf_windows`): Uses `IncrementalFeature` protocol
   - Incremental codon counting and rolling computation
   - Fallback to baseline `get_score()` when needed
   - Works for: All scores, including ENC

## References

### Primary metric references
- **CAI**: Sharp PM, Li WH (1987) The codon adaptation index - a measure of directional synonymous codon usage bias, and its potential applications. Nucleic Acids Research 15(3):1281-1295. https://doi.org/10.1093/nar/15.3.1281
- **ENC**: Wright F (1990) The 'effective number of codons' used in a gene. Gene 87(1):23-29. https://doi.org/10.1016/0378-1119(90)90491-9
- **FOP**: Ikemura T (1981) Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. Journal of Molecular Biology 151(3):389-409. https://doi.org/10.1016/0022-2836(81)90003-6
- **RSCU**: Sharp PM, Tuohy TM, Mosurski KR (1986) Codon usage in yeast: cluster analysis clearly differentiates highly and lowly expressed genes. Nucleic Acids Research 14(13):5125-5143. https://doi.org/10.1093/nar/14.13.5125
- **tAI**: dos Reis M, Savva R, Wernisch L (2004) Solving the riddle of codon usage preferences: a test for translational selection. Nucleic Acids Research 32(17):5036-5044. https://doi.org/10.1093/nar/gkh834
- **nTE**: Carbone A, Zinovyev A, Képès F (2003) Codon adaptation index as a measure of dominating codon bias. Bioinformatics 19(16):2005-2015. https://doi.org/10.1093/bioinformatics/btg272
- **CPB**: Gutman GA, Hatfield GW (1989) Nonrandom utilization of codon pairs in Escherichia coli. PNAS 86(10):3699-3703. https://doi.org/10.1073/pnas.86.10.3699

### Foundational codon bias research
- Ikemura T (1985) Codon usage and tRNA content in unicellular and multicellular organisms. Molecular Biology and Evolution 2(1):13-34. https://doi.org/10.1093/oxfordjournals.molbev.a040335
- Bulmer M (1991) The selection-mutation-drift theory of synonymous codon usage. Genetics 129(3):897-907. https://doi.org/10.1093/genetics/129.3.897

### Reviews and applications
- Plotkin JB, Kudla G (2011) Synonymous but not the same: the causes and consequences of codon bias. Nature Reviews Genetics 12(1):32-42. https://doi.org/10.1038/nrg2899
- Quax TE, Claassens NJ, Söll D, van der Oost J (2015) Codon Bias as a Means to Fine-Tune Gene Expression. Molecular Cell 59(2):149-161. https://doi.org/10.1016/j.molcel.2015.05.035

## Upstream library links

- **Official documentation**: https://codon-bias.readthedocs.io/en/latest/
- **PyPI package**: https://pypi.org/project/codon-bias/
- **API reference**: https://codon-bias.readthedocs.io/en/latest/api.html

## Examples

```python
from codonbias.scores import CodonAdaptationIndex, EffectiveNumberOfCodons
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.codon_bias import CodonBiasFeature

# Create reference sequence for CAI
ref_seq = "ATGATGATGATGATGATGATGATGATGATGATG"

# Initialize models
cai = CodonAdaptationIndex(ref_seq)
enc = EffectiveNumberOfCodons()

# Create feature with custom names
feature = CodonBiasFeature([cai, enc], names=["CAI", "ENC"])

# Compute on a sequence
test_seq = "ATGATGATGATGATGATG"
record = SeqRecord(Seq(test_seq), id="test")
result = feature(record)

print(result)  # {'CAI': 1.0, 'ENC': 20.0}
```

## Edge cases and validation

### Validated
- RNA sequences are automatically converted to DNA before analysis
- Protein sequences raise `ValueError`
- Window boundaries are validated against sequence length
- Model and name list lengths must match
- Empty sequences are handled by underlying codonbias models
- Non-triplet sequences (incomplete codons) are handled by codonbias logic

### Known limitations
- Rolling mode requires accessible model weights; falls back to baseline otherwise
- Sequences must be coding sequences in correct reading frame for meaningful results
- Very short sequences (<3 codons) may produce unreliable statistics

## Maintenance notes

### Dependencies
- `codonbias` package for all score implementations
- `biopython` for SeqRecord handling
- Internal `biotooler.core.seq_utils` for sequence conversion

### Future enhancements
- Support for custom codon tables (e.g., mitochondrial, ciliate)
- Batch computation across multiple sequences
- Additional codon bias metrics as they become available in codonbias
- Caching of reference sequence weights for repeated CAI computations
