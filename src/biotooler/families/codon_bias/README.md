# Codon Bias Feature Family

## What this family provides

The codon bias family provides feature computation for analyzing codon usage bias in DNA and RNA sequences. It wraps external codonbias package models (CAI, ENC, FOP, RSCU, RCBS/DCBS, tAI, nTE, CPB/CPS) to compute scalar features on biological sequences. The family supports both baseline per-window computation and optimized rolling/incremental computation for efficient analysis of large sequences.

## Intuition

Codon usage bias refers to the non-uniform usage of synonymous codons in coding sequences. Different organisms, tissues, and genes exhibit distinct codon preferences, which can affect translation efficiency, protein expression levels, and gene regulation. By quantifying codon bias through metrics like CAI (Codon Adaptation Index) and ENC (Effective Number of Codons), researchers can gain insights into gene expression patterns, evolutionary pressures, and optimal codon usage for heterologous protein expression.

## Mathematical formulation

The family computes various codon bias scores using the codonbias package:
- **CAI (Codon Adaptation Index)**: Measures adaptation to codon usage of highly expressed genes using geometric mean of relative adaptiveness values
- **ENC (Effective Number of Codons)**: Quantifies overall codon usage bias, ranges from 20 (extreme bias) to 61 (no bias)
- **FOP (Frequency of Optimal codons)**: Fraction of codons that are optimal codons in highly expressed genes
- **RSCU (Relative Synonymous Codon Usage)**: Ratio of observed codon frequency to expected frequency under uniform usage
- Other metrics: RCBS/DCBS (codon bias strength), tAI (tRNA adaptation index), nTE (normalized translation efficiency), CPB/CPS (codon pair bias/score)

Each model implements `get_score(sequence)` to return a scalar value. For rolling computation, the family maintains codon counts and incrementally updates scores using model-specific weights when accessible.

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

### Modes
1. **Baseline mode**: Calls `model.get_score(seq_str, slice=slice(start, end))` for each window
2. **Rolling mode**: Maintains codon counts and recomputes incrementally (when weights are accessible)

## References

- Sharp PM, Li WH (1987) The codon adaptation index - a measure of directional synonymous codon usage bias, and its potential applications. https://doi.org/10.1093/nar/15.3.1281
- Wright F (1990) The 'effective number of codons' used in a gene. https://doi.org/10.1016/0378-1119(90)90491-9
- Ikemura T (1981) Correlation between the abundance of Escherichia coli transfer RNAs and the occurrence of the respective codons in its protein genes. https://doi.org/10.1016/0022-2836(81)90003-6

## Upstream library links

- codonbias package: https://pypi.org/project/codon-bias/
- GitHub repository: https://github.com/Benjamin-Lee/CodonBias

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
