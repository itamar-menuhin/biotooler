# Codon Bias Features

The `CodonBiasFeature` class provides a unified interface for computing codon usage bias metrics using the external [codon-bias](https://pypi.org/project/codon-bias/) package. This feature wraps codonbias models to compute various codon usage statistics on DNA/RNA sequences.

## Overview

Codon usage bias refers to the non-uniform usage of synonymous codons in protein-coding sequences. Different organisms and genes exhibit distinct patterns of codon preference. The biotooler library provides a feature class that wraps the `codonbias` package to compute various bias metrics including:

- **CAI** (Codon Adaptation Index): Measures how similar a sequence's codon usage is to a reference set
- **ENC** (Effective Number of Codons): Quantifies the degree of codon bias
- **FOP** (Frequency of Optimal Codons): Fraction of optimal codons in a sequence
- **RSCU** (Relative Synonymous Codon Usage): Relative usage of synonymous codons
- **RCBS/DCBS** (Relative Codon Bias Score / Distance from Codon Bias Score)
- **tAI** (tRNA Adaptation Index): Measures translation efficiency based on tRNA availability
- **nTE** (Normalized Translation Efficiency)
- **CPB/CPS** (Codon Pair Bias/Score)

## Basic Usage

### Computing Scalar Features

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from codonbias.scores import CodonAdaptationIndex, EffectiveNumberOfCodons
from biotooler.features import CodonBiasFeature

# Create codonbias models
# CAI requires a reference sequence representing highly expressed genes
ref_seq = "ATGATGATGATGATGATGATGATGATG"
cai = CodonAdaptationIndex(ref_seq)

# ENC doesn't require a reference
enc = EffectiveNumberOfCodons()

# Create the feature with custom names
feature = CodonBiasFeature([cai, enc], names=["CAI", "ENC"])

# Compute features on a sequence
test_seq = "ATGATGATGATGATGATG"
record = SeqRecord(Seq(test_seq), id="test_gene")
result = feature(record)

print(result)
# {'CAI': 1.0, 'ENC': 61.0}
```

### Using Default Names

If you don't provide custom names, the feature uses the model class names:

```python
feature = CodonBiasFeature([cai, enc])
result = feature(record)

print(result.keys())
# dict_keys(['CodonAdaptationIndex', 'EffectiveNumberOfCodons'])
```

## Windowed Analysis

The `CodonBiasFeature` supports both baseline and incremental window computation modes when used with `FeatureSet.compute_orf_windows()`.

### Window Naming Convention

Windows are named with suffixes indicating the absolute nucleotide start position:
- `CAI_0`: CAI for window starting at position 0
- `CAI_3`: CAI for window starting at position 3
- `ENC_6`: ENC for window starting at position 6

### Baseline Window Mode

In baseline mode, each window is computed independently using `model.get_score(seq, slice=slice(start, end))`:

```python
from biotooler.features import FeatureSet

# Create models
ref_seq = "ATGATGATGATGATGATGATGATG"
cai = CodonAdaptationIndex(ref_seq)
enc = EffectiveNumberOfCodons()

# Create feature and feature set
feature = CodonBiasFeature([cai, enc], names=["CAI", "ENC"])
fs = FeatureSet(feature, name="codon_bias")

# Compute windowed features
test_seq = "ATGATGATGATGATGATGATGATGATGATG"  # 10 codons = 30 nt
record = SeqRecord(Seq(test_seq), id="test")

result = fs.compute_orf_windows(
    record,
    orf=(0, 30),
    window_nt=9,   # 3 codons per window
    step_nt=3      # step by 1 codon (must be multiple of 3)
)

print(result.columns.tolist())
# ['record_id', 'orf_start', 'orf_end', 
#  'codon_bias.CAI_0', 'codon_bias.CAI_3', 'codon_bias.CAI_6', ...,
#  'codon_bias.ENC_0', 'codon_bias.ENC_3', 'codon_bias.ENC_6', ...]
```

### Rolling (Incremental) Window Mode

For models with accessible weights (like CAI which has `log_weights`), the feature automatically uses incremental computation:

1. **Initialization**: Uses `codonbias.stats.CodonCounter` to count codons in the first window
2. **Stepping**: Updates codon counts by subtracting outgoing codons and adding incoming codons
3. **Scoring**: Computes scores using `codonbias.utils.geomean()` (for log weights) or `codonbias.utils.mean()` (for linear weights)

This provides significant performance improvements for large sequences with many overlapping windows.

```python
# The feature automatically detects if incremental mode is possible
# For CAI (has log_weights): uses incremental computation
# For ENC (no weights): falls back to baseline per-window calls

# Both produce the same results (within floating-point tolerance)
result_incremental = fs.compute_orf_windows(
    record, orf=(0, 30), window_nt=9, step_nt=3
)
```

### Weight Detection and Fallback

The feature tries to extract weights in this order:

1. `model.log_weights` (preferred) → uses `codonbias.utils.geomean()`
2. `model.weights` → uses `codonbias.utils.mean()`
3. `model.get_weights()` if callable → uses `codonbias.utils.mean()`
4. No accessible weights → falls back to baseline `get_score()` calls

Models without accessible weights still work but use baseline mode for each window.

## Multiple Models

You can compute multiple codon bias metrics simultaneously:

```python
from codonbias.scores import (
    CodonAdaptationIndex,
    EffectiveNumberOfCodons,
    FrequencyOfOptimalCodons,
)

# Create multiple models
ref_seq = "ATGATGATGATGATGATGATGATGATG"
cai = CodonAdaptationIndex(ref_seq)
enc = EffectiveNumberOfCodons()
fop = FrequencyOfOptimalCodons(ref_seq)

# Combine in single feature
feature = CodonBiasFeature(
    [cai, enc, fop],
    names=["CAI", "ENC", "FOP"]
)
fs = FeatureSet(feature, name="cb")

# All metrics computed together
result = fs.compute_orf_windows(
    record, orf=(0, 30), window_nt=12, step_nt=6
)

print([c for c in result.columns if c.startswith("cb.")])
# ['cb.CAI_0', 'cb.CAI_6', 'cb.CAI_12', ...,
#  'cb.ENC_0', 'cb.ENC_6', 'cb.ENC_12', ...,
#  'cb.FOP_0', 'cb.FOP_6', 'cb.FOP_12', ...]
```

## Validation and Special Handling

### Alphabet Validation

Codon bias features only work with DNA or RNA sequences:

```python
# This will raise ValueError
protein_record = SeqRecord(Seq("MKALVSWGR"), id="protein")
protein_record.annotations["molecule_type"] = "protein"

try:
    feature(protein_record)
except ValueError as e:
    print(e)
    # "Codon bias features apply only to DNA/RNA sequences, not protein sequences"
```

### RNA to DNA Conversion

RNA sequences (containing U/u) are automatically converted to DNA (T/t) before calling codonbias:

```python
# RNA sequence with U
rna_record = SeqRecord(Seq("AUGUGAUGUGAUG"), id="rna")
result_rna = feature(rna_record)

# DNA sequence with T
dna_record = SeqRecord(Seq("ATGATGATGATG"), id="dna")
result_dna = feature(dna_record)

# Results are identical after conversion
assert result_rna["CAI"] == result_dna["CAI"]
```

### Codon Window Requirements

When using `compute_orf_windows`, both `window_nt` and `step_nt` must be multiples of 3 (enforced by `iter_orf_codon_windows`):

```python
# This will raise ValueError
try:
    fs.compute_orf_windows(
        record, orf=(0, 30),
        window_nt=9,
        step_nt=2  # Not a multiple of 3!
    )
except ValueError as e:
    print(e)
    # "step_nt must be a multiple of 3 for codon-aligned windows"

# This works
result = fs.compute_orf_windows(
    record, orf=(0, 30),
    window_nt=9,
    step_nt=3  # Multiple of 3
)
```

## Performance Considerations

### When to Use Incremental Mode

Incremental mode provides the most benefit when:
- You have many overlapping windows (small `step_nt` relative to `window_nt`)
- You're using models with accessible weights (e.g., CAI, FOP)
- You're analyzing long sequences

### Baseline Mode Performance

Even when incremental mode is not available (models without weights), the feature still benefits from:
- Efficient sequence caching in biotooler
- Direct slice-based calls to codonbias models
- Minimal overhead for window management

## Examples with Different Models

### Using CAI with a Reference Set

```python
# CAI requires a reference sequence representing highly expressed genes
# In practice, you might use actual high-expression genes from your organism
highly_expressed_genes = "ATGATGATGATGATGATGATGATGATGATG"
cai = CodonAdaptationIndex(highly_expressed_genes)

feature = CodonBiasFeature([cai], names=["CAI"])
result = feature(record)
```

### Using ENC (No Reference Required)

```python
# ENC measures codon bias without needing a reference
enc = EffectiveNumberOfCodons()
feature = CodonBiasFeature([enc], names=["ENC"])

# ENC values typically range from 20 (extreme bias) to 61 (no bias)
result = feature(record)
print(f"ENC: {result['ENC']:.2f}")
```

### Using FOP with Optimal Codons

```python
# FOP requires a reference to determine optimal codons
ref_seq = "ATGATGATGATGATGATGATG"
fop = FrequencyOfOptimalCodons(ref_seq)

feature = CodonBiasFeature([fop], names=["FOP"])
result = feature(record)
print(f"FOP: {result['FOP']:.2%}")  # Display as percentage
```

## See Also

- [codon-bias package documentation](https://pypi.org/project/codon-bias/)
- [ORF Windows](./orf_windows.md) for general windowing concepts
- [Feature Sets](../README.md#feature-sets) for feature computation patterns
