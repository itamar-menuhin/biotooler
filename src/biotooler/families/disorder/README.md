# Disorder Family

## What this family provides

The disorder family provides protein intrinsic disorder prediction features. Intrinsically disordered regions (IDRs) are protein segments that lack stable 3D structure under physiological conditions.

**Two feature groups are available**:

### 1. Metapredict-based disorder prediction (biotooler[disorder])
- `DisorderProfileMetapredict`: Per-residue disorder probability scores (DISORDER_P)
- `DisorderDerivedScalars`: Summary statistics computed from DISORDER_P
  - DISORDER_FRAC: Fraction of disordered residues
  - DISORDER_LONGEST_IDR: Length of longest disordered region
  - DISORDER_MEAN: Mean disorder probability
  - DISORDER_P95: 95th percentile disorder probability

### 2. IDRPred consensus disorder prediction (biotooler[disorder-idrpred])
- `IDRPredConsensusMask`: Per-residue binary IDR membership mask (IDRPRED_IDR)
- `IDRPredDerivedScalars`: Summary statistics computed from IDRPRED_IDR
  - IDRPRED_FRAC_IDR: Fraction of residues in IDRs
  - IDRPRED_LONGEST_IDR_LEN: Length of longest IDR segment
  - IDRPRED_NUM_IDR_SEGMENTS: Number of distinct IDR segments

Both feature groups are **additive** and can be used together in the same analysis.

**Planned backend details**:

The disorder family is designed to support multiple disorder prediction backends. Currently implemented are metapredict (fast neural network predictor) and IDRPred (consensus predictor identifying long IDRs). Future backends may include IUPred3, ANCHOR2, and other methods.

## Intuition

Intrinsically disordered proteins (IDPs) and regions (IDRs) are functional protein segments that lack fixed 3D structure. Unlike classical structured proteins, these regions are highly flexible and dynamic, often adopting different conformations depending on their binding partners or environment.

**Key concepts**:
- **Disorder vs. Structure**: While traditional proteins fold into stable structures, disordered regions remain flexible and lack a unique native conformation
- **Functional importance**: IDRs play critical roles in signaling, regulation, and molecular recognition through conformational flexibility
- **Context-dependent behavior**: Disordered regions can undergo disorder-to-order transitions upon binding to partners

**Biological significance**:
- **Signaling and regulation**: Many signaling proteins contain IDRs that enable rapid conformational changes
- **Protein-protein interactions**: IDRs provide interaction interfaces with high specificity but low affinity
- **Post-translational modifications**: Disordered regions are often sites for phosphorylation, acetylation, and other modifications
- **Phase separation**: IDRs drive liquid-liquid phase separation forming membraneless organelles

**Applications**:
- Identifying potential regulatory regions in proteins
- Understanding protein-protein interaction mechanisms
- Designing therapeutic interventions targeting flexible regions
- Predicting protein aggregation propensity

## Mathematical formulation

### Metapredict (DisorderProfileMetapredict)

Metapredict uses a bidirectional recurrent neural network (brnn-LSTM) trained on DisProt database annotations to predict per-residue disorder scores.

**Model architecture**:
```
Input: Protein sequence (amino acids)
  ↓
Embedding layer (converts amino acids to vectors)
  ↓
Bidirectional LSTM layers (captures sequence context)
  ↓
Dense output layer
  ↓
Output: Per-residue disorder score (0.0 to 1.0)
```

**Output interpretation**:
- Disorder score: Probability that a residue is disordered (0.0 = structured, 1.0 = disordered)
- Typical threshold: 0.5 (residues with score ≥ 0.5 are predicted as disordered)
- Scores reflect local sequence composition and context

**Training data**:
- DisProt database: Manually curated disorder annotations from experimental evidence
- Positive examples: Regions confirmed disordered by NMR, X-ray crystallography, or other methods
- Negative examples: Well-structured regions from PDB structures

**Algorithm reference**: Emenecker et al. (2021) Bioinformatics. The metapredict model is described in detail in the reference below.

### IDRPred (IDRPredConsensusMask)

IDRPred is a consensus predictor that identifies long intrinsically disordered regions (IDRs) by combining multiple disorder prediction methods.

**Method**:
- Consensus approach combining multiple disorder predictors
- Focuses on identifying long, high-confidence IDRs (typically ≥30 residues)
- Returns binary mask indicating IDR membership per residue

**Output interpretation**:
- Binary mask: 0 = not in IDR, 1 = in IDR
- More conservative than individual predictors (identifies high-confidence regions)
- Segments are contiguous regions, not individual residues

**Algorithm reference**: The IDRPred consensus method is described at https://github.com/matthiasblum/idrpred

## Features and output schema

### Feature Group 1: Metapredict-based disorder prediction

Install with: `pip install "biotooler[disorder]"`

#### DisorderProfileMetapredict

Per-residue disorder prediction using metapredict neural network.

**Strict Output Contract**:
- **Feature key**: `DISORDER_P`
- **Type**: numpy array (dtype: float64)
- **Length**: L (where L = number of amino acids in protein sequence)
- **Value range**: [0.0, 1.0] (inclusive)
- **Position space**: Amino acid (residue) indexed, 0-based
- **Interpretation**: DISORDER_P[i] = probability that amino acid at position i is disordered
  - 0.0 = confidently structured
  - 1.0 = confidently disordered
  - Typical threshold: 0.5 (≥0.5 predicted as disordered)
- **Aggregation**: Mean disorder score across windows (when using windowing engine)

**Parameters**:
- `table`: NCBI genetic code table for translation (default: 1)
- `use_orf_if_present`: Use attached ORF if present (default: True)
- `strip_terminal_stop`: Strip terminal stop codon (default: True)
- `on_internal_stop`: Action for internal stops - "error" or "ignore" (default: "error")

#### DisorderDerivedScalars

Summary statistics computed from DISORDER_P without re-running predictor.

**Output**:
- `DISORDER_FRAC`: Fraction of residues with disorder probability >= threshold (scalar, float)
- `DISORDER_LONGEST_IDR`: Length of longest contiguous disordered region (scalar, int)
- `DISORDER_MEAN`: Mean disorder probability across all residues (scalar, float)
- `DISORDER_P95`: 95th percentile of disorder probabilities (scalar, float)

**Parameters**:
- `threshold`: Disorder probability threshold for determining disordered residues (default: 0.5)

**Usage pattern**: First compute DISORDER_P with DisorderProfileMetapredict, then compute derived features from the cached DISORDER_P values. This allows efficient computation of summary statistics at different thresholds without re-running the predictor.

### Feature Group 2: IDRPred consensus disorder prediction

Install with: `pip install "biotooler[disorder-idrpred]"`

Requires the `idrpred` command-line tool to be available on PATH.

#### IDRPredConsensusMask

Per-residue binary IDR membership mask using IDRPred consensus predictor.

**Strict Output Contract**:
- **Feature key**: `IDRPRED_IDR`
- **Type**: numpy array (dtype: float64)
- **Length**: L (where L = number of amino acids in protein sequence)
- **Value range**: {0.0, 1.0} (binary mask)
- **Position space**: Amino acid (residue) indexed, 0-based
- **Interpretation**: IDRPRED_IDR[i] = 1.0 if amino acid at position i is in an IDR, 0.0 otherwise
  - 0.0 = not in IDR
  - 1.0 = in IDR (high-confidence long disordered region)
- **Aggregation**: Mean across windows gives fraction of IDR residues (when using windowing engine)

**Parameters**:
- `table`: NCBI genetic code table for translation (default: 1)
- `use_orf_if_present`: Use attached ORF if present (default: True)
- `strip_terminal_stop`: Strip terminal stop codon (default: True)
- `on_internal_stop`: Action for internal stops - "error" or "ignore" (default: "error")

#### IDRPredDerivedScalars

Summary statistics computed from IDRPRED_IDR mask without re-running predictor.

**Output**:
- `IDRPRED_FRAC_IDR`: Fraction of residues in IDRs (scalar, float)
- `IDRPRED_LONGEST_IDR_LEN`: Length of longest contiguous IDR segment (scalar, int)
- `IDRPRED_NUM_IDR_SEGMENTS`: Number of distinct contiguous IDR segments (scalar, int)

**Parameters**: None (IDRPRED_IDR is already a binary mask)

**Usage pattern**: First compute IDRPRED_IDR with IDRPredConsensusMask, then compute derived features from the cached IDRPRED_IDR values.

### Additive feature groups

Both feature groups can be used together in the same analysis. Simply install both extras and use both features in your FeatureSet:

```bash
pip install "biotooler[disorder,disorder-idrpred]"
```

## Installation

### Metapredict backend

Install biotooler with disorder support:
```bash
pip install "biotooler[disorder]"
```

This installs metapredict (≥3.0) for fast neural network-based disorder prediction.

### IDRPred backend

Install biotooler with IDRPred support:
```bash
pip install "biotooler[disorder-idrpred]"
```

This installs the IDRPred package and makes the `idrpred` command-line tool available on PATH.

### Both backends

To use both feature groups in the same analysis:
```bash
pip install "biotooler[disorder,disorder-idrpred]"
```

### Checking installation

**Metapredict**:
```python
import metapredict
print(metapredict.__version__)
```

**IDRPred**:
```bash
idrpred --version
```

If metapredict or idrpred is not installed, attempting to use the respective features will raise `ImportError` with installation instructions.

## Windowing correctness

### Position space

All disorder features operate in **RESIDUE (amino acid) position space**:
- Position indices: 0-based amino acid positions
- Example: Sequence "MKTAY" has positions 0=M, 1=K, 2=T, 3=A, 4=Y
- Both `DisorderProfileMetapredict` and `IDRPredConsensusMask` return `PositionSpace.RESIDUE`

### Vector computation

Disorder features compute **per-residue vectors across the entire sequence** using upstream library APIs:

1. **Full-sequence input**: Pass complete protein sequence to predictor
2. **Per-residue values**: Models return probability or mask for each amino acid position
3. **Vector output**: numpy array with length equal to sequence length
4. **Window aggregation**: Windowing engine slices the full-context vector and applies aggregation

**Why full-sequence?** Context matters for disorder prediction. Neural networks and consensus predictors consider surrounding amino acids. Computing on substrings would lose context and produce incorrect results.

**Upstream API usage**:
- **Metapredict**: Calls `metapredict.predict_disorder(protein_seq)` which returns a per-residue numpy array of disorder probabilities
- **IDRPred**: Calls `idrpred` CLI tool via subprocess with FASTA input, parses TSV output to generate per-residue binary mask

Both upstream libraries are documented to require full sequences for accurate predictions. See the [metapredict documentation](https://metapredict.readthedocs.io/) and [IDRPred GitHub](https://github.com/matthiasblum/idrpred) for details on their prediction algorithms.

### Aggregation strategy

When using sliding windows with `FeatureSet.compute_windows()`, per-residue vectors are aggregated:

**For DISORDER_P** (metapredict):
- Aggregation: `np.mean(DISORDER_P[window_start:window_end])`
- Interpretation: Average disorder propensity across window region
- Wide-format columns: `DISORDER_P_0`, `DISORDER_P_50`, `DISORDER_P_100`, etc.

**For IDRPRED_IDR** (IDRPred):
- Aggregation: `np.mean(IDRPRED_IDR[window_start:window_end])`
- Interpretation: Fraction of residues in IDRs within window
- Wide-format columns: `IDRPRED_IDR_0`, `IDRPRED_IDR_50`, `IDRPRED_IDR_100`, etc.

**Example window naming**:
```python
# 3 windows at positions 0, 50, 100
# Output columns: "disorder.DISORDER_P_0", "disorder.DISORDER_P_50", "disorder.DISORDER_P_100"
# Each contains the mean disorder score for that window
```

### Testing approach

Windowing correctness for this family focuses on:
1. **Full-sequence computation**: Verify that disorder predictors (metapredict and IDRPred) receive complete protein sequences, not substrings, to preserve context
2. **Upstream library correctness**: 
   - For metapredict: Results should match direct `metapredict.predict_disorder()` calls on the same sequence
   - For IDRPred: Binary masks should match direct `idrpred` CLI output for the same sequence
3. **Translation correctness**: For DNA/RNA inputs, ensure proper translation to protein before disorder prediction
4. **ORF handling**: When ORFs are present, verify correct extraction and translation
5. **Position space mapping**: Vector indices correctly map to amino acid positions (0-based)
6. **Aggregation accuracy**: Mean aggregation over windows produces expected values
7. **Determinism**: Same input sequence produces identical predictions across multiple runs
8. **Edge cases**: Handle empty sequences, very short sequences, and terminal/internal stop codons

Tests validate that:
- `compute_vector()` is called once per sequence with the full sequence
- Vector lengths match protein sequence lengths
- Windowing engine correctly slices full vectors without recomputation
- Wide-format output columns follow naming convention (feature_position pattern)

See `tests/families/disorder_family/` for comprehensive test coverage.


## Examples

### Computing both metapredict and IDRPred features together

This example shows how to compute both DISORDER_P and IDRPRED_IDR in a single run:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict
from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask
from biotooler.features.sets import FeatureSet

# Create both features
metapredict_feature = DisorderProfileMetapredict()
idrpred_feature = IDRPredConsensusMask()

# Create a protein sequence record
protein = SeqRecord(Seq("MKALVSWGRPQMTESTDISORDERREGION"), id="protein1")
protein.annotations["molecule_type"] = "protein"

# Compute both disorder predictions
feature_set = FeatureSet({
    "metapredict": metapredict_feature,
    "idrpred": idrpred_feature,
})

# Compute in sliding windows (wide-format output)
results = feature_set.compute_windows(
    protein,
    window_size=10,  # 10 amino acid window
    step=5           # Step by 5 amino acids
)

# Results contain both predictions with columns:
# "metapredict.DISORDER_P_0", "metapredict.DISORDER_P_5", "metapredict.DISORDER_P_10", ...
# "idrpred.IDRPRED_IDR_0", "idrpred.IDRPRED_IDR_5", "idrpred.IDRPRED_IDR_10", ...
print(results.columns.tolist())
print(results[["metapredict.DISORDER_P_0", "idrpred.IDRPRED_IDR_0"]])
```

### Basic usage: Computing metapredict disorder predictions

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict
import numpy as np

# Create feature
feature = DisorderProfileMetapredict()

# Create protein sequence record
record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="protein1")
record.annotations["molecule_type"] = "protein"

# Compute disorder scores (one score per residue)
result = feature.compute_vector(record)
disorder_scores = result["DISORDER_P"]  # numpy array

print(f"Disorder scores: {disorder_scores}")
print(f"Length: {len(disorder_scores)}")  # 16 (number of residues)
print(f"Mean disorder: {disorder_scores.mean():.3f}")
print(f"Max disorder: {disorder_scores.max():.3f}")

# Identify disordered residues (threshold = 0.5)
disordered_residues = disorder_scores >= 0.5
print(f"Disordered positions: {np.where(disordered_residues)[0]}")
```

### Computing IDRPred consensus mask

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.disorder.idrpred_feature import IDRPredConsensusMask
import numpy as np

# Create feature
feature = IDRPredConsensusMask()

# Create protein sequence record
record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="protein1")
record.annotations["molecule_type"] = "protein"

# Compute IDR mask (binary mask per residue)
result = feature.compute_vector(record)
idr_mask = result["IDRPRED_IDR"]  # numpy array of 0s and 1s

print(f"IDR mask: {idr_mask}")
print(f"Fraction in IDRs: {idr_mask.mean():.2%}")
print(f"IDR positions: {np.where(idr_mask == 1.0)[0]}")
```

### Computing derived features without re-running predictor

```python
import numpy as np
from biotooler.families.disorder.derived import DisorderDerivedScalars
from biotooler.families.disorder.idrpred_derived import IDRPredDerivedScalars

# Assume DISORDER_P has already been computed and stored in record.annotations
record.annotations["DISORDER_P"] = result["DISORDER_P"]

# Compute metapredict-based summary statistics
disorder_derived = DisorderDerivedScalars(threshold=0.5)
disorder_summary = disorder_derived(record)

print(f"Metapredict summaries:")
print(f"  Fraction disordered: {disorder_summary['DISORDER_FRAC']:.2%}")
print(f"  Longest IDR: {disorder_summary['DISORDER_LONGEST_IDR']} residues")
print(f"  Mean disorder: {disorder_summary['DISORDER_MEAN']:.3f}")
print(f"  95th percentile: {disorder_summary['DISORDER_P95']:.3f}")

# Assume IDRPRED_IDR has also been computed
record.annotations["IDRPRED_IDR"] = result["IDRPRED_IDR"]

# Compute IDRPred-based summary statistics
idrpred_derived = IDRPredDerivedScalars()
idrpred_summary = idrpred_derived(record)

print(f"\nIDRPred summaries:")
print(f"  Fraction in IDRs: {idrpred_summary['IDRPRED_FRAC_IDR']:.2%}")
print(f"  Longest IDR: {idrpred_summary['IDRPRED_LONGEST_IDR_LEN']} residues")
print(f"  Number of IDR segments: {idrpred_summary['IDRPRED_NUM_IDR_SEGMENTS']}")
```

### Windowed analysis with FeatureSet

```python
from biotooler.features.sets import FeatureSet

# Create feature set with both disorder predictors
metapredict_feature = DisorderProfileMetapredict()
idrpred_feature = IDRPredConsensusMask()

feature_set = FeatureSet({
    "metapredict": metapredict_feature,
    "idrpred": idrpred_feature,
})

# Compute disorder in sliding windows (50 AA window, 25 AA step)
protein = SeqRecord(Seq("M" * 200), id="long_protein")  # 200 AA protein
protein.annotations["molecule_type"] = "protein"

results = feature_set.compute_windows(
    protein,
    window_size=50,
    step=25
)

# Wide-format output with columns like:
# "metapredict.DISORDER_P_0", "metapredict.DISORDER_P_25", "metapredict.DISORDER_P_50", ...
# "idrpred.IDRPRED_IDR_0", "idrpred.IDRPRED_IDR_25", ...

print(results.columns)
print(results[["metapredict.DISORDER_P_0", "metapredict.DISORDER_P_25",
               "idrpred.IDRPRED_IDR_0", "idrpred.IDRPRED_IDR_25"]])
```

## Edge cases and validation

The disorder family handles common edge cases:

1. **Empty sequences**: Returns empty arrays or zero-filled summaries
2. **DNA/RNA sequences**: Automatically translated via `ensure_protein_record()`
3. **Stop codons**: Terminal stops stripped by default; internal stops raise error
4. **Non-standard amino acids**: May raise `ValueError` if not supported by predictor
5. **Derived features**: Require corresponding vector (DISORDER_P or IDRPRED_IDR) in `record.annotations`

**Known limitations**:
- Very short sequences (< 10 residues): Limited context for prediction
- Transmembrane regions: May be incorrectly predicted as disordered (metapredict trained on soluble proteins)
- Low-complexity regions: Often predicted as disordered (biologically accurate in many cases)
- Post-translational modifications: Not accounted for (predictions based on primary sequence only)

## Maintenance notes

### Dependencies

**Required**:
- Python ≥ 3.11
- `biopython` - For SeqRecord and Seq objects
- `numpy` - For vector computations and aggregation

**Optional (metapredict group)**:
- `metapredict` ≥ 3.0 (required to use DisorderProfileMetapredict and DisorderDerivedScalars)
  - Install via: `pip install "biotooler[disorder]"`
  - Not included in default biotooler installation to keep dependencies lightweight

**Optional (IDRPred group)**:
- `idrpred` package (required to use IDRPredConsensusMask and IDRPredDerivedScalars)
  - Install via: `pip install "biotooler[disorder-idrpred]"`
  - Provides `idrpred` command-line tool on PATH
  - Repository: https://github.com/matthiasblum/idrpred

### Implementation status

**Implemented**:
- ✅ Family infrastructure (registry, imports, lazy loading)
- ✅ Integration helpers (require_metapredict(), require_idrpred_cli())
- ✅ DisorderProfileMetapredict feature (per-residue predictions)
- ✅ DisorderDerivedScalars feature (summary statistics from DISORDER_P)
- ✅ IDRPredConsensusMask feature (per-residue binary mask)
- ✅ IDRPredDerivedScalars feature (summary statistics from IDRPRED_IDR)
- ✅ Comprehensive test coverage
- ✅ Windowing engine integration

**Planned enhancements**:
- Support for additional backends (IUPred3, ANCHOR2, etc.)
- Batch prediction optimization for multiple sequences
- Additional derived features (disorder clusters, region boundaries)

### Design decisions

**1. Two feature groups (metapredict and IDRPred)**:
   - Rationale: Complementary approaches (neural network vs consensus predictor)
   - metapredict: Fast, continuous probability scores
   - IDRPred: Conservative, identifies high-confidence long IDRs
   - Both operate on same position space (residues) and can be used together

**2. Lazy imports**:
   - Rationale: Deep learning dependencies (PyTorch for metapredict) kept optional
   - Implementation: `require_metapredict()` and `require_idrpred_cli()` check on first use
   - Benefits: Fast biotooler import, lightweight default installation

**3. Full-sequence computation**:
   - Rationale: Disorder predictors use sequence context; substrings would reduce accuracy
   - Trade-off: Single computation per sequence regardless of number of windows
   - Windowing engine slices full vector after prediction

**4. Vector-based output**:
   - Rationale: Per-residue values provide maximum information
   - Users can aggregate using mean, max, or custom functions via windowing engine
   - Derived features compute summaries from vectors without re-running predictors

### Performance considerations

**Metapredict performance**:
- Time complexity: O(n) for sequence of length n (LSTM model)
- Memory: Scales with sequence length and model size (~100 MB for model)
- Typical timing: Fast inference (~1-10ms per sequence of length 100-1000)
- GPU acceleration: metapredict uses PyTorch; GPU can speed up batch predictions
- Caching: Results cached in `record.annotations["DISORDER_P"]` for derived features

**IDRPred performance**:
- Time complexity: O(n) for sequence of length n (subprocess overhead + prediction)
- Typical timing: Slower than metapredict (~100-1000ms per sequence)
- No GPU acceleration (external CLI tool)
- Caching: Results cached in `record.annotations["IDRPRED_IDR"]` for derived features

**DisorderDerivedScalars / IDRPredDerivedScalars**:
- Time complexity: O(n) for sequence of length n
- Memory: Minimal (no model loading, operates on cached vectors)
- Typical timing: Microseconds per sequence
- Caching benefit: Does not re-run predictor, reuses vector from annotations

## References

- **metapredict**: Emenecker, R.J., Griffith, D. & Holehouse, A.S. (2021) "Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure" *Bioinformatics* 37(26):5035-5037. https://doi.org/10.1093/bioinformatics/btab527
- **IDRPred**: Blum, M., et al. "IDRPred: A consensus-based disorder prediction method" https://github.com/matthiasblum/idrpred
- **DisProt database**: Quaglia, F., et al. (2022) "DisProt in 2022: improved quality and accessibility of protein intrinsic disorder annotation" *Nucleic Acids Research* 50:D480-D487. https://doi.org/10.1093/nar/gkab1082
- **Intrinsically disordered proteins review**: van der Lee, R., et al. (2014) "Classification of intrinsically disordered regions and proteins" *Chemical Reviews* 114(13):6589-6631. https://doi.org/10.1021/cr400525m

## Upstream library links

- **metapredict GitHub**: https://github.com/idptools/metapredict
- **metapredict Documentation**: https://metapredict.readthedocs.io/
- **metapredict PyPI**: https://pypi.org/project/metapredict/
- **IDRPred GitHub**: https://github.com/matthiasblum/idrpred
- **IDPTools Suite**: https://idptools-parrot.readthedocs.io/
- **DisProt Database**: https://disprot.org/
