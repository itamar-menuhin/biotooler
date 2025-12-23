# Disorder Family

## What this family provides

The disorder family provides protein intrinsic disorder prediction features using metapredict as the backend. Intrinsically disordered regions (IDRs) are protein segments that lack stable 3D structure under physiological conditions.

**Available Features**:
- `DisorderProfileMetapredict`: Per-residue disorder prediction scores (DISORDER_P)
- `DisorderDerivedScalars`: Summary statistics computed from DISORDER_P
  - DISORDER_FRAC: Fraction of disordered residues
  - DISORDER_LONGEST_IDR: Length of longest disordered region
  - DISORDER_MEAN: Mean disorder probability
  - DISORDER_P95: 95th percentile disorder probability

**Future Backends**:
- IUPred3 backend (planned)
- Consensus predictions from multiple methods (planned)

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

## Features and output schema

**Available Features**:

### 1. DisorderProfileMetapredict

Per-residue disorder prediction using metapredict neural network.

**Output**:
- `DISORDER_P`: Vector of disorder probabilities (0.0 to 1.0) for each amino acid position
- Aggregation: Mean disorder score across windows

**Parameters**:
- `table`: NCBI genetic code table for translation (default: 1)
- `use_orf_if_present`: Use attached ORF if present (default: True)
- `strip_terminal_stop`: Strip terminal stop codon (default: True)
- `on_internal_stop`: Action for internal stops - "error" or "ignore" (default: "error")

### 2. DisorderDerivedScalars

Summary statistics computed from DISORDER_P without re-running predictor.

**Output**:
- `DISORDER_FRAC`: Fraction of residues with disorder probability >= threshold
- `DISORDER_LONGEST_IDR`: Length of longest contiguous disordered region
- `DISORDER_MEAN`: Mean disorder probability across all residues
- `DISORDER_P95`: 95th percentile of disorder probabilities

**Parameters**:
- `threshold`: Disorder probability threshold for determining disordered residues (default: 0.5)

**Usage pattern**: First compute DISORDER_P with DisorderProfileMetapredict, then compute derived features from the cached DISORDER_P values. This allows efficient computation of summary statistics at different thresholds without re-running the predictor.

**Integration**: Features integrate with biotooler's windowing engine for sliding window analysis across protein sequences.

## References

- **metapredict**: Emenecker, R.J., Griffith, D. & Holehouse, A.S. (2021) "Metapredict: a fast, accurate, and easy-to-use predictor of consensus disorder and structure" *Bioinformatics* 37(26):5035-5037. https://doi.org/10.1093/bioinformatics/btab527
- **DisProt database**: Quaglia, F., et al. (2022) "DisProt in 2022: improved quality and accessibility of protein intrinsic disorder annotation" *Nucleic Acids Research* 50:D480-D487. https://doi.org/10.1093/nar/gkab1082
- **Intrinsically disordered proteins review**: van der Lee, R., et al. (2014) "Classification of intrinsically disordered regions and proteins" *Chemical Reviews* 114(13):6589-6631. https://doi.org/10.1021/cr400525m
- **IUPred3**: Erdős, G., et al. (2021) "IUPred3: prediction of protein disorder enhanced with unambiguous experimental annotation and visualization of evolutionary conservation" *Nucleic Acids Research* 49:W297-W303. https://doi.org/10.1093/nar/gkab408

## Upstream library links

- **metapredict GitHub**: https://github.com/idptools/metapredict
- **metapredict Documentation**: https://metapredict.readthedocs.io/
- **metapredict PyPI**: https://pypi.org/project/metapredict/
- **IDPTools Suite**: https://idptools-parrot.readthedocs.io/
- **DisProt Database**: https://disprot.org/

## Examples

### Basic usage: Computing disorder predictions

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict

# Create feature
feature = DisorderProfileMetapredict()

# Create protein sequence record
record = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="protein1")
record.annotations["molecule_type"] = "protein"

# Compute disorder scores (one score per residue)
result = feature.compute_vector(record)
disorder_scores = result["DISORDER_P"]  # numpy array

print(f"Disorder scores: {disorder_scores}")
print(f"Mean disorder: {disorder_scores.mean():.3f}")
```

### Computing derived features without re-running predictor

```python
from biotooler.families.disorder.derived import DisorderDerivedScalars

# Assume DISORDER_P has already been computed and stored in record.annotations
# (by DisorderProfileMetapredict or loaded from cache)

# Compute summary statistics
derived = DisorderDerivedScalars(threshold=0.5)
summary = derived(record)

print(f"Fraction disordered: {summary['DISORDER_FRAC']:.2%}")
print(f"Longest IDR: {summary['DISORDER_LONGEST_IDR']} residues")
print(f"Mean disorder: {summary['DISORDER_MEAN']:.3f}")
print(f"95th percentile: {summary['DISORDER_P95']:.3f}")
```

### Using different thresholds efficiently

```python
# Compute derived features at multiple thresholds
# without re-running metapredict
thresholds = [0.5, 0.6, 0.7]
for t in thresholds:
    derived = DisorderDerivedScalars(threshold=t)
    result = derived(record)
    print(f"Threshold {t}: {result['DISORDER_FRAC']:.2%} disordered")
```

## Edge cases and validation

**Status**: Edge case handling will be documented once features are implemented.

**Planned validation**:
- **Input validation**: Verify protein sequences (20 standard amino acids)
- **Sequence length**: Handle short sequences (minimum length requirements TBD)
- **Invalid characters**: Reject or handle non-standard amino acids
- **Empty sequences**: Return appropriate empty results

**Expected edge cases**:
- Very short sequences (< 10 residues): Limited context for disorder prediction
- Transmembrane regions: May be incorrectly predicted as disordered due to sequence composition
- Low-complexity regions: Often predicted as disordered (poly-Q, poly-A stretches)
- Modified residues: Standard metapredict does not account for post-translational modifications

## Windowing correctness

This section will describe windowing behavior once features are implemented.

### Position space

**Planned**: Features will operate at the **RESIDUE position space** (per-amino-acid level).
- Position indices: 0-based amino acid positions
- Example: Sequence "MKTAY" has positions 0=M, 1=K, 2=T, 3=A, 4=Y

### Vector computation

**Planned approach**: Per-residue disorder prediction

Disorder prediction will be computed as a **vector across the entire sequence**:

1. **Full-sequence input**: Pass complete protein sequence to metapredict
2. **Per-residue scores**: Model returns disorder probability for each amino acid position
3. **Vector output**: numpy array with length equal to sequence length
4. **Window aggregation**: When windowing is applied, disorder scores will be aggregated using mean or other appropriate functions

**Why full-sequence?** Context matters for disorder prediction. The neural network considers surrounding amino acids when predicting disorder at each position. Computing on substrings would lose important sequence context.

**Upstream API** (planned):
```python
# Simplified planned implementation
import metapredict
scores = metapredict.predict_disorder(protein_sequence)
# Returns array of disorder probabilities, one per residue
```

### Aggregation strategy

**Planned aggregation functions**:

**Mean disorder score** (primary):
- Aggregation: `np.mean(disorder_scores[window_start:window_end])`
- Interpretation: Average disorder propensity across a window region
- Use case: Identify regions with overall high or low disorder

**Max disorder score** (alternative):
- Aggregation: `np.max(disorder_scores[window_start:window_end])`  
- Interpretation: Peak disorder in window
- Use case: Detect presence of any highly disordered residue

**Fraction disordered** (alternative):
- Aggregation: `np.mean(disorder_scores[window_start:window_end] >= threshold)`
- Interpretation: Proportion of residues predicted as disordered
- Use case: Binary classification of windows

### Testing approach

**Planned tests** (once features are implemented):

1. **Upstream library correctness**:
   - Compare feature output to direct metapredict calls
   - Verify scores match for same input sequence
   - Test determinism (same input → same output)

2. **Position space correctness**:
   - Verify vector lengths match sequence lengths
   - Confirm window positions map correctly to residues
   - Test edge cases at sequence boundaries

3. **Full-context computation**:
   - Verify full sequence is used for prediction (not substrings)
   - Test that context is maintained for accurate predictions

4. **Edge case handling**:
   - Test empty sequences, very short sequences
   - Validate handling of non-standard amino acids
   - Check behavior with unusual sequence compositions

5. **Windowing engine integration**:
   - Verify compute_vector() is called once per sequence
   - Test aggregation functions work correctly
   - Validate window boundaries and partial windows

## Maintenance notes

### Dependencies

**Required**:
- Python ≥ 3.11
- `biopython` - For SeqRecord and Seq objects
- `numpy` - For vector computations and aggregation

**Optional**:
- `metapredict` ≥ 3.0 (required to use this family once features are implemented)
  - Install via: `pip install metapredict` or `pip install "biotooler[disorder]"`
  - Not included in default biotooler installation to keep dependencies lightweight

### Installation

**Install biotooler with disorder support**:
```bash
pip install "biotooler[disorder]"
```

**Install metapredict separately**:
```bash
pip install metapredict
```

**Check installation**:
```python
import metapredict  # Should not raise ImportError
print(metapredict.__version__)
```

If metapredict is not installed, attempting to use disorder features will raise `ImportError` with installation instructions.

### Implementation status

- ✅ Family infrastructure (registry, imports, lazy loading)
- ✅ Integration helper (`require_metapredict()`)
- ✅ DisorderProfileMetapredict feature (per-residue predictions)
- ✅ DisorderDerivedScalars feature (summary statistics)
- ✅ Comprehensive test coverage (35+ tests)
- ✅ Windowing engine integration

### Design decisions

**1. Metapredict as initial backend**:
   - Rationale: Fast, accurate, well-maintained, and easy to install
   - Alternative backends (IUPred3, consensus) planned for future releases

**2. Lazy imports**:
   - Rationale: metapredict has deep learning dependencies (PyTorch); lazy loading keeps biotooler import fast
   - Implementation: `require_metapredict()` helper checks and imports metapredict on first use

**3. Full-sequence computation** (planned):
   - Rationale: Disorder prediction models use sequence context; computing on substrings would reduce accuracy
   - Trade-off: Single computation per sequence regardless of number of windows

**4. Vector-based output** (planned):
   - Rationale: Per-residue disorder scores provide maximum information
   - Users can aggregate using mean, max, or custom functions via windowing engine

### Performance considerations

**Metapredict performance** (once implemented):
- Time complexity: O(n) for sequence of length n (LSTM model)
- Memory: Scales with sequence length and model size
- Typical timing: Fast inference (~1-10ms per sequence)
- GPU acceleration: metapredict uses PyTorch; GPU can speed up batch predictions

### Future work

**Planned enhancements**:
- Support for IUPred3 backend
- Consensus predictions combining multiple methods
- Integration with structure prediction features
- Batch prediction optimization for multiple sequences
- Additional derived features (disorder clusters, region boundaries)

### Troubleshooting

**Import Error: "No module named 'metapredict'"**
- Solution: Install metapredict with `pip install metapredict` or `pip install "biotooler[disorder]"`

**Feature not available error**:
- Current status: No features implemented yet in disorder family
- Features are coming in future updates
- The `require_metapredict()` helper is available for testing integration

### Version compatibility

- metapredict ≥ 3.0 recommended (latest stable version)
- biotooler ≥ 0.1.0 required for lazy import infrastructure
- Future feature implementations will specify exact version requirements
