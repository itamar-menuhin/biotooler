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
- IUPred3 backend (external binary/script)
- ANCHOR backend (context-dependent binding regions)
- Consensus predictions from multiple methods

**Planned backend details**:

### IUPred3
- **Type**: External command-line tool
- **Installation**: Download standalone binary/script from https://iupred3.elte.hu/
- **Interface**: Python subprocess wrapper calling the IUPred3 executable
- **Output**: Per-residue disorder scores (similar to metapredict)
- **Features**: Multiple disorder types (long, short, structured domains)
- **Implementation approach**:
  - Feature class: `DisorderProfileIUPred3`
  - Lazy import check for IUPred3 installation
  - Subprocess wrapper: `subprocess.run(["iupred3", "long", input_fasta])`
  - Parse output format (FASTA or tab-delimited)
  - Return DISORDER_P_IUPRED3 vector

### ANCHOR2
- **Type**: External command-line tool (part of IUPred3 package)
- **Purpose**: Predict context-dependent protein binding regions within disordered segments
- **Installation**: Bundled with IUPred3
- **Interface**: Similar subprocess wrapper pattern
- **Output**: Per-residue binding region scores
- **Features**: Identifies regions that can gain structure upon binding
- **Implementation approach**:
  - Feature class: `DisorderProfileAnchor`
  - Return ANCHOR_P vector (binding probability per residue)
  - Can be combined with disorder predictions for functional annotation

### Consensus Backend
- **Type**: Meta-predictor combining multiple backends
- **Purpose**: Improve accuracy by ensembling metapredict, IUPred3, and others
- **Implementation approach**:
  - Feature class: `DisorderProfileConsensus`
  - Requires multiple backends installed
  - Aggregation strategies:
    - Mean: `(metapredict + iupred3) / 2`
    - Voting: Threshold each, take majority
    - Weighted mean: Learned weights from benchmark datasets
  - Return DISORDER_P_CONSENSUS vector

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

### DNA/RNA Translation

Both `DisorderProfileMetapredict` and derived features accept DNA/RNA sequences as input. Translation to protein is handled automatically via the shared helper `biotooler.core.translation.ensure_protein_record()`:

**Translation behavior**:
1. **Molecule type detection**: Uses `record.annotations["molecule_type"]` (DNA/RNA/protein)
2. **ORF handling**: 
   - If `use_orf_if_present=True` (default): Uses attached ORF from `record.annotations["biotooler.orf"]`
   - If no attached ORF: Uses full sequence in frame 0
   - Manual ORF can be specified via `orf=(start, end)` parameter
3. **Genetic code**: Configurable via `table` parameter (default: 1 = standard code)
4. **Stop codon handling**:
   - Terminal stops: Stripped by default (`strip_terminal_stop=True`)
   - Internal stops: Raises error by default (`on_internal_stop="error"`)
5. **Output**: Protein SeqRecord with translated sequence and preserved metadata

**Example translation workflow**:
```python
from biotooler.core.translation import ensure_protein_record

# DNA record with ORF
dna_record = SeqRecord(Seq("ATGAAACGCTTA"), id="gene1")
dna_record.annotations["molecule_type"] = "DNA"
dna_record.annotations["biotooler.orf"] = (0, 12)  # Attached ORF

# Automatic translation (used internally by DisorderProfileMetapredict)
protein_record = ensure_protein_record(dna_record, table=1)
# protein_record.seq = "MKR" (translated from ORF)
```

This translation layer is shared across all protein-level feature families (disorder, protparam, etc.), ensuring consistent behavior.

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
print(f"Length: {len(disorder_scores)}")  # 16 (number of residues)
print(f"Mean disorder: {disorder_scores.mean():.3f}")
print(f"Max disorder: {disorder_scores.max():.3f}")

# Identify disordered residues (threshold = 0.5)
disordered_residues = disorder_scores >= 0.5
print(f"Disordered positions: {np.where(disordered_residues)[0]}")
```

### Automatic translation from DNA/RNA

```python
# DNA sequence - automatically translated to protein
dna_record = SeqRecord(Seq("ATGAAAGCCCTGGTGTCTTGGGGACGTCCACAAATG"), id="gene1")
dna_record.annotations["molecule_type"] = "DNA"

result = feature.compute_vector(dna_record)
disorder_scores = result["DISORDER_P"]

# Length will be len(DNA)/3 = 12 amino acids (after translation)
print(f"Translated protein length: {len(disorder_scores)}")
```

### Using ORF regions

```python
from biotooler.core.orf_store import attach_orf

# DNA sequence with specific ORF region
dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTCTTGGGGACGTCCACAAATGTAG"), id="gene2")
dna.annotations["molecule_type"] = "DNA"

# Attach ORF (start=0, end=36, length=12 codons)
dna = attach_orf(dna, (0, 36))

# DisorderProfileMetapredict will use attached ORF for translation
feature = DisorderProfileMetapredict(use_orf_if_present=True)
result = feature.compute_vector(dna)

# Disorder scores for the 12 amino acids encoded by the ORF
print(f"ORF disorder scores: {result['DISORDER_P']}")
```

### Computing derived features without re-running predictor

```python
from biotooler.families.disorder.derived import DisorderDerivedScalars

# Assume DISORDER_P has already been computed and stored in record.annotations
# (by DisorderProfileMetapredict or loaded from cache)
record.annotations["DISORDER_P"] = result["DISORDER_P"]

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
    result_t = derived(record)
    print(f"Threshold {t}: {result_t['DISORDER_FRAC']:.2%} disordered, "
          f"longest IDR = {result_t['DISORDER_LONGEST_IDR']}")
```

### Windowed analysis with FeatureSet

```python
from biotooler.features.sets import FeatureSet

# Create feature set with disorder predictor
feature = DisorderProfileMetapredict()
feature_set = FeatureSet({"disorder": feature})

# Compute disorder in sliding windows (50 AA window, 25 AA step)
protein = SeqRecord(Seq("M" * 200), id="long_protein")  # 200 AA protein
protein.annotations["molecule_type"] = "protein"

results = feature_set.compute_windows(
    protein,
    window_size=50,
    step=25
)

# Wide-format output with columns:
# "disorder.DISORDER_P_0" (window at position 0)
# "disorder.DISORDER_P_25" (window at position 25)
# "disorder.DISORDER_P_50" (window at position 50)
# ...
# "disorder.DISORDER_P_150" (last window)

print(results.columns)
print(results[["disorder.DISORDER_P_0", "disorder.DISORDER_P_25"]])
```

### Integration with multiple features

```python
from biotooler.families.disorder import DisorderProfileMetapredict, DisorderDerivedScalars
from biotooler.families.protparam import ProtParamFeature
from biotooler.features.sets import FeatureSet

# Combine disorder with other protein features
feature_set = FeatureSet({
    "disorder": DisorderProfileMetapredict(),
    "protparam": ProtParamFeature(),
})

protein = SeqRecord(Seq("MKALVSWGRPQMTEST"), id="test")
protein.annotations["molecule_type"] = "protein"

# Compute all features together
results = feature_set.compute_windows(protein, window_size=16, step=16)

# Results contain both disorder and protparam features
print(results[["disorder.DISORDER_P_0", 
               "protparam.molecular_weight_0",
               "protparam.gravy_0"]])
```

### Batch processing multiple sequences

```python
# Process multiple sequences efficiently
sequences = [
    SeqRecord(Seq("MKALVSWGRPQM"), id="seq1"),
    SeqRecord(Seq("TESTSEQUENCE"), id="seq2"),
    SeqRecord(Seq("DISORDEREDPROT"), id="seq3"),
]

feature = DisorderProfileMetapredict()
results = []

for seq in sequences:
    seq.annotations["molecule_type"] = "protein"
    result = feature.compute_vector(seq)
    results.append({
        "id": seq.id,
        "length": len(seq),
        "mean_disorder": result["DISORDER_P"].mean(),
        "max_disorder": result["DISORDER_P"].max(),
    })

import pandas as pd
df = pd.DataFrame(results)
print(df)
```

### Identifying highly disordered proteins

```python
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def classify_protein_disorder(sequence: str, threshold: float = 0.5) -> dict:
    """Classify protein disorder characteristics."""
    from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict
    
    record = SeqRecord(Seq(sequence), id="protein")
    record.annotations["molecule_type"] = "protein"
    
    # Predict disorder
    feature = DisorderProfileMetapredict()
    result = feature.compute_vector(record)
    disorder_p = result["DISORDER_P"]
    
    # Compute statistics
    disorder_frac = np.mean(disorder_p >= threshold)
    
    # Classify
    if disorder_frac > 0.4:
        classification = "Highly disordered (IDP)"
    elif disorder_frac > 0.2:
        classification = "Partially disordered"
    else:
        classification = "Mostly structured"
    
    return {
        "sequence": sequence,
        "length": len(sequence),
        "disorder_fraction": disorder_frac,
        "mean_disorder": disorder_p.mean(),
        "classification": classification,
    }

# Example usage
protein = "MKALVSWGRPQMTESTDISORDER"
result = classify_protein_disorder(protein)
print(f"{result['classification']}: {result['disorder_fraction']:.1%} disordered")
```

## Edge cases and validation

### Validated behavior

The disorder family has comprehensive edge case handling:

1. **Empty sequences**: 
   - `DisorderProfileMetapredict`: Returns empty array `np.array([], dtype=np.float64)`
   - `DisorderDerivedScalars`: Returns `DISORDER_FRAC=0.0`, `DISORDER_LONGEST_IDR=0`, `DISORDER_MEAN=nan`, `DISORDER_P95=nan`

2. **Protein sequences**: Passed directly to metapredict (no translation)

3. **DNA/RNA sequences**: Automatically translated via `ensure_protein_record()`
   - Uses attached ORF if present (`use_orf_if_present=True`)
   - Falls back to frame 0 full sequence if no ORF attached
   - Respects `table`, `strip_terminal_stop`, `on_internal_stop` parameters

4. **Stop codons**:
   - Terminal stops: Stripped by default (`strip_terminal_stop=True`)
   - Internal stops: Raises error by default (`on_internal_stop="error"`)

5. **Invalid amino acids**: metapredict may raise `ValueError` or `KeyError` for non-standard amino acids
   - Wrapped and re-raised as `ValueError` with helpful message

6. **Array validation**: 
   - Output length must match protein sequence length
   - Values must be float64 type
   - Returns immediately raise error if metapredict returns wrong length

7. **Derived features**:
   - Requires `DISORDER_P` in `record.annotations`
   - Raises `ValueError` with helpful message if not found
   - Validates that `DISORDER_P` is 1D array

### Known edge cases and limitations

1. **Very short sequences (< 10 residues)**: 
   - Limited context for disorder prediction
   - metapredict may produce less reliable scores
   - Still returns valid output

2. **Transmembrane regions**: 
   - May be incorrectly predicted as disordered due to hydrophobic composition
   - metapredict is trained on soluble proteins
   - Consider using specialized transmembrane predictors for membrane proteins

3. **Low-complexity regions**: 
   - Often predicted as disordered (poly-Q, poly-A stretches)
   - This is biologically accurate in many cases
   - May need domain-specific thresholds

4. **Modified residues**: 
   - Standard metapredict does not account for post-translational modifications
   - Predictions based only on primary sequence
   - Modifications like phosphorylation may affect disorder in vivo

5. **Non-standard amino acids**:
   - Selenocysteine (U) and pyrrolysine (O) may not be supported
   - metapredict trained on 20 standard amino acids
   - May raise `ValueError` if encountered

6. **Genetic code variations**:
   - Translation uses specified `table` parameter (default: 1 = standard code)
   - Mitochondrial, plastid codes supported via `table` parameter
   - Ensure correct table for organism being analyzed

### Performance characteristics

**metapredict (bidirectional LSTM)**:
- **Time complexity**: O(n) for sequence of length n
- **Memory**: Scales with sequence length and model size (~100 MB for model)
- **Typical timing**: Fast inference (~1-10ms per sequence of length 100-1000)
- **GPU acceleration**: metapredict uses PyTorch; GPU can speed up batch predictions
- **Caching**: Results can be cached in `record.annotations["DISORDER_P"]` for reuse by derived features

**DisorderDerivedScalars**:
- **Time complexity**: O(n) for sequence of length n
- **Memory**: Minimal (no model loading, operates on cached DISORDER_P)
- **Typical timing**: Microseconds per sequence
- **Caching benefit**: Does not re-run predictor, reuses DISORDER_P

## Windowing correctness

The disorder family implements the PositionalFeature protocol for correct windowing semantics.

### Position space

Features operate at the **RESIDUE (amino acid) position space**:
- Position indices: 0-based amino acid positions
- Example: Sequence "MKTAY" has positions 0=M, 1=K, 2=T, 3=A, 4=Y
- `DisorderProfileMetapredict.position_space` returns `PositionSpace.RESIDUE`

### Vector computation

Disorder prediction computes a **per-residue vector across the entire sequence**:

1. **Full-sequence input**: Pass complete protein sequence to metapredict
2. **Per-residue scores**: Model returns disorder probability for each amino acid position
3. **Vector output**: numpy array with length equal to sequence length
4. **Window aggregation**: Windowing engine slices the full-context vector and applies aggregation

**Why full-sequence?** Context matters for disorder prediction. The neural network (bidirectional LSTM) considers surrounding amino acids when predicting disorder at each position. Computing on substrings would lose important sequence context and produce incorrect results.

**Upstream API**:
```python
import metapredict

# metapredict.predict_disorder() returns per-residue probabilities
protein_seq = "MKALVSWGRPQMTEST"
disorder_probs = metapredict.predict_disorder(protein_seq)
# Returns: numpy array of length 16 (one value per residue)
# Values in range [0.0, 1.0]
```

**Implementation in `compute_vector`**:
```python
def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
    # Translate DNA/RNA to protein if needed
    protein_record = ensure_protein_record(record, ...)
    protein_seq = str(protein_record.seq)
    
    # Call metapredict on FULL sequence (not sliced)
    disorder_probs = metapredict.predict_disorder(protein_seq)
    
    # Return full-length vector
    return {"DISORDER_P": disorder_probs}
```

### Aggregation strategy

When windowing is applied, per-residue disorder values are aggregated using **mean**:

**Mean disorder score** (implemented):
- Aggregation: `np.mean(disorder_scores[window_start:window_end])`
- Interpretation: Average disorder propensity across a window region
- Use case: Identify regions with overall high or low disorder
- Specified in `vector_keys`: `{"DISORDER_P": AggregationSpec(aggregation_fn=np.mean)}`

**Alternative aggregations** (can be implemented by wrapping the feature):
- **Max disorder**: `np.max(disorder_scores[window])` - Peak disorder in window
- **Fraction disordered**: `np.mean(disorder_scores[window] >= 0.5)` - Binary classification

**Wide-format output**: When using windowing with `FeatureSet.compute_windows()`, output columns follow the pattern `FEATURE_<window_start>`:
```python
# Example: 3 windows at positions 0, 50, 100
# Output columns: "DISORDER_P_0", "DISORDER_P_50", "DISORDER_P_100"
# Each contains the mean disorder score for that window
```

### Testing approach

Windowing correctness validated through:

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
   - Validate that `compute_vector` is called once per sequence

4. **Edge case handling**:
   - Test empty sequences, very short sequences
   - Validate handling of non-standard amino acids
   - Check behavior with unusual sequence compositions

5. **Windowing engine integration**:
   - Verify `compute_vector()` is called once per sequence
   - Test aggregation functions work correctly
   - Validate window boundaries and partial windows
   - Check wide-format output column naming

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
- Support for ANCHOR2 backend for binding region prediction
- Consensus predictions combining multiple methods
- Integration with structure prediction features
- Batch prediction optimization for multiple sequences
- Additional derived features (disorder clusters, region boundaries)

## Contributing a New Backend

This section provides a step-by-step guide for contributors who want to add a new disorder prediction backend (e.g., IUPred3, ANCHOR2, or other tools).

### Overview

The disorder family uses a **lazy import pattern** to keep biotooler lightweight. Each backend is:
1. An optional dependency installed via extras (e.g., `pip install "biotooler[disorder,iupred]"`)
2. Imported only when used (via helper functions in `integration.py`)
3. Registered in the family's `get_features()` function

### Files to Add

To add a new backend (e.g., IUPred3), create the following files:

```
src/biotooler/families/disorder/
├── __init__.py                    # UPDATE: Add new feature to get_features()
├── integration.py                 # UPDATE: Add require_<backend>() helper
├── iupred3_backend.py             # NEW: Backend-specific feature class
└── README.md                      # UPDATE: Document the new backend

tests/families/disorder/
├── test_iupred3_backend.py        # NEW: Backend-specific tests
└── test_integration.py            # UPDATE: Add lazy import tests
```

### Step 1: Create Backend Feature Class

**File**: `src/biotooler/families/disorder/iupred3_backend.py`

**Template** (adapt from `metapredict_backend.py`):

```python
"""IUPred3 backend for disorder prediction."""

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.record import get_molecule_type
from biotooler.core.translation import ensure_protein_record
from biotooler.families.disorder.integration import require_iupred3
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class DisorderProfileIUPred3:
    """Compute per-residue intrinsic disorder using IUPred3.
    
    This feature implements the PositionalFeature protocol to compute protein
    intrinsic disorder predictions using the IUPred3 external tool.
    
    Args:
        mode: IUPred3 mode - "long" (long disorder), "short" (short disorder),
            or "glob" (globular domains) (default: "long")
        table: NCBI genetic code table for translation (default: 1)
        use_orf_if_present: Use attached ORF if present (default: True)
        strip_terminal_stop: Strip terminal stop codon (default: True)
        on_internal_stop: Action for internal stops - "error" or "ignore" (default: "error")
    """
    
    def __init__(
        self,
        *,
        mode: str = "long",
        table: int = 1,
        use_orf_if_present: bool = True,
        strip_terminal_stop: bool = True,
        on_internal_stop: str = "error",
    ):
        if mode not in ("long", "short", "glob"):
            raise ValueError(f"mode must be 'long', 'short', or 'glob', got {mode!r}")
        self.mode = mode
        self.table = table
        self.use_orf_if_present = use_orf_if_present
        self.strip_terminal_stop = strip_terminal_stop
        self.on_internal_stop = on_internal_stop
    
    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space."""
        return PositionSpace.RESIDUE
    
    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specifications for each feature key."""
        return {
            f"DISORDER_P_IUPRED3_{self.mode.upper()}": AggregationSpec(aggregation_fn=np.mean),
        }
    
    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-residue disorder using IUPred3.
        
        Args:
            record: SeqRecord containing DNA, RNA, or protein sequence
            **kwargs: Additional parameters (positions ignored, full computation always done)
        
        Returns:
            Dictionary mapping feature key to numpy array of disorder probabilities
        """
        # Lazy import IUPred3 wrapper
        iupred3 = require_iupred3()
        
        # Translate DNA/RNA to protein if needed
        mol_type = get_molecule_type(record)
        if mol_type.upper() in ("DNA", "RNA"):
            protein_record = ensure_protein_record(
                record,
                table=self.table,
                use_orf_if_present=self.use_orf_if_present,
                strip_terminal_stop=self.strip_terminal_stop,
                on_internal_stop=self.on_internal_stop,
            )
        else:
            protein_record = record
        
        protein_seq = str(protein_record.seq)
        
        # Handle empty sequence
        if len(protein_seq) == 0:
            return {f"DISORDER_P_IUPRED3_{self.mode.upper()}": np.array([], dtype=np.float64)}
        
        # Call IUPred3 (implementation depends on wrapper)
        # This is pseudocode - actual implementation depends on IUPred3 interface
        disorder_probs = iupred3.predict(protein_seq, mode=self.mode)
        
        # Ensure numpy array with correct length
        disorder_probs = np.array(disorder_probs, dtype=np.float64)
        if len(disorder_probs) != len(protein_seq):
            raise ValueError(
                f"IUPred3 returned array of length {len(disorder_probs)}, "
                f"but expected {len(protein_seq)} for record {record.id!r}"
            )
        
        return {f"DISORDER_P_IUPRED3_{self.mode.upper()}": disorder_probs}
```

### Step 2: Add Lazy Import Helper

**File**: `src/biotooler/families/disorder/integration.py`

**Add function**:

```python
def require_iupred3():
    """Require IUPred3 dependency with lazy loading.
    
    Returns:
        The iupred3 wrapper module
    
    Raises:
        ImportError: If IUPred3 is not installed, with installation instructions
    """
    return lazy_import(
        "iupred3",  # or "biotooler_iupred3" if we create a wrapper package
        extra="iupred",
        purpose="computing protein disorder features using IUPred3",
    )
```

### Step 3: Register in Family

**File**: `src/biotooler/families/disorder/__init__.py`

**Update `get_features()`**:

```python
def get_features() -> list[type]:
    """Get list of feature classes provided by this family."""
    from biotooler.families.disorder.derived import DisorderDerivedScalars
    from biotooler.families.disorder.metapredict_backend import DisorderProfileMetapredict
    from biotooler.families.disorder.iupred3_backend import DisorderProfileIUPred3  # NEW
    
    return [
        DisorderProfileMetapredict,
        DisorderProfileIUPred3,  # NEW
        DisorderDerivedScalars,
    ]
```

### Step 4: Update pyproject.toml

**File**: `pyproject.toml`

**Add optional dependency group**:

```toml
[project.optional-dependencies]
disorder = [
    "metapredict>=3.0",
]
iupred = [
    # If IUPred3 is packaged on PyPI:
    "iupred3>=3.0",
    # Or if we create a wrapper:
    # "biotooler-iupred3>=0.1.0",
]
# For installing both:
all-disorder = [
    "biotooler[disorder,iupred]",
]
```

**Note**: If IUPred3 is not available on PyPI, document manual installation steps in README.

### Step 5: Add Tests

**File**: `tests/families/disorder/test_iupred3_backend.py`

**Copy and adapt from** `tests/families/disorder/test_metapredict_backend.py` (if exists, or create):

```python
"""Tests for IUPred3 backend."""

import pytest
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Try importing IUPred3, skip tests if not available
pytest.importorskip("iupred3", reason="IUPred3 not installed")

from biotooler.families.disorder.iupred3_backend import DisorderProfileIUPred3


def test_iupred3_basic():
    """Test basic IUPred3 disorder prediction."""
    feature = DisorderProfileIUPred3(mode="long")
    
    # Simple protein sequence
    record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
    record.annotations["molecule_type"] = "protein"
    
    result = feature.compute_vector(record)
    
    # Check output structure
    assert "DISORDER_P_IUPRED3_LONG" in result
    disorder_p = result["DISORDER_P_IUPRED3_LONG"]
    
    # Check array properties
    assert isinstance(disorder_p, np.ndarray)
    assert len(disorder_p) == 12  # Length of sequence
    assert disorder_p.dtype == np.float64
    
    # Check value range [0, 1]
    assert np.all(disorder_p >= 0.0)
    assert np.all(disorder_p <= 1.0)


def test_iupred3_modes():
    """Test different IUPred3 modes."""
    record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
    record.annotations["molecule_type"] = "protein"
    
    for mode in ("long", "short", "glob"):
        feature = DisorderProfileIUPred3(mode=mode)
        result = feature.compute_vector(record)
        
        expected_key = f"DISORDER_P_IUPRED3_{mode.upper()}"
        assert expected_key in result
        assert len(result[expected_key]) == 12


def test_iupred3_dna_translation():
    """Test automatic DNA to protein translation."""
    feature = DisorderProfileIUPred3()
    
    # DNA sequence encoding "MKAL"
    dna = SeqRecord(Seq("ATGAAAGCCCTG"), id="dna1")
    dna.annotations["molecule_type"] = "DNA"
    
    result = feature.compute_vector(dna)
    disorder_p = result["DISORDER_P_IUPRED3_LONG"]
    
    # Should have 4 values (4 amino acids after translation)
    assert len(disorder_p) == 4


def test_iupred3_empty_sequence():
    """Test handling of empty sequences."""
    feature = DisorderProfileIUPred3()
    
    record = SeqRecord(Seq(""), id="empty")
    record.annotations["molecule_type"] = "protein"
    
    result = feature.compute_vector(record)
    disorder_p = result["DISORDER_P_IUPRED3_LONG"]
    
    assert len(disorder_p) == 0
    assert disorder_p.dtype == np.float64


def test_iupred3_position_space():
    """Test position space property."""
    from biotooler.features.aggregation import PositionSpace
    
    feature = DisorderProfileIUPred3()
    assert feature.position_space == PositionSpace.RESIDUE


def test_iupred3_vector_keys():
    """Test vector_keys property."""
    feature = DisorderProfileIUPred3(mode="long")
    vector_keys = feature.vector_keys
    
    assert "DISORDER_P_IUPRED3_LONG" in vector_keys
    assert vector_keys["DISORDER_P_IUPRED3_LONG"].aggregation_fn == np.mean


def test_iupred3_deterministic():
    """Test that IUPred3 produces deterministic results."""
    feature = DisorderProfileIUPred3()
    
    record = SeqRecord(Seq("MKALVSWGRPQM"), id="test")
    record.annotations["molecule_type"] = "protein"
    
    result1 = feature.compute_vector(record)
    result2 = feature.compute_vector(record)
    
    np.testing.assert_array_equal(
        result1["DISORDER_P_IUPRED3_LONG"],
        result2["DISORDER_P_IUPRED3_LONG"]
    )
```

**File**: `tests/families/disorder/test_integration.py`

**Add lazy import test**:

```python
def test_require_iupred3():
    """Test IUPred3 lazy import helper."""
    from biotooler.families.disorder.integration import require_iupred3
    
    # Should not raise if IUPred3 is installed
    iupred3 = require_iupred3()
    assert iupred3 is not None
```

### Step 6: Update Documentation

**File**: `src/biotooler/families/disorder/README.md`

Add sections:
1. Update "What this family provides" with new backend
2. Add backend-specific examples
3. Document installation: `pip install "biotooler[disorder,iupred]"`
4. Update "Future work" to mark IUPred3 as implemented

### Testing Your Backend

**Run tests**:
```bash
# Install with backend
pip install -e ".[disorder,iupred]"

# Run backend-specific tests
pytest tests/families/disorder/test_iupred3_backend.py -v

# Run all disorder tests
pytest tests/families/disorder/ -v
```

**Test lazy loading**:
```python
# Should not import backend on family import
from biotooler.families import disorder
# iupred3 not imported yet

# Should import backend only when feature is used
feature = disorder.DisorderProfileIUPred3()
result = feature.compute_vector(record)
# Now iupred3 is imported
```

**Test windowing**:
```python
from biotooler.features import FeatureSet

feature_set = FeatureSet({"iupred3": DisorderProfileIUPred3()})
results = feature_set.compute_windows(
    record,
    window_size=50,
    step=25
)
# Verify wide-format output with DISORDER_P_IUPRED3_LONG_0, DISORDER_P_IUPRED3_LONG_25, etc.
```

### Key Implementation Patterns

**1. Lazy Import Pattern**:
- ✅ DO: Import backend inside `compute_vector()` using `require_<backend>()`
- ❌ DON'T: Import backend at module top-level

**2. Translation Pattern**:
- ✅ DO: Use `ensure_protein_record()` with standard translation parameters
- ✅ DO: Support DNA/RNA/protein input consistently across backends
- ❌ DON'T: Implement custom translation logic

**3. PositionalFeature Protocol**:
- ✅ DO: Implement `position_space`, `vector_keys`, `compute_vector`
- ✅ DO: Return full-length vectors (length = sequence length)
- ✅ DO: Use appropriate aggregation (usually `np.mean` for disorder)
- ❌ DON'T: Slice sequences before computing (loses context)

**4. Error Handling**:
- ✅ DO: Catch backend-specific errors and raise `ValueError` with helpful messages
- ✅ DO: Validate output length matches input sequence length
- ✅ DO: Handle empty sequences gracefully

**5. Naming Conventions**:
- Feature class: `DisorderProfile<Backend>`
- Feature key: `DISORDER_P_<BACKEND>_<MODE>` (uppercase)
- Module file: `<backend>_backend.py` (lowercase)
- Test file: `test_<backend>_backend.py`

### Common Pitfalls

**❌ Importing backend at module level**:
```python
# BAD: Loads backend even if never used
import iupred3

class DisorderProfileIUPred3:
    def compute_vector(self, record):
        return iupred3.predict(...)
```

**✅ Lazy import pattern**:
```python
# GOOD: Backend loaded only when feature is used
class DisorderProfileIUPred3:
    def compute_vector(self, record):
        iupred3 = require_iupred3()  # Lazy import
        return iupred3.predict(...)
```

**❌ Slicing sequence before prediction**:
```python
# BAD: Loses context, incorrect results
def compute_vector(self, record, positions=None):
    if positions:
        seq = str(record.seq)[positions[0]:positions[-1]]  # Wrong!
    return {"DISORDER_P": predict(seq)}
```

**✅ Always use full sequence**:
```python
# GOOD: Predict on full sequence, let windowing engine slice vector
def compute_vector(self, record, **kwargs):
    seq = str(record.seq)  # Full sequence
    full_vector = predict(seq)
    return {"DISORDER_P": full_vector}  # Return full vector
```

### Getting Help

- Review existing backends: `metapredict_backend.py`, `derived.py`
- Review other families: `chimera/`, `protparam/`, `basic_stats/`
- Check tests for implementation patterns
- Ask maintainers: @itamar-menuhin (disorder family owner)

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
