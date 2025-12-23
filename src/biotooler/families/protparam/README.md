# ProtParam Family

## What this family provides

The protparam family provides feature computation for analyzing protein physicochemical properties using Bio.SeqUtils.ProtParam. It computes five fundamental parameters: molecular weight, isoelectric point (pI), GRAVY (Grand Average of Hydropathy), instability index, and aromaticity. The feature automatically translates DNA/RNA sequences to protein before computing parameters.

## Intuition

Protein physicochemical properties provide fundamental insights into protein behavior, stability, and function. These parameters can be computed directly from amino acid sequence and serve multiple purposes:

1. **Molecular weight**: Essential for protein identification, purification planning, and SDS-PAGE interpretation
2. **Isoelectric point (pI)**: Critical for protein separation techniques (IEF, ion exchange chromatography) and predicting protein behavior at different pH levels
3. **GRAVY (hydropathy)**: Indicates whether a protein is likely membrane-associated (hydrophobic) or soluble (hydrophilic)
4. **Instability index**: Predicts protein stability in vitro, useful for expression system planning
5. **Aromaticity**: Related to UV absorption properties and structural features like π-π stacking

These parameters are widely used in:
- **Protein identification and characterization**
- **Expression and purification strategy design**
- **Comparative proteomics and evolutionary studies**
- **Structure-function relationship analysis**

## Mathematical formulation

This family computes five physicochemical parameters using established algorithms from Bio.SeqUtils.ProtParam:

### 1. Molecular Weight

The sum of the atomic masses of all atoms in the protein, expressed in Daltons (Da).

**Formula**: 
```
MW = Σ (count_i × MW_i) - (n-1) × MW_H2O
```
Where:
- `count_i` = count of amino acid type i
- `MW_i` = molecular weight of amino acid i  
- `n` = number of amino acids
- `MW_H2O = 18.015` Da (water lost in peptide bond formation)

**Range**: Typically 5-500 kDa for most proteins

**Interpretation**: 
- Used for protein identification and purification planning
- SDS-PAGE migration distance is inversely proportional to log(MW)

### 2. Isoelectric Point (pI)

The pH at which a protein carries no net electrical charge.

**Formula**: 
Calculated iteratively by finding the pH where net charge equals zero:
```
Net charge = Σ (positive charges) - Σ (negative charges)
```

Ionizable groups with pKa values:
- N-terminus: pKa ≈ 9.6
- C-terminus: pKa ≈ 2.4
- Acidic residues (D, E): pKa ≈ 3.9-4.3
- Basic residues (K, R, H): pKa ≈ 6.0-12.5

**Range**: Typically 4-12 for most proteins

**Interpretation**:
- Proteins are positively charged at pH < pI
- Proteins are negatively charged at pH > pI
- Critical for isoelectric focusing (IEF) and ion exchange chromatography
- Affects protein solubility and stability at different pH levels

### 3. GRAVY (Grand Average of Hydropathy)

Average hydropathy value across all amino acids in the protein.

**Formula**:
```
GRAVY = Σ (H_i) / n
```
Where:
- `H_i` = hydropathy value of amino acid i (Kyte-Doolittle scale)
- `n` = total number of amino acids

Kyte-Doolittle scale (selected values):
- Most hydrophobic: Ile (+4.5), Val (+4.2), Leu (+3.8)
- Most hydrophilic: Arg (-4.5), Lys (-3.9), Asp (-3.5)

**Range**: Typically -2.0 to +2.0

**Interpretation**:
- **Negative values**: Hydrophilic (water-loving), typically soluble proteins
- **Positive values**: Hydrophobic (water-fearing), may be membrane proteins
- Membrane proteins often have GRAVY > 0
- Secreted/soluble proteins often have GRAVY < 0

**Reference**: Kyte & Doolittle (1982) J. Mol. Biol. 157:105-132

### 4. Instability Index

Estimate of protein stability in vitro based on destabilizing dipeptides.

**Formula**:
```
II = (10/L) × Σ (DIWV_ij)
```
Where:
- `L` = sequence length
- `DIWV_ij` = instability weight for dipeptide (amino acids i, j)

DIWV values derived from analysis of 12 unstable and 32 stable proteins.

**Range**: Typically 0-100+

**Interpretation**:
- **II < 40**: Protein predicted to be **stable** in vitro
- **II > 40**: Protein predicted to be **unstable** in vitro
- Based on test tube half-life analysis
- Not necessarily predictive of in vivo stability
- Useful for expression planning and purification strategy

**Reference**: Guruprasad et al. (1990) Protein Engineering 4:155-161

### 5. Aromaticity

Relative frequency of aromatic amino acids (Phe, Trp, Tyr).

**Formula**:
```
Aromaticity = (N_Phe + N_Trp + N_Tyr) / N_total
```
Where:
- `N_Phe`, `N_Trp`, `N_Tyr` = counts of aromatic amino acids
- `N_total` = total amino acids

**Range**: 0.0 to 1.0 (0% to 100%), typically 0.05-0.15 for most proteins

**Interpretation**:
- High aromaticity may indicate:
  - DNA/RNA binding proteins
  - Protein-protein interaction interfaces
  - Structural proteins with π-π stacking
- Aromatic residues contribute to:
  - UV absorption at 280 nm (protein quantification)
  - Protein folding through π-π interactions
  - Ligand binding sites

## Features and output schema

### Input
- No additional configuration required beyond translation options
- Optional translation parameters: `table`, `on_internal_stop`

### Output
Dictionary mapping feature names to scalar float values:
```python
{
    "molecular_weight": float,    # Molecular weight in Daltons
    "isoelectric_point": float,   # pI in pH units
    "gravy": float,               # Hydropathy index (dimensionless)
    "instability_index": float,   # Stability index (dimensionless)
    "aromaticity": float          # Fraction of aromatic residues (0.0-1.0)
}
```

### Translation Handling
- DNA/RNA sequences are automatically translated to protein
- Supports custom genetic code tables via `table` parameter
- Handles ORF regions: explicit, attached, or full sequence frame 0
- By default strips terminal stops and errors on internal stops

## References

### Primary References
- **Kyte-Doolittle Hydropathy**: Kyte, J. and Doolittle, R.F. (1982) "A simple method for displaying the hydropathic character of a protein." *Journal of Molecular Biology* 157:105-132. https://doi.org/10.1016/0022-2836(82)90515-0

- **Instability Index**: Guruprasad, K., Reddy, B.V.B., and Pandit, M.W. (1990) "Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence." *Protein Engineering* 4:155-161. https://doi.org/10.1093/protein/4.2.155

### Additional Resources
- Biopython ProtParam: https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html
- ExPASy ProtParam tool: https://web.expasy.org/protparam/
- AAindex database: https://www.genome.jp/aaindex/

## Upstream library links

This family uses Bio.SeqUtils.ProtParam from Biopython (no additional dependencies):
- **Biopython documentation**: https://biopython.org/
- **ProtParam module**: https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html

## Examples

### Basic usage with protein sequence

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.families.protparam import ProtParamFeature

# Create feature
feature = ProtParamFeature()

# Protein sequence
protein = SeqRecord(Seq("MKALVSWGR"), id="protein1")
protein.annotations["molecule_type"] = "protein"
result = feature(protein)

print(result["molecular_weight"])  # 1047.27 Da
print(result["isoelectric_point"]) # 11.00
print(result["gravy"])             # 0.1333 (slightly hydrophobic)
print(result["instability_index"]) # -2.63 (stable)
print(result["aromaticity"])       # 0.1111 (11.1% aromatic)
```

### Automatic translation of DNA sequences

```python
# DNA sequence - automatically translated
dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTCTTGGGGACGT"), id="dna1")
dna.annotations["molecule_type"] = "DNA"
result = feature(dna)

# Translation and parameter computation happen automatically
print(result["molecular_weight"])
```

### Custom translation options

```python
# Use bacterial genetic code
feature = ProtParamFeature(table=11)

# Explicit ORF region
from biotooler.core.orf_store import attach_orf
record = attach_orf(dna, (0, 27))
result = feature(record)

# Allow internal stop codons (if needed)
feature = ProtParamFeature(on_internal_stop="ignore")
```

### Example output interpretation

For the protein "MKALVSWGR":
```python
{
    'molecular_weight': 1047.27,    # Medium-small protein
    'isoelectric_point': 11.00,      # Basic protein (pI > 7)
    'gravy': 0.1333,                 # Slightly hydrophobic
    'instability_index': -2.63,      # Stable (< 40)
    'aromaticity': 0.1111            # 11.1% aromatic (W = 1/9)
}
```

**Interpretation**:
- Small peptide (~1 kDa)
- Highly basic (pI=11), will be positively charged at physiological pH
- Slightly hydrophobic (GRAVY > 0), but still likely soluble
- Very stable in vitro (II < 0)
- Normal aromatic content

## Edge cases and validation

### Validated behavior
1. **Protein sequences**: Passed directly to ProtParam (no translation)
2. **DNA/RNA sequences**: Automatically translated using `ensure_protein_record()`
3. **Empty protein sequences**: May raise errors from ProtParam
4. **Stop codons**: Stripped by default during translation
5. **Non-standard amino acids**: Handled by ProtParam's internal logic

### Known limitations
1. **Modified amino acids**: Does not handle post-translational modifications
2. **Protein complexes**: Calculations assume monomeric proteins
3. **Context-dependent properties**: pI and stability are sequence-based predictions
4. **Very short sequences**: Statistics may be unreliable for peptides < 5 residues

## Windowing correctness

### Position space

This family operates at the **protein sequence level** (amino acid sequences). Features are computed on complete protein sequences, not per-residue values. The underlying ProtParam algorithms analyze global sequence properties.

### Vector computation

Not applicable - this family does not use the PositionalFeature protocol. Instead, it computes protein physicochemical parameters by calling Biopython's `ProtParam.ProteinAnalysis` class with complete protein sequences:

1. **Protein input**: If input is already protein, analyze directly
2. **DNA/RNA input**: Translate to protein using `ensure_protein_record()` before analysis
3. **Upstream library call**: Create `ProteinAnalysis(protein_seq)` with full sequence
4. **Feature computation**: Call methods like `molecular_weight()`, `isoelectric_point()`, `gravy()`, etc.

The upstream ProtParam library requires full protein sequences to compute properties like isoelectric point (depends on charged residue distribution) and instability index (uses empirical dipeptide instability values). These are inherently sequence-level properties that cannot be meaningfully computed per-residue.

### Aggregation strategy

Not applicable - ProtParam features are naturally sequence-level properties:
- **Molecular weight**: Sum of all amino acid masses
- **Isoelectric point**: pH where net charge is zero (requires full sequence charge distribution)
- **GRAVY**: Average hydropathy across all residues
- **Instability index**: Weighted sum of dipeptide instability values
- **Aromaticity**: Fraction of aromatic residues (Phe, Trp, Tyr)

These features are computed by the upstream Biopython library and returned as single scalar values per sequence.

### Testing approach

Windowing correctness for this family focuses on:
1. **Full sequence usage**: Verify `ProteinAnalysis` is instantiated with complete protein sequences
2. **Translation correctness**: For DNA/RNA inputs, ensure proper translation to protein before ProtParam analysis
3. **ORF handling**: When ORFs are present, verify correct extraction and translation
4. **Upstream library correctness**: Results should match direct Biopython ProtParam calculations
5. **Edge cases**: Handle terminal stops, internal stops, and ambiguous amino acids correctly

Tests validate that windowing (if applied to protein features) uses complete protein sequences per window, and that ProtParam is called with proper full-sequence context for each analysis.

## Maintenance notes

### Dependencies
- `biopython` (includes Bio.SeqUtils.ProtParam) - core dependency
- Internal `biotooler.core.translation.ensure_protein_record()` for DNA/RNA handling

### Implementation notes
- Uses lazy imports to avoid loading ProtParam at module import time
- Translation handled by shared `ensure_protein_record()` helper
- All five parameters computed in single `ProtParamFeature` call
- No incremental/rolling computation support (protein-level analysis)

### Design decisions
1. **Single feature class**: All five parameters computed together (efficient)
2. **Automatic translation**: Simplifies API for mixed DNA/RNA/protein workflows
3. **No caching**: ProtParam calculations are fast enough not to require caching
4. **Lazy imports**: Keeps module loading lightweight
