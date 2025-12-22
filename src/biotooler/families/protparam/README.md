# ProtParam Family

Protein physicochemical parameter calculations using Bio.SeqUtils.ProtParam.

## Overview

This family provides the `ProtParamFeature` class for computing protein properties including:
- Molecular weight
- Isoelectric point (pI)
- GRAVY (Grand Average of Hydropathy)
- Instability index
- Aromaticity

The feature automatically translates DNA/RNA sequences to protein before computing parameters.

## Features

### ProtParamFeature

Computes five protein physicochemical parameters:

#### 1. Molecular Weight

**Description**: The sum of the atomic masses of all atoms in the protein, expressed in Daltons (Da) or grams per mole (g/mol).

**Calculation**: 
```
MW = Σ (count_i × MW_i) - (n-1) × MW_H2O
```
Where:
- `count_i` is the count of amino acid type i
- `MW_i` is the molecular weight of amino acid i
- `n` is the number of amino acids
- `MW_H2O = 18.015` Da (water lost in peptide bond formation)

**Interpretation**: 
- Useful for protein identification and purification
- SDS-PAGE migration distance is inversely proportional to log(MW)
- Typical range: 5-500 kDa for most proteins

**Key**: `molecular_weight` (float, in Daltons)

#### 2. Isoelectric Point (pI)

**Description**: The pH at which a protein carries no net electrical charge. At this pH, the protein has equal numbers of positive and negative charges.

**Calculation**: 
The pI is calculated iteratively by finding the pH where the net charge equals zero:
```
Net charge = Σ (positive charges) - Σ (negative charges)
```
Charges depend on pKa values of ionizable groups:
- N-terminus: pKa ≈ 9.6
- C-terminus: pKa ≈ 2.4
- Acidic residues (D, E): pKa ≈ 3.9-4.3
- Basic residues (K, R, H): pKa ≈ 6.0-12.5

**Interpretation**:
- Proteins are positively charged at pH < pI
- Proteins are negatively charged at pH > pI
- Used for isoelectric focusing (IEF) separation
- Affects protein solubility and stability
- Typical range: 4-12 for most proteins

**Key**: `isoelectric_point` (float, pH units)

#### 3. GRAVY (Grand Average of Hydropathy)

**Description**: A measure of the overall hydrophobicity or hydrophilicity of a protein, calculated as the average hydropathy value of all amino acids.

**Calculation**:
```
GRAVY = Σ (H_i) / n
```
Where:
- `H_i` is the hydropathy value of amino acid i (Kyte-Doolittle scale)
- `n` is the total number of amino acids

Kyte-Doolittle hydropathy scale (selected values):
- Most hydrophobic: Ile (+4.5), Val (+4.2), Leu (+3.8)
- Most hydrophilic: Arg (-4.5), Lys (-3.9), Asp (-3.5)

**Interpretation**:
- **Negative values**: Hydrophilic (water-loving), typically soluble proteins
- **Positive values**: Hydrophobic (water-fearing), may be membrane proteins
- Range typically: -2.0 to +2.0
- Membrane proteins often have GRAVY > 0
- Secreted/soluble proteins often have GRAVY < 0

**Key**: `gravy` (float, dimensionless)

**Reference**: Kyte, J. and Doolittle, R.F. (1982) A simple method for displaying the hydropathic character of a protein. J. Mol. Biol. 157:105-132

#### 4. Instability Index

**Description**: An estimate of protein stability in vitro, based on the occurrence of destabilizing dipeptides in the sequence.

**Calculation**:
```
II = (10/L) × Σ (DIWV_ij)
```
Where:
- `L` is the length of the sequence
- `DIWV_ij` is the instability weight value for dipeptide composed of amino acids i and j
- Sum is over all dipeptides in the sequence

The DIWV values are derived from the stability analysis of 12 unstable and 32 stable proteins.

**Interpretation**:
- **II < 40**: Protein is predicted to be **stable** in vitro
- **II > 40**: Protein is predicted to be **unstable** in vitro
- Based on statistical analysis of protein half-life in test tubes
- Not necessarily predictive of in vivo stability
- Useful for expression system selection and purification planning

**Key**: `instability_index` (float, dimensionless)

**Reference**: Guruprasad, K., Reddy, B.V.B., and Pandit, M.W. (1990) Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. Protein Engineering 4:155-161

#### 5. Aromaticity

**Description**: The relative frequency of aromatic amino acids (Phe, Trp, Tyr) in the protein sequence.

**Calculation**:
```
Aromaticity = (N_Phe + N_Trp + N_Tyr) / N_total
```
Where:
- `N_Phe`, `N_Trp`, `N_Tyr` are counts of phenylalanine, tryptophan, and tyrosine
- `N_total` is the total number of amino acids

**Interpretation**:
- Range: 0.0 to 1.0 (0% to 100%)
- Typical proteins: 0.05-0.15 (5-15%)
- High aromaticity may indicate:
  - DNA/RNA binding proteins
  - Protein-protein interaction interfaces
  - Structural proteins with π-π stacking
- Aromatic residues contribute to:
  - UV absorption at 280 nm (used for protein quantification)
  - Protein folding through π-π interactions
  - Ligand binding sites

**Key**: `aromaticity` (float, fraction 0.0-1.0)

## Usage

### Basic Usage

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

print(result["molecular_weight"])  # Molecular weight in Da
print(result["isoelectric_point"]) # pI
print(result["gravy"])             # GRAVY
print(result["instability_index"]) # Instability index
print(result["aromaticity"])       # Aromaticity
```

### DNA/RNA Sequences

The feature automatically translates DNA/RNA sequences:

```python
# DNA sequence
dna = SeqRecord(Seq("ATGAAAGCCCTGGTGTCTTGGGGACGT"), id="dna1")
dna.annotations["molecule_type"] = "DNA"
result = feature(dna)  # Automatically translated
```

### Translation Options

```python
# Custom translation table
feature = ProtParamFeature(table=11)  # Use bacterial genetic code

# Explicit ORF region
from biotooler.core.orf_store import attach_orf
record = attach_orf(dna, (0, 27))
result = feature(record)  # Uses attached ORF

# Allow internal stop codons
feature = ProtParamFeature(on_internal_stop="ignore")
```

## Example Output

For the protein sequence "MKALVSWGR":
```python
{
    'molecular_weight': 1047.27,    # Da
    'isoelectric_point': 11.00,      # pH
    'gravy': 0.1333,                 # Slightly hydrophobic
    'instability_index': -2.63,      # Stable (< 40)
    'aromaticity': 0.1111            # 11.1% aromatic residues (1/9 = W)
}
```

## Dependencies

This family has no additional dependencies beyond the core biotooler requirements:
- `biopython` (includes Bio.SeqUtils.ProtParam)

## Implementation Notes

- Uses lazy imports to avoid loading ProtParam at module import time
- Automatically handles DNA/RNA to protein translation via `ensure_protein_record()`
- By default, strips terminal stop codons and errors on internal stops
- All ProtParam calculations are performed on the translated protein sequence

## References

### Primary References
- **Kyte-Doolittle Hydropathy**: Kyte, J. and Doolittle, R.F. (1982) "A simple method for displaying the hydropathic character of a protein." *Journal of Molecular Biology* 157:105-132. [DOI: 10.1016/0022-2836(82)90515-0](https://doi.org/10.1016/0022-2836(82)90515-0)

- **Instability Index**: Guruprasad, K., Reddy, B.V.B., and Pandit, M.W. (1990) "Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence." *Protein Engineering* 4:155-161. [DOI: 10.1093/protein/4.2.155](https://doi.org/10.1093/protein/4.2.155)

### Additional Resources
- [Biopython ProtParam documentation](https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html)
- [ExPASy ProtParam tool](https://web.expasy.org/protparam/) - Web interface with similar calculations
- [Amino acid scales database (AAindex)](https://www.genome.jp/aaindex/) - Comprehensive collection of amino acid indices and scales
