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

- **molecular_weight**: Molecular weight in Daltons
- **isoelectric_point**: Theoretical isoelectric point (pI)
- **gravy**: Grand Average of Hydropathy (negative = hydrophilic, positive = hydrophobic)
- **instability_index**: Instability index (>40 indicates unstable protein in vivo)
- **aromaticity**: Fraction of aromatic amino acids (F, W, Y)

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

## Dependencies

This family has no additional dependencies beyond the core biotooler requirements:
- `biopython` (includes Bio.SeqUtils.ProtParam)

## Implementation Notes

- Uses lazy imports to avoid loading ProtParam at module import time
- Automatically handles DNA/RNA to protein translation via `ensure_protein_record()`
- By default, strips terminal stop codons and errors on internal stops
- All ProtParam calculations are performed on the translated protein sequence

## References

- Bio.SeqUtils.ProtParam documentation: https://biopython.org/docs/1.75/api/Bio.SeqUtils.ProtParam.html
- Protein parameter calculations based on standard bioinformatics methods
