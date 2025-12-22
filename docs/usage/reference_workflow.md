# Reference-Driven Workflow

This guide demonstrates how to use `ReferenceSequenceSet` as a shared reference API across biotooler, with a focus on codon bias analysis. The reference set provides a unified way to manage CDS and protein sequences that can be used by multiple tools including `codon_bias` and `chimera`.

## Overview

`ReferenceSequenceSet` is a container for reference CDS (DNA/RNA) and protein sequences that:
- Loads sequences from FASTA or CSV files
- Normalizes and validates sequences (uppercase, U→T conversion)
- Provides CDS and protein sequence access with consistent defaults
- Caches results for efficiency
- Supports translation with configurable genetic codes

## Loading Reference Sequences

### From FASTA Files

The most common way to create a reference set is from a FASTA file containing CDS sequences:

```python
from biotooler.core.reference_sequences import ReferenceSequenceSet

# Load CDS sequences from a FASTA file
ref_set = ReferenceSequenceSet.from_fasta("highly_expressed_genes.fasta")

# Access the loaded sequences
print(f"Loaded {len(ref_set.cds)} CDS sequences")
print(f"Sequence IDs: {list(ref_set.cds.keys())}")
```

You can also provide separate CDS and protein FASTA files:

```python
# Load both CDS and protein sequences
ref_set = ReferenceSequenceSet.from_fasta(
    cds_fasta_path="genes_cds.fasta",
    protein_fasta_path="genes_protein.fasta",
    genetic_code_table=1  # Standard genetic code (default)
)
```

### From CSV Files

For tabular data, use the CSV loader:

```python
# Load from CSV with column mapping
ref_set = ReferenceSequenceSet.from_csv(
    csv_path="reference_genes.csv",
    id_column="gene_id",
    cds_column="nucleotide_sequence",
    protein_column="protein_sequence",  # Optional
    genetic_code_table=1
)
```

### Direct Construction

You can also create a reference set programmatically:

```python
# Create from dictionaries
ref_set = ReferenceSequenceSet(
    cds={
        "gene1": "ATGATGATGATGATGATGATG",
        "gene2": "ATGGGCGGCTAA",
        "gene3": "ATGAAACTGTAG",
    },
    proteins=None,  # Will be derived from CDS if needed
    genetic_code_table=1
)
```

## Using Reference Sets with Codon Bias

The `CodonBiasFeature.from_reference()` factory method makes it easy to build codon bias models from a reference set.

### Basic Codon Bias Analysis

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.reference_sequences import ReferenceSequenceSet
from biotooler.families.codon_bias import CodonBiasFeature

# Load reference set from highly expressed genes
ref_set = ReferenceSequenceSet.from_fasta("highly_expressed_genes.fasta")

# Build codon bias feature with multiple metrics
feature = CodonBiasFeature.from_reference(
    ref_set,
    ["CAI", "ENC", "FOP"]  # Codon Adaptation Index, Effective Number of Codons, Frequency of Optimal Codons
)

# Analyze a test sequence
test_seq = "ATGATGATGATGATGATG"
record = SeqRecord(Seq(test_seq), id="test_gene")
result = feature(record)

print(result)
# {'CAI': 1.0, 'ENC': 20.0, 'FOP': 1.0}
```

### Automatic Reference Handling

The `from_reference()` method automatically:
1. Concatenates CDS sequences from the reference set
2. Detects which scores require `ref_seq` parameter (e.g., CAI, FOP)
3. Passes the concatenated reference to scores that need it
4. Instantiates scores without reference when not needed (e.g., ENC)

```python
# Mix scores that do and don't require reference sequences
feature = CodonBiasFeature.from_reference(
    ref_set,
    ["CAI", "ENC"]  # CAI needs ref_seq, ENC doesn't
)
# Both work correctly
```

### Score-Specific Parameters

Pass custom parameters to individual scores using `score_kwargs`:

```python
feature = CodonBiasFeature.from_reference(
    ref_set,
    ["CAI", "ENC"],
    score_kwargs={
        "CAI": {"genetic_code": 11},  # Bacterial genetic code
        "ENC": {"bg_correction": True, "robust": False}
    }
)
```

## CDS-Only Workflow

For workflows that work directly with CDS sequences, use the normalized CDS strings:

```python
from biotooler.core.reference_sequences import ReferenceSequenceSet

# Load reference sequences
ref_set = ReferenceSequenceSet.from_fasta("reference_genes.fasta")

# Get normalized CDS strings
# - Converted to uppercase
# - U replaced with T (RNA to DNA)
# - Validated to be multiple of 3
cds_sequences = ref_set.cds_strings(require_multiple_of_three=True)

# Use for codon bias models
concatenated_cds = "".join(cds_sequences)
print(f"Total CDS length: {len(concatenated_cds)} nt")

# Build codon bias feature directly from reference
feature = CodonBiasFeature.from_reference(ref_set, ["CAI", "FOP"])
```

**Default CDS normalization:**
- Sequences are converted to uppercase
- U is replaced with T (RNA sequences converted to DNA)
- `require_multiple_of_three=True` (default): validates length divisible by 3

## CDS→Protein Derived Workflow

When protein sequences are needed but not provided, the reference set can derive them from CDS:

```python
from biotooler.core.reference_sequences import ReferenceSequenceSet

# Create reference set with CDS only
ref_set = ReferenceSequenceSet.from_fasta(
    "reference_genes.fasta",
    genetic_code_table=1  # Standard genetic code
)

# Derive protein sequences from CDS
protein_sequences = ref_set.protein_strings(
    strip_terminal_stop=False,      # Keep terminal stop codon (*) by default
    error_on_internal_stop=True     # Error if internal stop codons found (default)
)

print(f"Translated {len(protein_sequences)} proteins")
for i, protein in enumerate(protein_sequences[:3]):
    print(f"Protein {i}: {protein[:50]}...")  # Show first 50 amino acids
```

### Translation Parameters

The `protein_strings()` method provides control over translation:

**Default behavior:**
```python
# Default: keep terminal stops, error on internal stops
proteins = ref_set.protein_strings()
# Returns: ["MKG*", "MVLA*", "MFRT*"]
```

**Strip terminal stop codons:**
```python
# Remove terminal stop codons for downstream analysis
proteins = ref_set.protein_strings(strip_terminal_stop=True)
# Returns: ["MKG", "MVLA", "MFRT"]
```

**Allow internal stop codons:**
```python
# Useful for incomplete ORFs or pseudogenes
proteins = ref_set.protein_strings(error_on_internal_stop=False)
# Returns proteins even if they contain internal stops
```

### Custom Genetic Code Tables

Use alternative genetic codes for different organisms:

```python
# Bacterial/Archaeal/Plant Plastid code (table 11)
ref_set = ReferenceSequenceSet.from_fasta(
    "bacterial_genes.fasta",
    genetic_code_table=11
)

# Yeast mitochondrial code (table 3)
ref_set_mito = ReferenceSequenceSet.from_fasta(
    "yeast_mito_genes.fasta",
    genetic_code_table=3
)

# Pass custom code to codon bias models
feature = CodonBiasFeature.from_reference(
    ref_set,
    ["CAI"],
    score_kwargs={"CAI": {"genetic_code": 11}}
)
```

Common genetic code tables:
- **1** (default): Standard genetic code
- **11**: Bacterial, Archaeal, and Plant Plastid code
- **2**: Vertebrate mitochondrial code
- **3**: Yeast mitochondrial code
- **4**: Mold, Protozoan, and Coelenterate mitochondrial code

See [NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for full list.

## Default Settings Summary

### Translation Defaults

| Parameter | Default Value | Purpose |
|-----------|--------------|---------|
| `genetic_code_table` | `1` | Standard genetic code for translation |
| `strip_terminal_stop` | `False` | Keep terminal stop codons (*) in proteins |
| `error_on_internal_stop` | `True` | Raise error if internal stop codons found |

### CDS Normalization

| Operation | Behavior |
|-----------|----------|
| Case normalization | Convert to uppercase |
| RNA→DNA conversion | Replace U with T |
| Length validation | Check divisible by 3 when `require_multiple_of_three=True` (default) |

### ORF Policy

When working with ORFs and reference sequences:

- **ORF finding**: Use `find_orf_candidates()` to enumerate all in-frame start/stop codon pairs
- **ORF coordinates**: 0-based and end-exclusive (Python slice convention)
- **Stop codons**: ATG start codon required, stop codons (TAA, TAG, TGA) mark end
- **Frame**: ORFs are in-frame (start and end positions differ by multiple of 3)

```python
from biotooler.core.orf_candidates import find_orf_candidates

# Find all ORF candidates in a sequence
seq = "ATGAAATAGATGCCCTAGTAA"
orfs = find_orf_candidates(seq)
print(orfs)  # [(0, 9), (0, 18), (0, 21), (9, 18), (9, 21)]

# Select specific ORF for analysis
from biotooler.core.orf_store import select_orf_by_index
orf = select_orf_by_index(orfs, 2)  # Select longest ORF at index 2
```

## Complete Reference Workflow Example

Here's a complete workflow combining reference loading, codon bias analysis, and ORF handling:

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from biotooler.core.reference_sequences import ReferenceSequenceSet
from biotooler.core.orf_candidates import find_orf_candidates
from biotooler.families.codon_bias import CodonBiasFeature
from biotooler.features import FeatureSet

# 1. Load reference sequences from FASTA
ref_set = ReferenceSequenceSet.from_fasta(
    "highly_expressed_genes.fasta",
    genetic_code_table=1
)

print(f"Loaded {len(ref_set.cds)} reference sequences")

# 2. Build codon bias feature from reference
feature = CodonBiasFeature.from_reference(
    ref_set,
    ["CAI", "ENC", "FOP"],
    names=["CodonAdaptIndex", "EffNumCodons", "FreqOptimalCodons"]
)

# 3. Analyze a test sequence
test_seq = "ATGATGATGATGATGATGATGATGATGATGTAA"  # 33 nt with stop codon
record = SeqRecord(Seq(test_seq), id="test_gene")

# 3a. Compute global features
global_result = feature(record)
print(f"Global CAI: {global_result['CodonAdaptIndex']:.3f}")

# 3b. Find ORFs and compute windowed features
orfs = find_orf_candidates(record)
print(f"Found {len(orfs)} ORF candidates")

# Use the longest ORF
longest_orf = max(orfs, key=lambda x: x[1] - x[0])

# Create feature set for windowed analysis
fs = FeatureSet(feature, name="codon_bias")

# Compute features in sliding windows across the ORF
windowed_result = fs.compute_orf_windows(
    record,
    orf=longest_orf,
    window_nt=9,   # 3 codons per window
    step_nt=3      # Step by 1 codon (must be multiple of 3)
)

print(windowed_result)
# DataFrame with columns:
# - record_id, orf_start, orf_end
# - codon_bias.CodonAdaptIndex_0, codon_bias.CodonAdaptIndex_3, ...
# - codon_bias.EffNumCodons_0, codon_bias.EffNumCodons_3, ...
# - codon_bias.FreqOptimalCodons_0, codon_bias.FreqOptimalCodons_3, ...

# 4. Optional: Derive protein sequences for other analyses
proteins = ref_set.protein_strings(strip_terminal_stop=True)
print(f"Derived {len(proteins)} protein sequences")
```

## Validation and Error Handling

The reference set provides comprehensive validation:

```python
# Empty reference set
try:
    empty_ref = ReferenceSequenceSet(cds={})
    feature = CodonBiasFeature.from_reference(empty_ref, ["CAI"])
except ValueError as e:
    print(e)
    # "ReferenceSequenceSet must contain CDS sequences to build codon bias models"

# Invalid score identifier
try:
    feature = CodonBiasFeature.from_reference(ref_set, ["INVALID"])
except ValueError as e:
    print(e)
    # "Cannot resolve score identifier 'INVALID'. Expected one of: CAI, ENC, ..."

# CDS not multiple of 3
try:
    bad_ref = ReferenceSequenceSet(cds={"gene1": "ATGAA"})  # 5 nt, not divisible by 3
    cds = bad_ref.cds_strings(require_multiple_of_three=True)
except ValueError as e:
    print(e)
    # "CDS sequence 'gene1' has length 5 which is not a multiple of 3"

# Internal stop codons
try:
    bad_ref = ReferenceSequenceSet(cds={"gene1": "ATGTAAAAATAA"})  # Internal stop
    proteins = bad_ref.protein_strings(error_on_internal_stop=True)
except ValueError as e:
    print(e)
    # "CDS sequence 'gene1' (index: 0) contains internal stop codon at protein position 1"
```

## Caching and Performance

The reference set caches normalized sequences for efficiency:

```python
ref_set = ReferenceSequenceSet.from_fasta("genes.fasta")

# First call: normalizes and caches
cds1 = ref_set.cds_strings()  # Normalization performed

# Subsequent calls: returns cached result
cds2 = ref_set.cds_strings()  # Instant (uses cache)

# Same for proteins
proteins1 = ref_set.protein_strings()  # Translation performed and cached
proteins2 = ref_set.protein_strings()  # Instant (uses cache)
```

The cache is automatically validated when stricter requirements are requested:

```python
# First call without validation
cds = ref_set.cds_strings(require_multiple_of_three=False)

# Second call with validation: checks cached sequences
cds_validated = ref_set.cds_strings(require_multiple_of_three=True)
```

## See Also

- [Codon Bias Features](./codon_bias.md) for detailed codon bias usage and windowing
- [ORF Handling](./orfs.md) for finding and selecting ORFs
- [ORF Windows](./orf_windows.md) for windowed feature computation concepts
- [NCBI Genetic Codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi) for translation table reference
- [codon-bias package](https://pypi.org/project/codon-bias/) for underlying codon bias algorithms
