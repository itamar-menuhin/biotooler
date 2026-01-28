# Frequently Asked Questions

## Installation & Setup

### Q: Which Python versions are supported?

Python 3.11 and 3.12 are currently supported. Python 3.10 and earlier are not supported due to the use of newer type hinting features.

### Q: Why am I getting ImportError for a feature family?

Feature families with optional dependencies are not installed by default. Install the specific extra you need:

```bash
pip install "biotooler[codon_bias]"  # For codon bias features
pip install "biotooler[disorder]"     # For disorder prediction
```

See the [Installation Guide](installation.md) for all available extras.

### Q: Can I install all features at once?

Yes, but be aware this installs many heavy dependencies:

```bash
pip install "biotooler[codon_bias,disorder,viennarna,chimera]"
```

We recommend installing only what you need to keep your environment lean.

## Usage

### Q: What sequence formats are supported?

biotooler works with BioPython `SeqRecord` objects. You can load sequences from any format BioPython supports:

```python
from Bio import SeqIO

# From FASTA
records = SeqIO.parse("sequences.fasta", "fasta")

# From GenBank
records = SeqIO.parse("sequences.gb", "genbank")

# Then use with biotooler
from biotooler.families.basic_stats import BasicStatsFeature
feature = BasicStatsFeature()
for record in records:
    result = feature(record)
```

### Q: Can I use biotooler with DNA, RNA, and protein sequences?

Yes! Most features automatically detect sequence type. Some features are specific to certain types:

- **DNA/RNA only**: Codon bias features, ViennaRNA
- **Protein only**: Disorder prediction, some ProtParam features
- **All types**: Basic stats (adjusts output based on type)

### Q: What's the difference between `compute_windows` and `compute_orf_windows`?

- **`compute_windows`**: Generic sliding windows, works with any sequence type, any step size
- **`compute_orf_windows`**: Codon-aligned windows, DNA/RNA only, step must be multiple of 3

```python
# Generic windows - any step size
fs.compute_windows(record, window_size=10, step=1)

# ORF windows - codon-aligned
fs.compute_orf_windows(record, orf=(0, 99), window_nt=9, step_nt=3)
```

### Q: How do I handle sequences with ambiguous nucleotides?

biotooler handles standard ambiguous nucleotides:

- **N** in DNA/RNA: Counted separately, excluded from GC content
- **X** in proteins: Counted separately in amino acid composition

Other IUPAC ambiguity codes (R, Y, S, W, K, M) are not explicitly handled by most features.

### Q: Can I compute features on subsequences?

Yes, use the `region` parameter:

```python
# Compute only on positions 100-200
result = fs.compute_windows(
    record,
    window_size=10,
    step=5,
    region=(100, 200)
)
```

## Performance

### Q: How can I speed up sliding window analysis?

1. Use **incremental features** when available (like BasicStatsFeature)
2. Increase **step size** to reduce overlapping windows
3. Use **larger windows** when appropriate
4. Process sequences in **batches**

```python
# Slow: step=1 with small window
fs.compute_windows(record, window_size=10, step=1)  # Many windows

# Faster: larger step
fs.compute_windows(record, window_size=10, step=10)  # Fewer windows
```

### Q: What's incremental computation?

Some features support incremental computation, which reuses calculations from overlapping windows:

- **BasicStatsFeature**: Incremental ✅
- **CodonBiasFeature**: Rolling mode ✅
- **DisorderProfileMetapredict**: Vectorized ✅

This can be 10-100x faster for highly overlapping windows.

### Q: How much memory does windowing use?

Memory usage depends on:

- Number of windows (determined by sequence length, window size, step)
- Number of features computed
- Feature type (scalar features are lightweight, vectors are heavier)

For very large sequences, consider processing in chunks.

## Feature Families

### Q: Which codon bias metric should I use?

Depends on your question:

- **CAI**: Codon Adaptation Index - overall adaptation to reference set
- **ENC**: Effective Number of Codons - overall bias magnitude
- **FOP**: Frequency of Optimal Codons - bias toward specific "optimal" codons
- **tAI**: tRNA Adaptation Index - considers tRNA availability

See the [Codon Bias documentation](../families/codon_bias.md) for details.

### Q: How accurate is disorder prediction?

Metapredict achieves ~80% accuracy on DisProt dataset. Accuracy varies by protein type and disorder definition. For critical applications, validate predictions experimentally or with multiple predictors.

### Q: Can I use my own reference sequences for codon bias?

Yes! Provide your own reference sequences:

```python
from codonbias.scores import CodonAdaptationIndex
from biotooler.families.codon_bias import CodonBiasFeature

# Use your highly expressed genes as reference
ref_sequences = "ATGATGATG..."  # Your reference genes
cai = CodonAdaptationIndex(ref_sequences)
feature = CodonBiasFeature([cai], names=["CAI_custom"])
```

## Development

### Q: How do I add a new feature family?

Use the scaffolder script:

```bash
python scripts/new_family.py --name my_feature --owner @your_github_username
```

Then follow the [Adding a Family guide](../dev/adding_a_family.md).

### Q: Can I contribute to biotooler?

Yes! Contributions are welcome. See [CONTRIBUTING.md](https://github.com/itamar-menuhin/biotooler/blob/main/CONTRIBUTING.md) for guidelines.

### Q: How do I run tests?

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run all tests
pytest

# Run specific test file
pytest tests/families/test_basic_stats.py

# Run with coverage
pytest --cov=biotooler
```

### Q: What's the development workflow?

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Run tests, linting, and type checking
5. Open a pull request

```bash
# Check code quality
ruff check .
ruff format .
pyright

# Run tests
pytest
```

## Troubleshooting

### Q: I'm getting "ORF index out of range"

This means you're trying to access an ORF that doesn't exist:

```python
candidates = find_orf_candidates(record)
# candidates = [(0, 9), (9, 18)]  # Only 2 ORFs

# This will fail:
orf = select_orf_by_index(candidates, 5)  # Index 5 doesn't exist

# Use valid index:
orf = select_orf_by_index(candidates, 0)  # OK
```

### Q: Why is my window computation slow?

- Check if you're using step=1 with large sequences
- Verify you're using incremental features if available
- Consider if you need all windows (maybe increase step size)
- Profile your code to identify bottlenecks

### Q: How do I report a bug?

Open an issue on GitHub with:

1. biotooler version (`import biotooler; print(biotooler.__version__)`)
2. Python version
3. Minimal reproducible example
4. Expected vs actual behavior

Use the bug report template for best results.

## Still Have Questions?

- Check the [User Guide](../usage/windows.md)
- Browse [API Reference](../api/core.md)
- Open a [GitHub Discussion](https://github.com/itamar-menuhin/biotooler/discussions)
- File an [Issue](https://github.com/itamar-menuhin/biotooler/issues)