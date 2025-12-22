# Contributing to biotooler

Thank you for your interest in contributing to biotooler! This document provides guidelines for contributing new feature families and other improvements.

## Contributing a New Feature Family

Biotooler organizes features into isolated family submodules under `src/biotooler/families/<family>/`. This structure enables:
- Clear ownership and maintenance responsibilities via CODEOWNERS
- Isolated dependencies (heavy dependencies are only loaded when explicitly used)
- Easy contribution workflow for external contributors

### Quick Start: Using the Family Scaffolder

We provide a scaffolder script that creates a new feature family with the correct structure:

```bash
python scripts/new_family.py --name my_feature --owner @your_github_username
```

With optional description:
```bash
python scripts/new_family.py --name my_feature --owner @your_github_username --extra "Brief description of the feature"
```

This will create:
- `src/biotooler/families/my_feature/` - Your feature family directory
  - `__init__.py` - Package initialization
  - `my_feature.py` - Feature implementation skeleton
  - `integration.py` - Instructions for lazy import integration
  - `README.md` - Documentation template with required sections
- `tests/families/test_my_feature_smoke.py` - Basic smoke tests
- Updated `.github/CODEOWNERS` - Adds you as the owner of your feature family

### Implementation Steps

After running the scaffolder:

1. **Implement your feature** in `src/biotooler/families/<name>/<name>.py`:
   - Fill in the `__call__` method for basic feature computation
   - Optionally implement `init_state`, `step_state`, and `emit` for incremental computation over sliding windows
   - Add appropriate error handling and validation

2. **Complete the README** in `src/biotooler/families/<name>/README.md`:
   - Replace all TODO sections with actual content
   - Include clear usage examples
   - Document the API thoroughly
   - Add references to papers or algorithms if applicable

3. **Write tests** in `tests/families/test_<name>.py`:
   - Add comprehensive unit tests for your feature
   - Test edge cases and error conditions
   - Include integration tests if applicable
   - The scaffolder creates a basic smoke test to get you started

4. **(Optional) Add lazy import** to `src/biotooler/features/__init__.py`:
   - Follow the instructions in your family's `integration.py` file
   - This makes your feature accessible via `from biotooler.features import YourFeature`
   - Ensures heavy dependencies are only loaded when explicitly used

5. **Test your implementation**:
   ```bash
   # Run all tests
   python -m pytest tests/
   
   # Run only your family's tests
   python -m pytest tests/families/test_<name>.py
   
   # Check code style
   ruff check src/ tests/
   ruff format src/ tests/
   
   # Type checking
   pyright src/
   ```

### Feature Interface

Features should implement one or both of these interfaces:

#### Basic Feature Interface
```python
def __call__(self, record: SeqRecord) -> dict[str, Scalar]:
    """Compute features for the entire sequence."""
    pass
```

#### Incremental Feature Interface (Optional)
For better performance with sliding windows:
```python
def init_state(self, record, *, orf, window_start, window_end, **kwargs) -> dict:
    """Initialize state for the first window."""
    pass

def step_state(self, state, *, out_start, out_end, in_start, in_end, **kwargs) -> None:
    """Update state for the next window."""
    pass

def emit(self, state) -> dict[str, Scalar]:
    """Emit feature values from current state."""
    pass
```

See existing families (e.g., `codon_bias`) for complete examples.

### Code Style and Quality

- Follow PEP 8 style guidelines
- Use type hints for all function signatures
- Write clear docstrings for all public APIs
- Keep imports lightweight (use lazy imports for heavy dependencies)
- Run linters before submitting: `ruff check` and `ruff format`
- Ensure type checking passes: `pyright`

### Minimal Changes Philosophy

When contributing:
- Make the smallest possible changes to achieve your goal
- Don't refactor existing code unless necessary for your feature
- Don't modify other feature families
- Keep your feature isolated in its own family directory

### Documentation Requirements

Your family's README.md must include these sections:
- **Overview**: What the feature does and why it's useful
- **Features**: List of specific features provided
- **Installation**: Any optional dependencies needed
- **Usage**: Complete working examples
- **API Reference**: Detailed documentation of the public API
- **Implementation Details**: Algorithms, performance considerations
- **Contributing**: Family-specific guidelines
- **References**: Papers, algorithms, external packages
- **License**: License information

### CODEOWNERS

When you create a new family, you're automatically added as the owner in `.github/CODEOWNERS`. This means:
- You'll be automatically requested for review on PRs affecting your family
- You have primary responsibility for maintaining your feature family
- You can approve changes to your family's code

### Questions?

If you have questions or need help:
- Open an issue on GitHub
- Check existing feature families for examples
- Review the biotooler documentation

Thank you for contributing to biotooler!
