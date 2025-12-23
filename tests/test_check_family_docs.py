"""Tests for family documentation checker."""

import tempfile
from pathlib import Path

from scripts.check_family_docs import check_readme


# Helper to generate complete windowing correctness section
# For families WITHOUT upstream libraries
WINDOWING_SECTION_NO_UPSTREAM = """
## Windowing correctness

### Position space
This family uses RESIDUE position space.

### Vector computation
The compute_vector method processes the full sequence.

### Aggregation strategy
We use np.mean for averaging.

### Testing approach
Tests validate full-context computation.
"""

# For families WITH upstream libraries  
WINDOWING_SECTION_WITH_UPSTREAM = """
## Windowing correctness

### Position space
This family uses RESIDUE position space.

### Vector computation
The compute_vector method calls the upstream library API with the full sequence.

### Aggregation strategy
We use np.mean for averaging.

### Testing approach
Tests validate full-context computation.
"""


class TestCheckFamilyDocs:
    """Tests for check_family_docs script."""

    def test_missing_readme(self):
        """Test that missing README is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 1
            assert "Missing README.md" in errors[0]

    def test_missing_heading(self):
        """Test that missing required heading is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
- Link: https://example.com

## Examples
Content

## Edge cases and validation
Content

## Windowing correctness

### Position space
Content

### Vector computation
The upstream library API is used to compute vectors.

### Aggregation strategy
Content

### Testing approach
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Missing "## Maintenance notes"
            assert len(errors) == 1
            assert "Maintenance notes" in errors[0]

    def test_references_without_link(self):
        """Test that References section without http link is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
Some text without a link

## Upstream library links
- Link: https://example.com

## Examples
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_WITH_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 1
            assert "References" in errors[0]
            assert "http(s) link" in errors[0]

    def test_upstream_without_link_or_na(self):
        """Test that Upstream library links without link or N/A is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
Some text without a link

## Examples
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_NO_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 1
            assert "Upstream library links" in errors[0]

    def test_upstream_with_na(self):
        """Test that N/A in Upstream library links is accepted."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
N/A - This family has no upstream library

## Examples
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_NO_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 0

    def test_complete_readme(self):
        """Test that complete README passes validation."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Paper: https://example.com/paper

## Upstream library links
- Library: https://example.com/library

## Examples
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_WITH_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 0

    def test_heading_substring_not_matched(self):
        """Test that heading substrings don't match (e.g., '## Examples_REMOVED')."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
- Link: https://example.com

## Examples_REMOVED
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_WITH_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Should detect missing "## Examples" even though "## Examples_REMOVED" exists
            assert len(errors) == 1
            assert "Examples" in errors[0]

    def test_upstream_na_word_boundary(self):
        """Test that N/A detection uses word boundaries (e.g., 'BANANA' should not match)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
This family uses BANANA library

## Examples
Content

## Edge cases and validation
Content

""" + WINDOWING_SECTION_NO_UPSTREAM + """

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Should fail because 'BANANA' contains 'NA' but not as a word
            assert len(errors) == 1
            assert "Upstream library links" in errors[0]

    def test_missing_windowing_correctness(self):
        """Test that missing Windowing correctness section is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
- Link: https://example.com

## Examples
Content

## Edge cases and validation
Content

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Missing "## Windowing correctness"
            assert len(errors) == 1
            assert "Windowing correctness" in errors[0]

    def test_windowing_correctness_missing_subsections(self):
        """Test that Windowing correctness section must have required subsections."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
N/A - no upstream library

## Examples
Content

## Edge cases and validation
Content

## Windowing correctness
This section exists but is incomplete.

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Should have 4 errors - one for each missing subsection
            assert len(errors) == 4
            assert any("Position space" in e for e in errors)
            assert any("Vector computation" in e for e in errors)
            assert any("Aggregation strategy" in e for e in errors)
            assert any("Testing approach" in e for e in errors)

    def test_windowing_correctness_with_all_subsections(self):
        """Test that complete Windowing correctness section passes."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
N/A - no upstream library

## Examples
Content

## Edge cases and validation
Content

## Windowing correctness

### Position space
This family uses RESIDUE position space.

### Vector computation
The compute_vector method processes the full sequence.

### Aggregation strategy
We use np.mean for averaging.

### Testing approach
Tests validate full-context computation.

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 0

    def test_windowing_correctness_with_upstream_library(self):
        """Test that windowing section mentions upstream library when appropriate."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
- MyLibrary: https://github.com/example/mylibrary

## Examples
Content

## Edge cases and validation
Content

## Windowing correctness

### Position space
This family uses RESIDUE position space.

### Vector computation
The compute_vector method processes the full sequence without mentioning how.

### Aggregation strategy
We use np.mean for averaging.

### Testing approach
Tests validate full-context computation.

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Should have 1 error - Vector computation doesn't mention upstream library
            assert len(errors) == 1
            assert "Vector computation" in errors[0]
            assert "upstream" in errors[0].lower() or "library" in errors[0].lower()

    def test_windowing_correctness_with_upstream_library_properly_documented(self):
        """Test that windowing section properly documents upstream library usage."""
        with tempfile.TemporaryDirectory() as tmpdir:
            readme_path = Path(tmpdir) / "README.md"
            content = """
# Test Family

## What this family provides
Content

## Intuition
Content

## Mathematical formulation
Content

## Features and output schema
Content

## References
- Link: https://example.com

## Upstream library links
- MyLibrary: https://github.com/example/mylibrary

## Examples
Content

## Edge cases and validation
Content

## Windowing correctness

### Position space
This family uses RESIDUE position space.

### Vector computation
The compute_vector method calls the upstream library API with the full sequence.

### Aggregation strategy
We use np.mean for averaging.

### Testing approach
Tests validate full-context computation.

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            assert len(errors) == 0
