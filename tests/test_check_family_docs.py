"""Tests for family documentation checker."""

import tempfile
from pathlib import Path

from scripts.check_family_docs import check_readme


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

## Maintenance notes
Content
"""
            readme_path.write_text(content)
            errors = check_readme(readme_path, "test_family")
            # Should detect missing "## Examples" even though "## Examples_REMOVED" exists
            assert len(errors) == 1
            assert "Examples" in errors[0]
