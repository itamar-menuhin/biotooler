#!/usr/bin/env python3
"""Check that each family has a complete README with required sections."""

import re
import sys
from pathlib import Path

REQUIRED_HEADINGS = [
    "## What this family provides",
    "## Intuition",
    "## Mathematical formulation",
    "## Features and output schema",
    "## References",
    "## Upstream library links",
    "## Examples",
    "## Edge cases and validation",
    "## Windowing correctness",
    "## Maintenance notes",
]


def _extract_section(content: str, heading: str) -> str | None:
    """Extract content between a heading and the next heading or end of file."""
    match = re.search(
        rf"^{re.escape(heading)}(?!#)\s*(.*?)(?=^##(?!#)|\Z)", content, re.DOTALL | re.MULTILINE
    )
    return match.group(1) if match else None


def check_readme(readme_path: Path, family_name: str) -> list[str]:
    """Check a single README for completeness. Returns list of errors."""
    errors = []

    if not readme_path.exists():
        errors.append(f"Family '{family_name}': Missing README.md")
        return errors

    content = readme_path.read_text()

    # Check for required headings (with word boundary to ensure exact match)
    for heading in REQUIRED_HEADINGS:
        # Use regex to match heading at start of line with optional whitespace after
        pattern = r"^" + re.escape(heading) + r"(?:\s|$)"
        if not re.search(pattern, content, re.MULTILINE):
            errors.append(f"Family '{family_name}': Missing heading '{heading}'")

    # Check References section has at least one http(s) link
    references_section = _extract_section(content, "## References")
    if references_section and not re.search(r"https?://", references_section):
        errors.append(
            f"Family '{family_name}': 'References' section must contain "
            "at least one http(s) link"
        )

    # Check Upstream library links section has at least one http(s) link or "N/A"
    upstream_section = _extract_section(content, "## Upstream library links")
    if upstream_section:
        has_link = re.search(r"https?://", upstream_section)
        has_na = re.search(r"\bN/A\b", upstream_section)
        if not has_link and not has_na:
            errors.append(
                f"Family '{family_name}': 'Upstream library links' section must contain "
                "at least one http(s) link or 'N/A'"
            )

    # Check Windowing correctness section has required subsections
    windowing_section = _extract_section(content, "## Windowing correctness")
    if windowing_section:
        required_subsections = [
            "### Position space",
            "### Vector computation",
            "### Aggregation strategy",
            "### Testing approach",
        ]
        for subsection in required_subsections:
            if subsection not in windowing_section:
                errors.append(
                    f"Family '{family_name}': 'Windowing correctness' section must "
                    f"include '{subsection}' subsection"
                )
        
        # Check that windowing section references upstream APIs when appropriate
        # If Upstream library links has actual links (not N/A), windowing should mention them
        if upstream_section:
            has_upstream_link = re.search(r"https?://", upstream_section)
            has_na = re.search(r"\bN/A\b", upstream_section)
            if has_upstream_link and not has_na:
                # Family uses upstream library - check windowing docs reference it
                vector_computation_match = re.search(
                    r"### Vector computation\s+(.*?)(?=###|\Z)", 
                    windowing_section, 
                    re.DOTALL
                )
                if vector_computation_match:
                    vector_computation_text = vector_computation_match.group(1)
                    # Should mention library, API, or link to documentation
                    has_library_ref = (
                        re.search(r"\blibrary\b", vector_computation_text, re.IGNORECASE) or
                        re.search(r"\bAPI\b", vector_computation_text) or
                        re.search(r"https?://", vector_computation_text) or
                        re.search(r"\bupstream\b", vector_computation_text, re.IGNORECASE)
                    )
                    if not has_library_ref:
                        errors.append(
                            f"Family '{family_name}': 'Windowing correctness > Vector computation' "
                            "section should explain how upstream library APIs are used "
                            "(mention 'library', 'API', 'upstream', or include documentation links)"
                        )

    return errors


def main() -> int:
    """Main entry point. Returns exit code."""
    repo_root = Path(__file__).parent.parent
    families_dir = repo_root / "src" / "biotooler" / "families"

    if not families_dir.exists():
        print(f"Error: Families directory not found at {families_dir}")
        return 1

    all_errors = []

    # Find all family directories (exclude __pycache__ and similar)
    family_dirs = [
        d for d in families_dir.iterdir()
        if d.is_dir() and not d.name.startswith("_") and not d.name.startswith(".")
    ]

    if not family_dirs:
        print("Warning: No family directories found")
        return 0

    for family_dir in sorted(family_dirs):
        family_name = family_dir.name
        readme_path = family_dir / "README.md"
        errors = check_readme(readme_path, family_name)
        all_errors.extend(errors)

    if all_errors:
        print("Family documentation validation failed:\n")
        for error in all_errors:
            print(f"  ❌ {error}")
        print(f"\nTotal errors: {len(all_errors)}")
        return 1

    print("✅ All family documentation is complete")
    return 0


if __name__ == "__main__":
    sys.exit(main())
