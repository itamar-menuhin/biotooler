#!/usr/bin/env python3
"""Scaffolder script for creating new feature families in biotooler.

This script creates a new feature family from the template skeleton, replacing
placeholders with the specified family name and owner information. It also
creates a smoke test and updates the CODEOWNERS file.

Usage:
    python scripts/new_family.py --name my_feature --owner @username
    python scripts/new_family.py --name my_feature --owner @username \\
        --extra "Additional description"
"""

import argparse
import re
from pathlib import Path


def validate_family_name(name: str) -> str:
    """Validate that the family name is a valid Python module name.

    Args:
        name: The proposed family name

    Returns:
        The validated name (lowercase, with underscores)

    Raises:
        ValueError: If the name is not a valid Python module name
    """
    # Convert to lowercase and replace hyphens with underscores
    name = name.lower().replace("-", "_")

    # Check that it's a valid Python identifier
    if not re.match(r"^[a-z][a-z0-9_]*$", name):
        raise ValueError(
            f"Invalid family name: {name}. "
            "Must start with a letter and contain only lowercase letters, digits, and underscores."
        )

    # Check that it's not a Python keyword
    import keyword

    if keyword.iskeyword(name):
        raise ValueError(f"Invalid family name: {name}. Cannot be a Python keyword.")

    return name


def to_title_case(name: str) -> str:
    """Convert snake_case name to TitleCase.

    Args:
        name: Snake case name (e.g., 'my_feature')

    Returns:
        Title case name (e.g., 'MyFeature')
    """
    return "".join(word.capitalize() for word in name.split("_"))


def replace_placeholders(content: str, replacements: dict[str, str]) -> str:
    """Replace placeholders in content with actual values.

    Args:
        content: File content with placeholders
        replacements: Dictionary mapping placeholder names to replacement values

    Returns:
        Content with placeholders replaced
    """
    for placeholder, value in replacements.items():
        content = content.replace(f"{{{{{placeholder}}}}}", value)
    return content


def create_family(name: str, owner: str, extra: str = "") -> None:
    """Create a new feature family from the template.

    Args:
        name: Family name (will be converted to valid Python module name)
        owner: GitHub username or team for CODEOWNERS (e.g., '@username')
        extra: Additional description text (optional)
    """
    # Get repository root (assuming script is in scripts/)
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    # Validate and normalize the family name
    family_name = validate_family_name(name)
    family_name_title = to_title_case(family_name)

    print(f"Creating new feature family: {family_name}")
    print(f"  Title case: {family_name_title}")
    print(f"  Owner: {owner}")

    # Set up paths
    template_dir = repo_root / "templates" / "family_skeleton"
    family_dir = repo_root / "src" / "biotooler" / "families" / family_name
    test_dir = repo_root / "tests" / "families"
    codeowners_path = repo_root / ".github" / "CODEOWNERS"

    # Check if family already exists
    if family_dir.exists():
        raise ValueError(
            f"Family '{family_name}' already exists at {family_dir}. "
            "Please choose a different name or remove the existing family first."
        )

    # Prepare placeholder replacements
    description = extra if extra else f"Feature family for {family_name} analysis"
    replacements = {
        "FAMILY_NAME": family_name,
        "FAMILY_NAME_TITLE": family_name_title,
        "FAMILY_DESCRIPTION": description,
    }

    # Create family directory
    print(f"\nCreating family directory: {family_dir}")
    family_dir.mkdir(parents=True, exist_ok=True)

    # Copy and process template files
    for template_file in template_dir.iterdir():
        if template_file.is_file() and not template_file.name.startswith("."):
            # Determine target filename
            if template_file.name == "feature.py":
                target_name = f"{family_name}.py"
            else:
                target_name = template_file.name

            target_path = family_dir / target_name

            # Read template content
            content = template_file.read_text()

            # Replace placeholders
            content = replace_placeholders(content, replacements)

            # Write to target
            target_path.write_text(content)
            print(f"  Created: {target_path.relative_to(repo_root)}")

    # Create test directory if it doesn't exist
    test_dir.mkdir(parents=True, exist_ok=True)

    # Create smoke test
    test_file = test_dir / f"test_{family_name}_smoke.py"
    test_content = f'''"""Smoke tests for {family_name} feature family."""

from pathlib import Path


def test_{family_name}_imports():
    """Test that {family_name} family can be imported."""
    from biotooler.families.{family_name} import {family_name_title}Feature

    assert {family_name_title}Feature is not None


def test_{family_name}_readme_exists():
    """Test that README.md exists for {family_name} family."""
    # Get the family directory path
    import biotooler.families.{family_name}

    family_module_path = Path(biotooler.families.{family_name}.__file__).parent
    readme_path = family_module_path / "README.md"

    assert readme_path.exists(), f"README.md not found at {{readme_path}}"
    assert readme_path.stat().st_size > 0, "README.md is empty"
'''

    test_file.write_text(test_content)
    print(f"  Created: {test_file.relative_to(repo_root)}")

    # Update CODEOWNERS file
    print(f"\nUpdating CODEOWNERS: {codeowners_path.relative_to(repo_root)}")

    # Ensure owner starts with @ or is a valid team reference
    if not owner.startswith("@"):
        print(
            f"  Warning: Owner '{owner}' doesn't start with '@'. "
            "Adding '@' prefix for GitHub username format."
        )
        owner = f"@{owner}"

    # Create CODEOWNERS entry
    codeowners_entry = f"/src/biotooler/families/{family_name}/ {owner}\n"

    if codeowners_path.exists():
        # Append to existing file
        with codeowners_path.open("a") as f:
            f.write(codeowners_entry)
        print(f"  Appended entry for {family_name}")
    else:
        # Create new CODEOWNERS file
        codeowners_content = f"""# CODEOWNERS file for biotooler feature families
#
# Each feature family has its own directory and designated owners
# This enables isolated contributions and clear ownership

{codeowners_entry}"""
        codeowners_path.write_text(codeowners_content)
        print(f"  Created new CODEOWNERS file with entry for {family_name}")

    # Print success message and next steps
    print(f"\n✓ Successfully created feature family '{family_name}'!")
    print("\nNext steps:")
    family_rel = family_dir.relative_to(repo_root)
    print(f"  1. Edit {family_rel}/{family_name}.py to implement your feature")
    print(f"  2. Update {family_rel}/README.md with documentation")
    print(f"  3. Add tests in tests/families/test_{family_name}.py")
    print(
        "  4. (Optional) Add lazy import to src/biotooler/features/__init__.py"
        f" following instructions in {family_rel}/integration.py"
    )
    print("  5. Run tests to verify your implementation")
    print("\nSee CONTRIBUTING.md for more details on contributing new features.")


def main():
    """Main entry point for the scaffolder script."""
    parser = argparse.ArgumentParser(
        description="Create a new feature family from template",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/new_family.py --name gc_content --owner @john_doe
  python scripts/new_family.py --name motif_finder --owner @team/bioinformatics \\
    --extra "TFBS finder"
        """,
    )

    parser.add_argument(
        "--name",
        required=True,
        help="Name of the feature family (will be converted to snake_case)",
    )

    parser.add_argument(
        "--owner",
        required=True,
        help="GitHub username or team for CODEOWNERS (e.g., @username or @org/team)",
    )

    parser.add_argument(
        "--extra",
        default="",
        help="Additional description text for the feature family (optional)",
    )

    args = parser.parse_args()

    try:
        create_family(args.name, args.owner, args.extra)
    except Exception as e:
        print(f"\n✗ Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
