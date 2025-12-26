# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]
### Added
- Release automation via GitHub Actions for building and publishing to PyPI.
- Packaging metadata improvements (setuptools `src` layout, license, classifiers).
- Release documentation, including steps for tagging and publishing.
- Packaging polish: explicit `LICENSE` metadata, project URLs, and documented optional extras.
### Fixed
- Corrected `project.license` to use an SPDX-compatible string for PyPI validation.
- Removed legacy license classifier in favor of the license expression required by setuptools/PEP 639.

## [0.1.0]
### Added
- Initial public release with core bioinformatics utilities, feature families, and windowing utilities.
