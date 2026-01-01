# Releasing biotooler

The repository now includes automation for publishing releases to PyPI. Follow these steps to cut a release:

1. **Update version metadata**
   - Bump `version` in `pyproject.toml`.
   - Mirror the version in `src/biotooler/__init__.py` so `__version__` matches the release tag.
   - Add a new section to [CHANGELOG.md](../../CHANGELOG.md) describing the release.
2. **Verify changes**
   - Run `pytest`, `ruff check .`, and `pyright` locally.
   - Ensure any optional dependencies required for tests are available (e.g., `codon-bias`).
3. **Create a GitHub release**
   - Push a git tag that matches the version (for example, `v0.1.1`).
   - Draft and publish a GitHub release for that tag. Publishing triggers the `Publish` workflow.
   - Copy the relevant [CHANGELOG](../../CHANGELOG.md) entry into the release notes.
4. **Configure secrets**
   - Add a `PYPI_API_TOKEN` repository secret with a valid PyPI token for the `biotooler` project.
5. **Check workflow output**
   - The workflow builds source and wheel distributions using `python -m build`.
   - Build artifacts are attached to the workflow run as `dist-artifacts`.
   - The final step publishes the release to PyPI using the configured token.

If you want to test the pipeline without publishing, trigger the workflow manually via **Actions → Publish → Run workflow**, then cancel before the publish step or temporarily unset the `PYPI_API_TOKEN` secret.
