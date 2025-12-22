"""Tests for lazy_import helper."""

import pytest

from biotooler.core.lazy_import import lazy_import


class TestLazyImport:
    """Tests for lazy_import function."""

    def test_lazy_import_existing_returns_module(self):
        """Test that lazy_import returns a module for an existing module."""
        # Import a standard library module that definitely exists
        module = lazy_import("json", extra="test", purpose="testing")
        assert module.__name__ == "json"
        # Verify it's a module type
        assert hasattr(module, "loads")
        assert hasattr(module, "dumps")

    def test_lazy_import_missing_raises_with_hint(self):
        """Test that lazy_import raises ImportError with install hint for missing module."""
        with pytest.raises(ImportError) as exc_info:
            lazy_import(
                "nonexistent_module_xyz123",
                extra="test_extra",
                purpose="testing lazy imports",
            )

        # Check the error message contains the required information
        error_message = str(exc_info.value)
        assert "nonexistent_module_xyz123" in error_message
        assert "testing lazy imports" in error_message
        assert 'pip install "biotooler[test_extra]"' in error_message

    def test_lazy_import_preserves_original_error(self):
        """Test that lazy_import chains the original ImportError."""
        with pytest.raises(ImportError) as exc_info:
            lazy_import(
                "nonexistent_module_xyz456",
                extra="test",
                purpose="testing",
            )

        # Check that the original error is chained
        assert exc_info.value.__cause__ is not None
        assert isinstance(exc_info.value.__cause__, ImportError)
