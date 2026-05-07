"""Smoke test for the util_legacy thin shim.

Verifies that NotebookUtil can be imported and instantiated from util_legacy,
confirming the shim correctly delegates to KBUtilLib's mixin classes.
"""
import sys
import os

# Ensure we're testing the local module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


def test_import_and_instantiate():
    """NotebookUtil should import and instantiate without errors."""
    from util_legacy import NotebookUtil

    instance = NotebookUtil()
    assert instance is not None
    # Verify it has the expected custom methods
    assert hasattr(instance, "add_dgoa_reaction")
    assert hasattr(instance, "reaction_equation_with_names")
    assert hasattr(instance, "get_model_and_simulate")
    assert hasattr(instance, "find_significant_differences")
    assert hasattr(instance, "translate_expression_gene_ids")
    assert hasattr(instance, "classify_fva_flux")
    assert hasattr(instance, "build_gene_reaction_map")
    # Verify inherited KBUtilLib methods are accessible
    assert hasattr(instance, "get_object")
    assert hasattr(instance, "get_media")
    assert hasattr(instance, "save")
    assert hasattr(instance, "load")
    print("PASS: NotebookUtil imported and instantiated successfully from util_legacy shim")


def test_module_level_util():
    """The module-level `util` instance should exist."""
    from util_legacy import util

    assert util is not None
    assert hasattr(util, "add_dgoa_reaction")
    print("PASS: Module-level util instance exists")


if __name__ == "__main__":
    test_import_and_instantiate()
    test_module_level_util()
    print("\nAll util_legacy shim tests passed.")
