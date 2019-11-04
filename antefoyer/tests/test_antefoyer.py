"""
Unit and regression test for the antefoyer package.
"""

# Import package, test suite, and other packages as needed
import antefoyer
import pytest
import sys

def test_antefoyer_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "antefoyer" in sys.modules
