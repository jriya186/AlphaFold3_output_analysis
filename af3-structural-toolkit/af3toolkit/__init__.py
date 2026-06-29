"""
af3toolkit
==========

Post-AlphaFold3 structural analysis utilities for LDL-receptor-family /
Reelin complexes (ApoER2, LDLR, VLDLR).

Core idea: AF3 constructs are usually spliced fragments (not full-length
proteins), so raw PDB residue numbering does NOT correspond to the real
position in the full-length protein. This package re-maps local PDB
positions back to full-sequence positions, exon/repeat identity, and
(for receptors) structural domain -- while keeping receptor-side and
Reelin-side confidence levels clearly distinguished.
"""

from .pipeline import analyze_complex
from .confidence import get_best_models_by_iptm_ptm

__all__ = ["analyze_complex", "get_best_models_by_iptm_ptm"]

__version__ = "0.1.0"
