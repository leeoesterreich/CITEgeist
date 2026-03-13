"""Ensure CITEgeist/CITEgeist/ is on sys.path for model imports.

The repo has two model/ dirs: an outer legacy one (CITEgeist/model/) and the
active development one (CITEgeist/CITEgeist/model/). When pytest runs from the
repo root, the outer model/ gets imported first. This conftest forces the inner
model package to be used instead.
"""
import os
import sys

_inner_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, _inner_dir)

# Force re-resolution of the model package to the inner dir
if 'model' in sys.modules:
    del sys.modules['model']
