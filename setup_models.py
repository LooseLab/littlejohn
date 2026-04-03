#!/usr/bin/env python3
"""
Backward-compatible entry point for downloading ROBIN model assets.

Preferred (after ``pip install -e .``)::

    robin utils update-models
    robin utils update-clinvar
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
_SRC = _REPO_ROOT / "src"
if _SRC.is_dir():
    sys.path.insert(0, str(_SRC))


def main() -> None:
    msg = "setup_models.py is deprecated; use: robin utils update-models && robin utils update-clinvar"
    warnings.warn(msg, DeprecationWarning, stacklevel=1)
    print(msg, file=sys.stderr)
    from robin.utils.clinvar_manager import update_clinvar_if_newer
    from robin.utils.model_checker import get_models_directory
    from robin.utils.model_updater import update_models

    models_dir = get_models_directory(_REPO_ROOT)
    print(f"Models directory: {models_dir}")
    github_token = os.getenv("GITHUB_TOKEN")
    if not github_token:
        print(
            "Note: GITHUB_TOKEN is not set; private GitHub assets may fail to download.",
        )

    ok, msgs = update_models(models_dir=models_dir, github_token=github_token)
    for m in msgs:
        print(m)
    if not ok:
        sys.exit(1)

    try:
        updated = update_clinvar_if_newer(download_if_missing=True)
        print("ClinVar: updated." if updated else "ClinVar: already up to date or present.")
    except Exception as e:
        print(f"ClinVar setup failed: {e}", file=sys.stderr)
        print("Retry with: robin utils update-clinvar", file=sys.stderr)
        sys.exit(1)

    print("\nModel and ClinVar assets are ready.")


if __name__ == "__main__":
    main()
