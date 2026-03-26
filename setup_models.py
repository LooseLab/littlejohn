#!/usr/bin/env python3
"""
Setup script to download all required model assets for Robin.
This script should be run after cloning the repository to get the model files.
"""

import os
import sys
from pathlib import Path

# Add the project root to the Python path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from robin.utils.model_checker import get_required_models
from robin.utils.model_updater import update_models


def main():
    """Download all required model assets"""
    models_dir = project_root / "src" / "robin" / "models"
    
    print("Setting up Robin model assets...")
    print(f"Models will be downloaded to: {models_dir}")
    
    # Ensure models directory exists
    models_dir.mkdir(parents=True, exist_ok=True)
    
    github_token = os.getenv('GITHUB_TOKEN')
    if not github_token:
        print("Warning: No GITHUB_TOKEN environment variable set.")
        print("This may cause issues if the repository is private.")
        print("Set GITHUB_TOKEN to your personal access token if needed.")

    ok, msgs = update_models(models_dir=models_dir, github_token=github_token)
    for m in msgs:
        print(m)
    if not ok:
        sys.exit(1)

    print("\nAll model assets downloaded successfully!")
    print("You can now run Robin without Git LFS.")


if __name__ == "__main__":
    main()
