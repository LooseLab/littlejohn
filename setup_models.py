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

from scripts.fetch_asset import fetch_asset


def main():
    """Download all required model assets"""
    models_dir = project_root / "src" / "robin" / "models"
    
    print("Setting up Robin model assets...")
    print(f"Models will be downloaded to: {models_dir}")
    
    # Ensure models directory exists
    models_dir.mkdir(parents=True, exist_ok=True)
    
    assets = [
        ("general_model", "general.zip"),
        ("capper_model", "Capper_et_al_NN.pkl"),
        ("pancan_model", "pancan_devel_v5i_NN.pkl")
    ]
    
    github_token = os.getenv('GITHUB_TOKEN')
    if not github_token:
        print("Warning: No GITHUB_TOKEN environment variable set.")
        print("This may cause issues if the repository is private.")
        print("Set GITHUB_TOKEN to your personal access token if needed.")
    
    for asset_name, filename in assets:
        target_path = models_dir / filename
        
        if target_path.exists():
            print(f"Skipping {filename} - already exists")
            continue
        
        try:
            print(f"\nDownloading {asset_name}...")
            fetch_asset(asset_name, str(target_path), github_token)
        except Exception as e:
            print(f"Failed to download {asset_name}: {e}")
            sys.exit(1)
    
    print("\nAll model assets downloaded successfully!")
    print("You can now run Robin without Git LFS.")


if __name__ == "__main__":
    main()
