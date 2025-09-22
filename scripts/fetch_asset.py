#!/usr/bin/env python3
"""
Asset fetcher for Robin project models.
Downloads and verifies release assets using the assets.json manifest.
"""

import json
import hashlib
import os
import sys
import urllib.request
import urllib.error
from pathlib import Path
from typing import Dict, Any, Optional


def load_assets_manifest() -> Dict[str, Any]:
    """Load the assets manifest from assets.json"""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    assets_file = project_root / "assets.json"
    
    if not assets_file.exists():
        raise FileNotFoundError(f"assets.json not found at {assets_file}")
    
    with open(assets_file, 'r') as f:
        return json.load(f)


def calculate_sha256(file_path: Path) -> str:
    """Calculate SHA256 hash of a file"""
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            sha256_hash.update(chunk)
    return sha256_hash.hexdigest()


def download_file(url: str, target_path: Path, github_token: Optional[str] = None) -> None:
    """Download a file from URL to target path"""
    headers = {}
    if github_token:
        headers["Authorization"] = f"Bearer {github_token}"
    
    request = urllib.request.Request(url, headers=headers)
    
    try:
        with urllib.request.urlopen(request) as response:
            with open(target_path, 'wb') as f:
                f.write(response.read())
    except urllib.error.HTTPError as e:
        if e.code == 401:
            raise RuntimeError("Authentication failed. Check your GitHub token.")
        elif e.code == 404:
            raise RuntimeError(f"Asset not found at {url}")
        else:
            raise RuntimeError(f"HTTP error {e.code}: {e.reason}")


def fetch_asset(asset_name: str, target_path: str, github_token: Optional[str] = None) -> None:
    """Fetch and verify a release asset"""
    manifest = load_assets_manifest()
    
    if asset_name not in manifest["assets"]:
        available = ", ".join(manifest["assets"].keys())
        raise ValueError(f"Asset '{asset_name}' not found. Available: {available}")
    
    asset_info = manifest["assets"][asset_name]
    asset_url = asset_info["url"]
    expected_sha256 = asset_info["sha256"]
    asset_filename = asset_info["name"]
    
    target_path = Path(target_path)
    
    print(f"Fetching {asset_filename}...")
    
    # Create target directory if it doesn't exist
    target_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Download the asset
    download_file(asset_url, target_path, github_token)
    
    # Verify SHA256 checksum
    print("Verifying checksum...")
    calculated_sha256 = calculate_sha256(target_path)
    
    if calculated_sha256 != expected_sha256:
        print(f"Error: SHA256 checksum mismatch!")
        print(f"Expected: {expected_sha256}")
        print(f"Got:      {calculated_sha256}")
        target_path.unlink()  # Remove the corrupted file
        sys.exit(1)
    
    print(f"Successfully downloaded and verified {asset_filename}")
    print(f"File: {target_path}")
    print(f"SHA256: {expected_sha256}")


def main():
    """Main entry point"""
    if len(sys.argv) < 3:
        print("Usage: python scripts/fetch_asset.py <asset_name> <target_path> [github_token]")
        print("Asset names: general_model, capper_model, pancan_model")
        sys.exit(1)
    
    asset_name = sys.argv[1]
    target_path = sys.argv[2]
    github_token = sys.argv[3] if len(sys.argv) > 3 else os.getenv('GITHUB_TOKEN')
    
    try:
        fetch_asset(asset_name, target_path, github_token)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
