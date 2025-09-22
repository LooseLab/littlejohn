#!/bin/bash

# fetch_asset.sh - Download and verify release assets
# Usage: ./scripts/fetch_asset.sh <asset_name> <target_path> [github_token]

set -e

ASSET_NAME="$1"
TARGET_PATH="$2"
GITHUB_TOKEN="${3:-$GITHUB_TOKEN}"

if [ -z "$ASSET_NAME" ] || [ -z "$TARGET_PATH" ]; then
    echo "Usage: $0 <asset_name> <target_path> [github_token]"
    echo "Asset names: general_model, capper_model, pancan_model"
    exit 1
fi

# Load assets manifest
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
ASSETS_MANIFEST="$PROJECT_ROOT/assets.json"

if [ ! -f "$ASSETS_MANIFEST" ]; then
    echo "Error: assets.json not found at $ASSETS_MANIFEST"
    exit 1
fi

# Extract asset info from manifest
ASSET_URL=$(jq -r ".assets.$ASSET_NAME.url" "$ASSETS_MANIFEST")
ASSET_SHA256=$(jq -r ".assets.$ASSET_NAME.sha256" "$ASSETS_MANIFEST")
ASSET_NAME_FILE=$(jq -r ".assets.$ASSET_NAME.name" "$ASSETS_MANIFEST")

if [ "$ASSET_URL" = "null" ] || [ "$ASSET_SHA256" = "null" ]; then
    echo "Error: Asset '$ASSET_NAME' not found in manifest"
    exit 1
fi

echo "Fetching $ASSET_NAME_FILE..."

# Create target directory if it doesn't exist
mkdir -p "$(dirname "$TARGET_PATH")"

# Download the asset
if [ -n "$GITHUB_TOKEN" ]; then
    curl -L -H "Authorization: Bearer $GITHUB_TOKEN" "$ASSET_URL" -o "$TARGET_PATH"
else
    curl -L "$ASSET_URL" -o "$TARGET_PATH"
fi

# Verify SHA256 checksum
echo "Verifying checksum..."
CALCULATED_SHA256=$(shasum -a 256 "$TARGET_PATH" | cut -d' ' -f1)

if [ "$CALCULATED_SHA256" != "$ASSET_SHA256" ]; then
    echo "Error: SHA256 checksum mismatch!"
    echo "Expected: $ASSET_SHA256"
    echo "Got:      $CALCULATED_SHA256"
    rm -f "$TARGET_PATH"
    exit 1
fi

echo "Successfully downloaded and verified $ASSET_NAME_FILE"
echo "File: $TARGET_PATH"
echo "SHA256: $ASSET_SHA256"
