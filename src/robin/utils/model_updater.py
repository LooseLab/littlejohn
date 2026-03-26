from __future__ import annotations

import hashlib
import json
import os
import urllib.error
import urllib.request
from pathlib import Path
from typing import Optional, Tuple, List, Dict, Any


def _default_repo_root_guess() -> Optional[Path]:
    """
    Best-effort guess for repo root when running from source checkout.
    When installed as a package, this may not exist.
    """
    try:
        here = Path(__file__).resolve()
        # .../src/robin/utils/model_updater.py -> .../ (repo root)
        return here.parents[3]
    except Exception:
        return None


def _resolve_assets_manifest_path(manifest_path: Optional[str] = None) -> Optional[Path]:
    """
    Resolve assets manifest path.

    Priority:
    1) explicit argument
    2) env ROBIN_ASSETS_MANIFEST
    3) repo-root assets.json (dev checkout)
    4) packaged robin.resources/assets.json (installed package)
    """
    if manifest_path:
        return Path(manifest_path).expanduser().resolve()

    env = os.getenv("ROBIN_ASSETS_MANIFEST", "").strip()
    if env:
        return Path(env).expanduser().resolve()

    root = _default_repo_root_guess()
    if root:
        candidate = root / "assets.json"
        if candidate.exists():
            return candidate.resolve()

    # Installed-package fallback: look for robin.resources/assets.json
    try:
        import importlib.resources as importlib_resources

        p = importlib_resources.files("robin.resources").joinpath("assets.json")
        # `p` may be a Traversable; try converting to filesystem path.
        try:
            fs = Path(p)  # type: ignore[arg-type]
            if fs.exists():
                return fs.resolve()
        except TypeError:
            # Not a real filesystem path (e.g. zipped); handle below by reading content.
            return None
    except Exception:
        pass

    return None


def _load_assets_manifest(manifest_path: Optional[str]) -> Dict[str, Any]:
    """
    Load assets manifest JSON.
    If manifest_path is None, attempt to load packaged robin.resources/assets.json.
    """
    if manifest_path:
        mp = Path(manifest_path).expanduser().resolve()
        with mp.open("r", encoding="utf-8") as f:
            return json.load(f)

    # Packaged fallback
    import importlib.resources as importlib_resources

    txt = importlib_resources.files("robin.resources").joinpath("assets.json").read_text(encoding="utf-8")
    return json.loads(txt)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _download(url: str, target_path: Path, github_token: Optional[str]) -> None:
    headers: Dict[str, str] = {}
    if github_token:
        headers["Authorization"] = f"Bearer {github_token}"
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req) as resp:
        target_path.parent.mkdir(parents=True, exist_ok=True)
        with target_path.open("wb") as f:
            f.write(resp.read())


def update_models(
    models_dir: Path,
    github_token: Optional[str] = None,
    manifest_path: Optional[str] = None,
    overwrite: bool = False,
) -> Tuple[bool, List[str]]:
    """
    Download/verify required model assets into models_dir.

    Returns:
        (ok, messages)
    """
    messages: List[str] = []
    token = github_token if github_token is not None else os.getenv("GITHUB_TOKEN")

    # Lazy import so importing `robin` doesn't require heavy deps.
    try:
        from robin.utils.model_checker import get_required_models
    except Exception as e:
        return False, [f"Could not load required models list: {e}"]

    # Resolve manifest path (or fall back to packaged resources)
    mp = _resolve_assets_manifest_path(manifest_path)
    try:
        manifest = _load_assets_manifest(str(mp) if mp else None)
    except Exception:
        return False, [
            "Could not locate assets manifest (assets.json).",
            "Pass --manifest /path/to/assets.json or set ROBIN_ASSETS_MANIFEST=/path/to/assets.json.",
        ]

    models_dir = Path(models_dir).expanduser().resolve()
    models_dir.mkdir(parents=True, exist_ok=True)

    required = get_required_models()
    for asset_key, filename in required:
        if asset_key not in manifest.get("assets", {}):
            messages.append(f"Asset key not found in manifest: {asset_key}")
            return False, messages

        asset_info = manifest["assets"][asset_key]
        url = str(asset_info.get("url") or "")
        expected_sha256 = str(asset_info.get("sha256") or "")
        if not url or not expected_sha256:
            messages.append(f"Manifest entry incomplete for {asset_key} (missing url/sha256).")
            return False, messages

        target_path = models_dir / filename
        if target_path.exists() and not overwrite:
            messages.append(f"Skipping {filename} (already exists).")
            continue
        try:
            _download(url, target_path, token)
            got = _sha256_file(target_path)
            if got != expected_sha256:
                try:
                    target_path.unlink()
                except Exception:
                    pass
                messages.append(f"Checksum mismatch for {filename} (expected {expected_sha256}, got {got}).")
                return False, messages
            messages.append(f"Downloaded {filename}.")
        except Exception as e:
            messages.append(f"Failed to download {filename}: {e}")
            return False, messages

    return True, messages

