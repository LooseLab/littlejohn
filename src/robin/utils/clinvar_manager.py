from __future__ import annotations

import gzip
import logging
import os
import shutil
import tempfile
import urllib.error
import urllib.request
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Optional

logger = logging.getLogger("robin.clinvar")

CLINVAR_VCF_GZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)

CLINVAR_VCF_GZ_NAME = "clinvar.vcf.gz"
CLINVAR_VCF_NAME = "clinvar.vcf"


def get_resources_dir() -> Path:
    """
    Return the on-disk directory backing the `robin.resources` package.

    Notes:
    - This is expected to be writable in the typical ROBIN dev / workstation setup.
    - If it is not writable, the downloader will raise.
    """

    # Avoid importing `robin` at runtime.
    # `robin/__init__.py` imports analysis modules and can fail if optional
    # dependencies are not installed in minimal environments.
    this_file = Path(__file__).resolve()
    resources_dir = this_file.parent.parent / "resources"
    if not resources_dir.exists():
        raise RuntimeError(f"ROBIN resources directory not found: {resources_dir}")
    return resources_dir


def _file_ok(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _download_url_to_file(url: str, target_path: Path, *, timeout_s: int = 600) -> None:
    """
    Download `url` to `target_path` atomically (temp file + rename).
    """

    target_path.parent.mkdir(parents=True, exist_ok=True)

    # Use a temp file in the same directory so rename is atomic.
    with tempfile.NamedTemporaryFile(
        mode="wb", suffix=".part", prefix=target_path.name + ".", dir=str(target_path.parent), delete=False
    ) as tmp:
        tmp_path = Path(tmp.name)
    try:
        msg = f"Downloading ClinVar from {url}"
        print(msg)
        logger.info("Downloading ClinVar from %s", url)
        req = urllib.request.Request(url, headers={"User-Agent": "robin-clinvar/1.0"})
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            # Stream to disk.
            with open(tmp_path, "wb") as f:
                shutil.copyfileobj(resp, f)

        tmp_path.replace(target_path)
        print(f"ClinVar download complete: {target_path}")
        logger.info("Downloaded to %s", target_path)
    except Exception:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass
        raise


def _decompress_gz_to_vcf(gz_path: Path, vcf_path: Path) -> None:
    """
    Decompress `*.vcf.gz` to `*.vcf`.
    """

    vcf_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Decompressing {gz_path.name} -> {vcf_path.name}")
    logger.info("Decompressing %s -> %s", gz_path, vcf_path)
    with gzip.open(str(gz_path), "rb") as f_in:
        with open(str(vcf_path), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def _compress_vcf_to_gz(vcf_path: Path, gz_path: Path) -> None:
    """
    Compress `*.vcf` to `*.vcf.gz`.
    """

    gz_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Compressing {vcf_path.name} -> {gz_path.name}")
    logger.info("Compressing %s -> %s", vcf_path, gz_path)
    with open(str(vcf_path), "rb") as f_in:
        with gzip.open(str(gz_path), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def ensure_clinvar_files(
    *,
    resources_dir: Optional[Path] = None,
    download_if_missing: bool = True,
    url: str = CLINVAR_VCF_GZ_URL,
) -> None:
    """
    Ensure ClinVar files exist in `robin.resources`:
    - `clinvar.vcf.gz` (required by lightweight gene analysis)
    - `clinvar.vcf` (required by snpEff/snpSift steps)

    Strategy:
    - If one format is present, create the other locally (compress/decompress).
    - If neither exists, download `clinvar.vcf.gz`, then generate `clinvar.vcf`.
    """

    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    vcf_path = resources_dir / CLINVAR_VCF_NAME

    # If everything is already present, do nothing.
    if _file_ok(gz_path) and _file_ok(vcf_path):
        return

    gz_ok = _file_ok(gz_path)
    vcf_ok = _file_ok(vcf_path)

    print(
        "ClinVar status check: "
        f"{gz_path.name}={'OK' if gz_ok else 'MISSING'}; "
        f"{vcf_path.name}={'OK' if vcf_ok else 'MISSING'}"
    )

    if not gz_ok and not vcf_ok:
        if not download_if_missing:
            raise FileNotFoundError(
                f"Missing ClinVar files in {resources_dir}: "
                f"expected {gz_path.name} and {vcf_path.name}"
            )
        _download_url_to_file(url, gz_path)
        _decompress_gz_to_vcf(gz_path, vcf_path)
        return

    if not gz_ok and vcf_ok:
        # Avoid network if we already have the (uncompressed) VCF.
        _compress_vcf_to_gz(vcf_path, gz_path)

    if not vcf_ok and gz_ok:
        _decompress_gz_to_vcf(gz_path, vcf_path)


def _get_remote_last_modified(url: str, *, timeout_s: int = 30) -> Optional[float]:
    """
    Best-effort: fetch HTTP Last-Modified as a unix timestamp (seconds).
    Returns None when header is unavailable or request fails.
    """

    req = urllib.request.Request(url, headers={"User-Agent": "robin-clinvar/1.0"})
    req.method = "HEAD"
    try:
        with urllib.request.urlopen(req, timeout=timeout_s) as resp:
            last_modified = resp.headers.get("Last-Modified")
    except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError):
        return None

    if not last_modified:
        return None

    try:
        dt = parsedate_to_datetime(last_modified)
        return dt.timestamp()
    except Exception:
        return None


def update_clinvar_if_newer(
    *,
    resources_dir: Optional[Path] = None,
    url: str = CLINVAR_VCF_GZ_URL,
    download_if_missing: bool = True,
    timeout_s: int = 600,
) -> bool:
    """
    Update ClinVar if NCBI reports a newer `Last-Modified` than the local file.

    Returns:
        True if an update was performed, else False.
    """

    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    vcf_path = resources_dir / CLINVAR_VCF_NAME

    # Ensure we have at least one format locally so we can generate the other
    # after an update.
    ensure_clinvar_files(
        resources_dir=resources_dir,
        download_if_missing=download_if_missing,
        url=url,
    )

    if not _file_ok(gz_path):
        # Should not happen due to ensure_clinvar_files, but keep defensive.
        if not download_if_missing:
            return False
        _download_url_to_file(url, gz_path, timeout_s=timeout_s)

    local_mtime = gz_path.stat().st_mtime
    remote_mtime = _get_remote_last_modified(url)

    if remote_mtime is None:
        logger.info(
            "ClinVar update skipped: remote Last-Modified unavailable for %s", url
        )
        print("ClinVar update skipped: remote Last-Modified unavailable.")
        return False

    # Add a small tolerance to avoid re-downloading due to timestamp rounding.
    if remote_mtime <= local_mtime + 1:
        logger.info("ClinVar already up to date (local=%s, remote=%s).", local_mtime, remote_mtime)
        print("ClinVar already up to date.")
        return False

    # Download newer gz to temp then replace.
    tmp_path = Path(str(gz_path) + ".new")
    try:
        _download_url_to_file(url, tmp_path, timeout_s=timeout_s)
        tmp_path.replace(gz_path)
        # Keep the uncompressed copy in sync.
        _decompress_gz_to_vcf(gz_path, vcf_path)
    finally:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass

    logger.info("ClinVar updated successfully.")
    return True

