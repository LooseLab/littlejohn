from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request
from email.utils import parsedate_to_datetime
from pathlib import Path
from typing import Optional
import pysam

logger = logging.getLogger("robin.clinvar")

CLINVAR_VCF_GZ_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"
)

CLINVAR_VCF_GZ_NAME = "clinvar.vcf.gz"
CLINVAR_VCF_GZ_TBI_NAME = "clinvar.vcf.gz.tbi"


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
            total: Optional[int] = None
            try:
                cl = resp.headers.get("Content-Length")
                if cl:
                    total = int(cl)
            except Exception:
                total = None

            # Import click lazily to keep module usable in non-CLI contexts.
            try:
                import click  # type: ignore
            except Exception:
                click = None  # type: ignore

            chunk_size = 1024 * 1024  # 1 MiB
            downloaded = 0

            if click is not None and total and total > 0:
                with click.progressbar(
                    length=total,
                    label="Downloading ClinVar",
                    show_eta=True,
                    show_percent=True,
                ) as bar:
                    with open(tmp_path, "wb") as f:
                        while True:
                            chunk = resp.read(chunk_size)
                            if not chunk:
                                break
                            f.write(chunk)
                            downloaded += len(chunk)
                            bar.update(len(chunk))
            else:
                # Fallback: stream and emit sparse size updates.
                last_reported_mb = -1
                with open(tmp_path, "wb") as f:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        f.write(chunk)
                        downloaded += len(chunk)
                        if click is not None:
                            mb = downloaded // (1024 * 1024)
                            if mb // 16 != last_reported_mb // 16:
                                last_reported_mb = mb
                                click.echo(f"Downloading ClinVar: {mb} MiB downloaded...")

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


def _ensure_tabix_index(gz_path: Path, tbi_path: Path) -> None:
    """
    Ensure a tabix index exists for a bgzipped VCF.
    """
    if _file_ok(tbi_path):
        return

    print(f"Creating tabix index for {gz_path.name}")
    logger.info("Creating tabix index for %s", gz_path)
    last_err: Optional[BaseException] = None
    try:
        pysam.tabix_index(
            str(gz_path),
            preset="vcf",
            force=True,
            keep_original=True,
        )
    except Exception as exc:
        last_err = exc
        logger.warning("pysam.tabix_index failed for %s: %s", gz_path, exc)
        tabix_bin = shutil.which("tabix")
        if tabix_bin:
            try:
                subprocess.run(
                    [tabix_bin, "-p", "vcf", str(gz_path)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
            except subprocess.CalledProcessError as sub_exc:
                stderr = (sub_exc.stderr or "").strip()
                stdout = (sub_exc.stdout or "").strip()
                detail = stderr or stdout or str(sub_exc)
                raise RuntimeError(
                    f"Tabix indexing failed for {gz_path} (pysam: {exc!r}; "
                    f"tabix CLI: {detail})"
                ) from sub_exc
        else:
            raise RuntimeError(
                f"Tabix indexing failed for {gz_path}: {exc!r}. "
                "Install htslib (tabix/bgzip) or fix the VCF: it must be "
                "block-gzipped (bgzip), not plain gzip, and a valid sorted VCF."
            ) from exc

    if not _file_ok(tbi_path):
        hint = (
            f"Tabix index was not created at {tbi_path}. "
            "If the VCF was recompressed with gzip instead of bgzip, run: "
            f"bgzip -d {gz_path.name} && bgzip {gz_path.with_suffix('').name} "
            f"&& tabix -p vcf {gz_path.name}"
        )
        if last_err:
            raise RuntimeError(hint) from last_err
        raise RuntimeError(hint)


def ensure_clinvar_files(
    *,
    resources_dir: Optional[Path] = None,
    download_if_missing: bool = True,
    url: str = CLINVAR_VCF_GZ_URL,
) -> None:
    """
    Ensure required ClinVar files exist in `robin.resources`:
    - `clinvar.vcf.gz`
    - `clinvar.vcf.gz.tbi`

    Strategy:
    - Download `clinvar.vcf.gz` when missing.
    - Create tabix index (`.tbi`) when missing.
    """

    resources_dir = resources_dir or get_resources_dir()
    gz_path = resources_dir / CLINVAR_VCF_GZ_NAME
    tbi_path = resources_dir / CLINVAR_VCF_GZ_TBI_NAME

    # If everything is already present, do nothing.
    if _file_ok(gz_path) and _file_ok(tbi_path):
        return

    gz_ok = _file_ok(gz_path)
    tbi_ok = _file_ok(tbi_path)

    print(
        "ClinVar status check: "
        f"{gz_path.name}={'OK' if gz_ok else 'MISSING'}; "
        f"{tbi_path.name}={'OK' if tbi_ok else 'MISSING'}"
    )

    if not gz_ok:
        if not download_if_missing:
            raise FileNotFoundError(
                f"Missing ClinVar file in {resources_dir}: expected {gz_path.name}"
            )
        _download_url_to_file(url, gz_path)
    _ensure_tabix_index(gz_path, tbi_path)


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
    tbi_path = resources_dir / CLINVAR_VCF_GZ_TBI_NAME

    # Ensure local ClinVar artifacts exist.
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
        _ensure_tabix_index(gz_path, tbi_path)
    finally:
        try:
            if tmp_path.exists():
                tmp_path.unlink()
        except OSError:
            pass

    logger.info("ClinVar updated successfully.")
    return True

