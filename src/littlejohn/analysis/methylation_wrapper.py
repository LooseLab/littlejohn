# methylartist_locus_capture.py
import sys
import shutil
import runpy
from contextlib import contextmanager, redirect_stderr
import os
import io
from matplotlib.figure import Figure
from typing import List, Optional, Dict, Any

@contextmanager
def _patch_fig_savefig(capture_dict):
    orig = Figure.savefig
    def _capture(self, *args, **kwargs):
        capture_dict["fig"] = self
        # swallow file I/O
        return None
    Figure.savefig = _capture
    try:
        yield
    finally:
        Figure.savefig = orig

@contextmanager
def _patch_argv(argv):
    old = sys.argv
    sys.argv = argv
    try:
        yield
    finally:
        sys.argv = old

@contextmanager
def _suppress_methylartist_warnings():
    """Temporarily suppress WARNING-level logs and stderr noise from methylartist."""
    import logging
    previous_disable = logging.root.manager.disable
    # Disable WARNING and below (keeps ERROR/CRITICAL visible)
    logging.disable(logging.WARNING)
    devnull = open(os.devnull, "w")
    try:
        with redirect_stderr(devnull):
            yield
    finally:
        logging.disable(previous_disable)
        try:
            devnull.close()
        except Exception:
            pass

def _find_cli_or_raise() -> str:
    cli = shutil.which("methylartist")
    if not cli:
        raise RuntimeError(
            "Couldn't find 'methylartist' on PATH. "
            "Activate the environment where you installed it, "
            "or add the env's bin/Scripts directory to PATH."
        )
    return cli

def _choose_main_axis(fig: Figure):
    axes = [ax for ax in getattr(fig, "axes", []) if hasattr(ax, "get_xlim")]
    if not axes:
        return fig.gca()
    spans = []
    for ax in axes:
        try:
            x0, x1 = ax.get_xlim()
            spans.append((ax, abs(x1 - x0)))
        except Exception:
            continue
    if not spans:
        return fig.gca()
    spans.sort(key=lambda t: t[1], reverse=True)
    return spans[0][0]

def locus_figure(
    interval: str,
    bam_path: str,
    motif: str = "CG",
    mods: str = "m",
    extra_cli: Optional[List[str]] = None,
    site_rows: Optional[List[Dict[str, Any]]] = None,
):
    """
    Run the installed 'methylartist locus' command but return the Matplotlib Figure
    object instead of writing an image to disk.

    Parameters
    ----------
    interval : str          e.g. 'chr10:129466536-129467536'
    bam_path : str          path to BAM with MM/ML tags
    motif : str             e.g. 'CG' (default)
    mods : str              e.g. 'm' for 5mC (default)
    extra_cli : list[str]   any extra CLI flags, e.g. ["--minqual","20","--reads","2000"]

    Returns
    -------
    matplotlib.figure.Figure
    """
    cli = _find_cli_or_raise()

    argv = [
        "methylartist", "locus",
        "-i", str(interval),
        "-b", str(bam_path),
        "-o", "ignored.png",   # required by CLI; we'll intercept savefig anyway
        "--motif", str(motif),
        "--mods", str(mods),
    ]
    if extra_cli:
        argv.extend([str(a) for a in extra_cli])

    # Final safeguard: ensure all argv entries are plain strings
    argv = [str(a) for a in argv]

    captured = {"fig": None}

    with _patch_fig_savefig(captured), _patch_argv(argv), _suppress_methylartist_warnings():
        # Execute the installed CLI script as __main__
        runpy.run_path(cli, run_name="__main__")

    if captured["fig"] is None:
        raise RuntimeError("methylartist locus did not produce a figure.")

    fig = captured["fig"]

    # Optionally overlay vertical lines at CpG site positions
    if site_rows:
        try:
            ax = _choose_main_axis(fig)
            for row in site_rows:
                pos_str = str(row.get("pos", ""))
                parts = [p.strip() for p in pos_str.split("/") if p.strip().isdigit()]
                for p in parts:
                    try:
                        x = float(p)
                        ax.axvline(x, color="crimson", linestyle="--", linewidth=0.8, alpha=0.5)
                    except Exception:
                        continue
            fig.tight_layout()
        except Exception:
            pass

    return fig


# --- Optional: persistence helpers (same-Python/Matplotlib only) ---
def save_figure_pickle(fig: Figure, path: str) -> None:
    """Save a Matplotlib Figure object for later re-use (not portable across versions)."""
    import pickle
    with open(path, "wb") as f:
        pickle.dump(fig, f)

def load_figure_pickle(path: str) -> Figure:
    """Reload a pickled Matplotlib Figure object."""
    import pickle
    with open(path, "rb") as f:
        return pickle.load(f)

# --- Optional: portable-ish JSON via mpld3 (fidelity may vary) ---
def save_figure_mpld3(fig: Figure, path: str) -> None:
    import mpld3
    with open(path, "w") as f:
        f.write(mpld3.fig_to_json(fig))

def load_figure_mpld3(path: str) -> Figure:
    import mpld3
    with open(path, "r") as f:
        return mpld3.json_to_fig(f.read())


