# methylartist_locus_capture.py
import sys
import shutil
import runpy
from contextlib import contextmanager, redirect_stderr
import os

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
        "methylartist",
        "locus",
        "-i",
        str(interval),
        "-b",
        str(bam_path),
        "-o",
        "ignored.png",  # required by CLI; we'll intercept savefig anyway
        "--motif",
        str(motif),
        "--mods",
        str(mods),
    ]
    if extra_cli:
        argv.extend([str(a) for a in extra_cli])

    # Final safeguard: ensure all argv entries are plain strings
    argv = [str(a) for a in argv]

    captured = {"fig": None}

    with _patch_fig_savefig(captured), _patch_argv(
        argv
    ), _suppress_methylartist_warnings():
        # Execute the installed CLI script as __main__
        runpy.run_path(cli, run_name="__main__")

    if captured["fig"] is None:
        raise RuntimeError("methylartist locus did not produce a figure.")

    fig = captured["fig"]

    # Always add MGMT CpG site annotations if we're in the MGMT interval
    # Check if interval matches MGMT region (chr10:129466536-129467536)
    try:
        import logging
        
        # Determine which sites to annotate
        mgmt_sites = []
        if "chr10" in interval and "129466" in interval:
            # Use site_rows if available, otherwise use default positions
            if site_rows:
                for row in site_rows:
                    pos_str = str(row.get("pos", ""))
                    site_label = str(row.get("site", ""))
                    # Extract just the number (e.g., "1" from "Site 1" or "1" from "Site 1 (CpG ...)")
                    site_num = site_label.split("(")[0].strip().replace("Site ", "").strip()
                    if not site_num:
                        # Fallback: try to extract number from the label
                        import re
                        match = re.search(r'\d+', site_label)
                        site_num = match.group(0) if match else site_label.strip()
                    
                    parts = [p.strip() for p in pos_str.split("/") if p.strip().isdigit()]
                    if len(parts) >= 2:
                        try:
                            p1 = float(parts[0])
                            p2 = float(parts[1])
                            mgmt_sites.append((p1, p2, site_num))
                        except (ValueError, IndexError):
                            continue
            
            # If no sites from site_rows, use default MGMT positions
            if not mgmt_sites:
                mgmt_sites = [
                    (129467255, 129467256, "1"),
                    (129467258, 129467259, "2"),
                    (129467262, 129467263, "3"),
                    (129467272, 129467273, "4"),
                ]
            
            # Get all axes in the figure
            axes = [ax for ax in fig.get_axes() if hasattr(ax, 'get_xlim')]
            
            if not axes:
                axes = [fig.gca()]
            
            # Parse the interval to get the start position (methylartist uses relative coordinates)
            interval_start = None
            try:
                # Parse interval like "chr10:129466536-129467536"
                if ":" in interval:
                    chrom, coords = interval.split(":", 1)
                    if "-" in coords:
                        start_str, end_str = coords.split("-", 1)
                        interval_start = int(start_str)
                        interval_end = int(end_str)
            except Exception:
                pass
            
            # Convert absolute genomic positions to relative positions (relative to interval start)
            if interval_start:
                relative_mgmt_sites = []
                for p1, p2, label in mgmt_sites:
                    rel_p1 = p1 - interval_start
                    rel_p2 = p2 - interval_start
                    relative_mgmt_sites.append((rel_p1, rel_p2, label))
                mgmt_sites = relative_mgmt_sites
            
            # Find the main axis for text annotations (only calculate once)
            main_ax = _choose_main_axis(fig)
            
            # Add vertical lines to ALL axes (all subplots)
            label_positions = []
            for ax_idx, ax in enumerate(axes):
                try:
                    xlim = ax.get_xlim()
                    
                    # For each MGMT site, draw vertical lines at both CpG positions
                    for p1, p2, label in mgmt_sites:
                        # Only draw if positions are within the visible x-axis range
                        if xlim[0] <= p1 <= xlim[1] or xlim[0] <= p2 <= xlim[1]:
                            # Draw dashed vertical lines - more subtle, lower zorder so data points appear on top
                            ax.axvline(
                                p1, 
                                color="crimson", 
                                linestyle="--",
                                linewidth=1.0,
                                alpha=0.5,
                                zorder=5
                            )
                            ax.axvline(
                                p2, 
                                color="crimson", 
                                linestyle="--",
                                linewidth=1.0,
                                alpha=0.5,
                                zorder=5
                            )
                            
                            # Add text annotation on the main axis only
                            if ax == main_ax:
                                ylim = ax.get_ylim()
                                y_range = ylim[1] - ylim[0]
                                center_x = (p1 + p2) / 2.0
                                
                                # Calculate vertical position with spacing to prevent overlap
                                # Stagger labels vertically: each site gets a different y position
                                site_index = next((i for i, (_, _, l) in enumerate(mgmt_sites) if l == label), 0)
                                base_y_offset = 0.03
                                spacing = 0.03
                                annotation_y_offset = base_y_offset + (site_index * spacing)
                                annotation_y = ylim[1] - (annotation_y_offset * y_range)
                                
                                # Check for horizontal overlap with existing labels
                                overlap_threshold = 50
                                has_overlap = False
                                for existing_x, existing_y in label_positions:
                                    if abs(center_x - existing_x) < overlap_threshold:
                                        has_overlap = True
                                        # If horizontal overlap, move vertically further down
                                        if abs(annotation_y - existing_y) < spacing * y_range:
                                            annotation_y -= spacing * y_range
                                        break
                                
                                if not has_overlap:
                                    label_positions.append((center_x, annotation_y))
                                
                                ax.text(
                                    center_x,
                                    annotation_y,
                                    label,
                                    ha="center",
                                    va="bottom",
                                    fontsize=9,
                                    fontweight="bold",
                                    color="crimson",
                                    bbox=dict(
                                        boxstyle="round,pad=0.3",
                                        facecolor="white",
                                        edgecolor="crimson",
                                        alpha=0.9,
                                        linewidth=1.5,
                                    ),
                                    zorder=10,
                                )
                except Exception as e:
                    logging.debug(f"[MGMT] Failed to add annotations to axis: {e}")
                    continue
            
            # Use constrained_layout instead of tight_layout to avoid warnings
            try:
                fig.set_constrained_layout(True)
            except Exception:
                # Fallback to tight_layout if constrained_layout fails
                try:
                    import warnings
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        fig.tight_layout()
                except Exception:
                    pass
            
            logging.debug(f"[MGMT] Successfully added annotations for {len(mgmt_sites)} sites")
        else:
            logging.debug(f"[MGMT] Interval {interval} is not MGMT region, skipping annotations")
    except Exception as e:
        import logging
        logging.exception(f"[MGMT] Failed to add annotations: {e}")

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
