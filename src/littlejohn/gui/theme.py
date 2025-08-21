"""
Module: theme

This module defines the theme and layout for the entire application, ensuring a consistent look and feel across all pages.

It includes:

- A context manager `frame` to create a custom page frame with navigation, header, and footer.
- Utility functions to handle dark mode and remote access toggling.
- Paths to external resources like images, HTML, and CSS files for styling.

Functions:

- frame(navtitle: str): Context manager for creating a consistent page layout with header, footer, and navigation.
- cleanup_and_exit(): Handles cleanup operations before shutting down the application.
- dark_mode(event: events.ValueChangeEventArguments): Toggles dark mode based on the event argument.
- use_on_air(args: events.ValueChangeEventArguments): Toggles remote access based on the event argument.

Constants:

- IMAGEFILE: Path to the image file used in the header and footer.
- HEADER_HTML: HTML content for the header.
- STYLE_CSS: CSS styles for the application.

External Dependencies:

- contextlib.contextmanager
- nicegui (ui, app, events, core, air)
- pathlib.Path
- robin.images
- os
- psutil
- platform
"""

from contextlib import contextmanager
from packaging import version
import requests
import asyncio
import logging
import subprocess
import importlib.metadata


from nicegui import ui, app, events, core, run
import nicegui.air

from pathlib import Path

from littlejohn.gui import images
from robin import __about__
from robin.core.state import state

import os
import psutil
import platform


def get_version_from_github():
    response = requests.get(
        "https://raw.githubusercontent.com/LooseLab/ROBIN/main/src/robin/__about__.py"
    )
    response.raise_for_status()
    remote_version_str = None
    for line in response.text.split("\n"):
        if line.startswith("__version__"):
            remote_version_str = line.split("=")[1].strip().strip('"').strip("'")
            break
    return remote_version_str


async def check_version():
    """
    Check the current version against the remote version on GitHub.
    Shows a notification or dialog to the user about their version status.
    """
    # Check if version has already been checked in this session
    if app.storage.tab.get("version_checked", False):
        return

    try:
        remote_version_str = await run.io_bound(get_version_from_github)

        if not remote_version_str:
            with ui.dialog() as dialog, ui.card():
                ui.label("Version Check Error").classes("text-h6")
                ui.label("Could not determine remote version. Please check manually.")
                ui.button("OK", on_click=dialog.close)
            dialog.open()
            return

        local_version = version.parse(__about__.__version__)
        remote_version = version.parse(remote_version_str)

        if local_version == remote_version:
            ui.notify("Your ROBIN installation is up to date!", type="positive")
        elif local_version < remote_version:
            with ui.dialog() as dialog, ui.card():
                ui.label("Update Available!").classes("text-h6")
                ui.label(f"Your version: {local_version}")
                ui.label(f"Latest version: {remote_version}")
                ui.label("Would you like to visit the GitHub repository to update?")
                with ui.row():
                    ui.button("Continue with current version", on_click=dialog.close)
                    ui.button(
                        "Visit GitHub",
                        on_click=lambda: ui.open("https://github.com/LooseLab/ROBIN"),
                    ).classes("bg-primary")
            dialog.open()
        else:
            with ui.dialog() as dialog, ui.card():
                ui.label("Development Version").classes("text-h6")
                ui.label(f"You are running a development version ({local_version}).")
                ui.label(f"Latest release: {remote_version}")
                ui.label(
                    "This version may be unstable and is only for testing purposes. It is not recommended for production use."
                )
                ui.label("Please consider using the latest release instead.")
                ui.button("OK", on_click=dialog.close)
            dialog.open()

    except requests.RequestException:
        with ui.dialog() as dialog, ui.card():
            ui.label("Connection Error").classes("text-h6")
            ui.label("Could not check for updates.")
            ui.label(
                "Either you are not connected to the internet or you cannot access https://www.github.com/looselab/robin."
            )
            ui.label("Please manually check for updates.")
            ui.button("OK", on_click=dialog.close)
        dialog.open()
    except Exception as e:
        with ui.dialog() as dialog, ui.card():
            ui.label("Error").classes("text-h6")
            ui.label(f"Error checking version: {str(e)}")
            ui.button("OK", on_click=dialog.close)
        dialog.open()

    # Mark version as checked for this session
    app.storage.tab["version_checked"] = True


# Define the path to the image file used in the header and footer
IMAGEFILE = os.path.join(
    os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
)

# Module-level variables
quitdialog = None

MENU_BREAKPOINT = 1200

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()


@contextmanager
def frame(navtitle: str, batphone=False, smalltitle=None):
    """
    Context manager to create a custom page frame with consistent styling and behavior across all pages.

    Args:
        navtitle (str): The title to display in the navigation header.
        batphone (bool): Whether to show the BATMAN mode title.
        smalltitle (str): The title to display on small screens.

    Yields:
        None
    """
    global quitdialog
    if batphone:
        navtitle = f"BATMAN & {navtitle}"

    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.2.0/dist/igv.min.js"></script>'
    )
    ui.add_head_html(HEADER_HTML + f"<style>{STYLE_CSS}</style>")
    ui.add_head_html(
        """
        <script>
        function emitSize() {
            emitEvent('resize', {
                width: document.body.offsetWidth,
                height: document.body.offsetHeight,
            });
        }
        window.onload = emitSize;
        window.onresize = emitSize;
        </script>
    """
    )

    # Create disclaimer dialog that appears on first visit
    async def show_disclaimer():
        # await ui.context.client.connected()
        if not app.storage.tab.get("disclaimer_acknowledged", False):
            with ui.dialog().props(
                "persistent"
            ) as disclaimer_dialog, ui.card().classes("w-160"):
                ui.label("DISCLAIMER").classes("text-h5 text-weight-bold q-mb-md")
                ui.label(
                    "This tool and the data generated by it are intended for research use only and should not be used for "
                    "direct diagnostic purposes. The methylation-based classifications and other analyses provided here may "
                    "be considered by neuropathologists as supplementary information in the context of comprehensive "
                    "diagnostic assessment, which should include clinical history, radiological findings, and complete "
                    "histopathological and molecular evaluation. The final interpretation and diagnosis should always be "
                    "made by qualified healthcare professionals based on all available information."
                ).classes("text-body1 q-mb-md")

                if batphone:
                    ui.label("BATMAN Mode").classes("text-h5 text-weight-bold q-mb-md")
                    ui.label(
                        "You are running this tool in BATMAN mode. "
                        "This is a beta version of the tool and may not be fully functional. "
                        "BATMAN means: Breakpoint Adaptive Targeting alongside Methylation Analysis on Nanopore. "
                        "This means that the target regions will be updated in real-time based on detected breakpoints. "
                        "This code only works with ReadFish at this time. "
                    )

                def acknowledge():
                    app.storage.tab["disclaimer_acknowledged"] = True
                    disclaimer_dialog.close()

                ui.button("I agree", on_click=acknowledge).props("color=primary")
            disclaimer_dialog.open()

    ui.timer(0.5, show_disclaimer, once=True)

    # Add version check timer

    ui.timer(1.0, check_version, once=True)

    # Create a persistent dialog for quitting the app
    quitdialog = ui.dialog().props("persistent")

    async def quit_app():
        quitdialog.close()
        await cleanup_and_exit()

    with quitdialog, ui.card():
        ui.label(
            "Quitting the app will stop running methylation analysis. Are you sure?"
        )
        ui.label("If you want to keep analysis running, click Cancel.")
        ui.label(
            "You can safely close this window and analysis will keep running in the background."
        )
        ui.button("Cancel", on_click=quitdialog.close).props("outline").classes(
            "shadow-lg"
        )
        ui.button("Really Quit", icon="logout", on_click=quit_app).props(
            "outline"
        ).classes("shadow-lg")

    # Create a header with navigation title and menu
    header_classes = "items-center duration-200 p-0 px-4 no-wrap"
    if batphone:
        header_classes += " batphone"

    with ui.header(elevated=True).classes(header_classes):
        with ui.grid(columns=2).style("width: 100%"):
            with ui.row().classes(
                f"max-[{MENU_BREAKPOINT}px]:hidden items-center align-left"
            ):
                ui.html(navtitle).classes("shadows-into").style(
                    "font-size: 150%; font-weight: 300"
                ).tailwind("drop-shadow", "font-bold")
            with ui.row().classes(
                f"min-[{MENU_BREAKPOINT+1}px]:hidden items-center align-left"
            ):
                ui.html(smalltitle).style("font-size: 150%; font-weight: 300").tailwind(
                    "drop-shadow", "font-bold"
                )
            with ui.row().classes("ml-auto align-top"):
                with ui.row().classes("items-center m-auto"):
                    ui.label(f"Viewing: {platform.node()}").classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.label("CPU").classes(f"max-[{MENU_BREAKPOINT}px]:hidden")
                    cpu_activity = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.label("RAM").classes(f"max-[{MENU_BREAKPOINT}px]:hidden")
                    ram_utilisation = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )

                    # Create a data model for system metrics
                    class SystemMetrics:
                        def __init__(self):
                            self.cpu = 0
                            self.ram = 0

                    metrics = SystemMetrics()

                    # Bind the progress indicators to the model
                    cpu_activity.bind_value_from(metrics, "cpu")
                    ram_utilisation.bind_value_from(metrics, "ram")

                    # Single timer to update both metrics
                    def update_metrics():
                        metrics.cpu = round(
                            psutil.getloadavg()[1] / os.cpu_count() * 100, 1
                        )
                        metrics.ram = round(psutil.virtual_memory()[2], 1)

                    ui.timer(1.0, update_metrics)

                    with ui.button(icon="menu"):
                        with ui.menu() as menu:
                            ui.menu_item("Home", lambda: ui.navigate.to("/"))
                            ui.menu_item("Live Data", lambda: ui.navigate.to("/live"))
                            ui.menu_item(
                                "Browse Historic Data",
                                lambda: ui.navigate.to("/browse"),
                            )
                            ui.menu_item(
                                "Workflow",
                                lambda: ui.navigate.to("/workflow"),
                            )
                            ui.separator()
                            ui.switch("Allow Remote Access").classes(
                                "ml-4 bg-transparent"
                            ).props('color="black"').bind_value(
                                app.storage.general, "use_on_air"
                            )
                            ui.separator()
                            ui.switch("Dark Mode").classes("ml-4 bg-transparent").props(
                                'color="black"'
                            ).bind_value(app.storage.browser, "dark_mode")
                            ui.dark_mode().bind_value(app.storage.browser, "dark_mode")
                            ui.separator()
                            ui.menu_item("Close", menu.close)
                            ui.button("Quit", icon="logout", on_click=quitdialog.open)
                    ui.image(IMAGEFILE).style("width: 50px")

    with ui.column().classes("w-full h-full") as main_content:
        pass

    # Activity monitor removed

    # Create a footer with useful information and quit button
    footer_classes = "items-center"
    if batphone:
        footer_classes += " batphone"
    with ui.footer().classes(footer_classes):
        with ui.dialog() as dialog, ui.card():
            ui.label("Links").tailwind("text-2xl font-bold font-italic drop-shadow")
            ui.separator()
            ui.link("Code on GitHub", "https://github.com/looselab/robin")
            ui.link(
                "Rapid CNS2 Paper",
                "https://link.springer.com/article/10.1007/s00401-022-02415-6",
            )
            ui.link(
                "Sturgeon Classifier",
                "https://www.nature.com/articles/s41586-023-06615-2",
            )
            ui.link(
                "Protocol",
                "https://www.protocols.io/view/intra-operative-nanopore-sequencing-to-classify-br-c65qzg5w",
            )
            ui.link("Oxford Nanopore", "https://nanoporetech.com/")
            ui.link("epi2me labs", "https://labs.epi2me.io/")
            ui.link("Looselab", "https://looselab.github.io/")
            ui.button("Close", on_click=dialog.close)
        ui.image(IMAGEFILE).style("width: 40px")
        ui.colors(primary="#555")
        ui.button("Links", on_click=dialog.open)

        with ui.button(icon="info"):
            with ui.menu() as menu:
                ui.label().bind_text_from(
                    app, "urls", backward=lambda n: f"Available urls: {n}"
                )
                ui.label("Version: " + __about__.__version__)
        ui.label(
            "Some aspects of this application are ©Looselab - all analyses provided for research use only."
        ).classes(f"max-[{MENU_BREAKPOINT}px]:hidden").tailwind("text-sm font-italic")
        ui.label("©Looselab").classes(f"min-[{MENU_BREAKPOINT+1}px]:hidden").tailwind(
            "text-sm font-italic"
        )
        ui.label("Not for diagnostic use.").classes(
            f"min-[{MENU_BREAKPOINT+1}px]:hidden"
        ).tailwind("text-sm font-italic")

    with main_content:
        yield


async def cleanup_and_exit():
    """
    Handle any necessary cleanup operations before exiting the application and then shut down the application.

    Returns:
        None

    Example:
        >>> cleanup_and_exit()
        None
    """
    logging.info("User initiated shutdown via UI")

    # Create and show shutdown modal
    with ui.dialog().props("persistent") as shutdown_dialog, ui.card().classes("w-96"):
        ui.label("Shutting Down").classes("text-h5 text-weight-bold q-mb-md")
        ui.label("ROBIN is shutting down. Please wait while we clean up...").classes(
            "text-body1 q-mb-md"
        )
        with ui.row().classes("w-full justify-center"):
            ui.spinner(size="lg", color="primary")

    # Close the quit dialog if it's open
    if quitdialog:
        quitdialog.close()

    # Show the shutdown dialog
    shutdown_dialog.open()
    logging.info("Shutdown dialog opened")

    # Set shutdown event
    state.shutdown_event = True
    logging.info("Shutdown event set")

    # Wait a moment to ensure the dialog is visible
    await asyncio.sleep(1.0)
    logging.info("Initial wait complete")

    # Perform cleanup
    logging.info("Performing cleanup operations...")
    logging.info("Shutting down ROBIN... from theme.py")
    logging.info(
        "Here we need to do some very graceful shutdown to make sure we don't leave any threads running and we don't leave any files open."
    )

    # Wait for cleanup to complete
    while state.shutdown_event and state.get_running_process_count() > 0:
        await asyncio.sleep(3.0)
    logging.info("Cleanup wait complete")

    # Close the shutdown dialog
    shutdown_dialog.close()
    logging.info("Shutdown dialog closed")

    # Shutdown the application
    logging.info("Application shutdown initiated")
    app.shutdown()


def use_on_air(args: events.ValueChangeEventArguments):
    """
    Enable or disable remote access based on the value of the event argument.

    Args:
        args (events.ValueChangeEventArguments): The event argument containing the value for remote access toggle.

    Returns:
        None

    Example:
        >>> args = events.ValueChangeEventArguments(value=True)
        >>> use_on_air(args)
        None
    """
    if args.value:
        if core.air is None:
            core.air = nicegui.air.Air("")
        nicegui.air.connect()
    else:
        nicegui.air.disconnect()


@ui.page("/")
def my_page():
    with frame(
        "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
        smalltitle="<strong>R.O.B.I.N</strong>",
    ):
        ui.label("Welcome to the Application")


@ui.page("/workflow")
def workflow_page():
    """Display the ROBIN workflow diagram."""
    with frame(
        "ROBIN Workflow",
        smalltitle="Workflow",
    ):
        # Get versions for the diagram
        # modkit_version = get_modkit_version()
        sturgeon_version = get_sturgeon_version()
        crossnn_version = get_crossnn_version()
        cnv_from_bam_version = get_cnv_from_bam_version()

        ui.mermaid(
            f"""
flowchart TD
    %% Style definitions
    classDef minKNOW fill:#f5fafd80,stroke:#339af0,stroke-width:1px,color:#1c7ed6,font-size:14px,font-weight:500
    classDef robin fill:#f6fcf7,stroke:#495057,stroke-width:1px,color:#388e3c,font-size:14px,font-weight:500
    classDef classifier fill:#fff0f699,stroke:#e64980,stroke-width:1px,color:#c2255c,font-size:14px,font-weight:500
    classDef analysis fill:#fff4e699,stroke:#fd7e14,stroke-width:1px,color:#e8590c,font-size:14px,font-weight:500
    classDef output fill:#ebfbee99,stroke:#40c057,stroke-width:1px,color:#2b8a3e,font-size:14px,font-weight:500
    classDef headerLabel fill:#ffffff00,stroke:#ffffff00,color:#222,font-size:18px,font-weight:700

    %% Custom header nodes
    robinLabel["<b>R.O.B.I.N v{__about__.__version__}</b>"]:::headerLabel
    minKNOWLabel["<b>MinKNOW Pipeline</b>"]:::headerLabel

    %% MinKNOW pipeline
    subgraph MinKNOW[" "]
        sequencing["Sequencing"]
        adaptive["Adaptive Sampling"]
        alignment["Alignment"]
        bam["BAM Files"]
    end
    minKNOWLabel -.-> MinKNOW

    %% R.O.B.I.N. pipeline
    subgraph ROBIN[" "]
        runInfo["Extract Run Information"]
        groupRuns["Group by Run"]
        mergeBam["Merge BAM Files"]
        methylation["Extract Methylation<br>(modkit retired)"]
        individualBam["Process BAMs Individually"]
        cnv["CNV Analysis<br>(CNV from BAM v{cnv_from_bam_version})"]
        fusion["Fusion Analysis"]
        mgmt["MGMT Analysis"]
        coverage["Coverage Analysis"]
        snv["Optional SNP/V Analysis<br>(ClairS-To)"]
        crossnnCNS["CrossNN - CNS<br>({crossnn_version})"]
        crossnnPan["CrossNN - PanCancer<br>({crossnn_version})"]
        sturgeon["Sturgeon<br>({sturgeon_version})"]
        randomForest["Random Forest"]
        integrate["Integrate Results"]
        report["Report Generation"]
    end
    robinLabel -.-> ROBIN

    %% Connections
    sequencing --> adaptive
    adaptive --> alignment
    alignment --> bam
    bam --> runInfo
    runInfo --> groupRuns

    %% Split after groupRuns
    groupRuns -->|"For methylation"| mergeBam
    groupRuns -->|"For other analyses"| individualBam

    mergeBam --> methylation
    methylation --> crossnnCNS & crossnnPan & sturgeon & randomForest
    crossnnCNS --> integrate
    crossnnPan --> integrate
    sturgeon --> integrate
    randomForest --> integrate

    individualBam --> cnv & fusion & mgmt & coverage
    coverage --> snv

    integrate --> report
    cnv --> report
    fusion --> report
    mgmt --> report
    snv --> report

    %% Styling
    class sequencing,adaptive,alignment,bam minKNOW
    class runInfo,groupRuns,mergeBam,individualBam,integrate robin
    class crossnnCNS,crossnnPan,sturgeon,randomForest classifier
    class cnv,fusion,mgmt,coverage,snv analysis
    class report output
    %% Subgraph background coloring
    style MinKNOW fill:#f5fafd80,stroke:#339af0,stroke-width:2px
    style ROBIN fill:#f6fcf7,stroke:#388e3c,stroke-width:2px
""",
            config={
                "theme": "redux",
                "look": "neo",
                "flowchart": {"curve": "basis", "defaultRenderer": "elk"},
            },
        ).classes("w-full")


def get_modkit_version():
    """Get the version of modkit installed."""
    try:
        result = subprocess.run(["modkit", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "unknown"


def get_sturgeon_version():
    """Get the version of sturgeon installed."""
    try:
        result = subprocess.run(
            ["sturgeon", "--version"], capture_output=True, text=True
        )
        if result.returncode == 0:
            return result.stdout.strip()
        return "unknown"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        try:
            return importlib.metadata.version("sturgeon")
        except importlib.metadata.PackageNotFoundError:
            return "unknown"


def get_crossnn_version():
    """Get the version of CrossNN installed."""
    try:
        return importlib.metadata.version("nanoDX")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/nanoDX",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


def get_cnv_from_bam_version():
    """Get the version of CNV from BAM installed."""
    try:
        return importlib.metadata.version("cnv_from_bam")
    except importlib.metadata.PackageNotFoundError:
        try:
            result = subprocess.run(
                ["git", "describe", "--tags"],
                cwd=os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "submodules/cnv_from_bam",
                ),
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                return result.stdout.strip()
            return "unknown"
        except (subprocess.SubprocessError, FileNotFoundError, OSError):
            return "unknown"


def main():
    """
    Main function to test the theme by creating a simple page using the frame context manager.

    Example:
        >>> main()
        None
    """
    # Add some custom CSS because - why not!
    ui.add_css(
        """
        .shadows-into light-regular {
            font-family: "Shadows Into Light", cursive;
            font-weight: 800;
            font-style: normal;
        }
    """
    )
    # Register some fonts that we might need later on (guarded)
    # try:
    fonts_dir = Path(__file__).parent / "fonts"
    if fonts_dir.exists():
        app.add_static_files("/fonts", str(fonts_dir))
    # except Exception:
    #    pass
    ui.run(storage_secret="robin")


if __name__ in {"__main__", "__mp_main__"}:
    # import doctest
    # doctest.testmod()
    main()


def get_process_ram_usage():
    """Get RAM usage for Python and R processes."""
    python_ram = 0
    r_ram = 0

    for proc in psutil.process_iter(["pid", "name", "memory_info"]):
        try:
            # Get process name and memory info
            name = proc.info["name"].lower()
            mem_info = proc.info["memory_info"]

            # Skip if memory_info is None
            if mem_info is None:
                continue

            # Calculate RAM usage in GB
            ram_gb = mem_info.rss / (1024 * 1024 * 1024)  # Convert bytes to GB

            # Only count processes that are actually using memory
            if ram_gb > 0:
                if "python" in name:
                    python_ram += ram_gb
                elif "r" in name or "Rscript" in name:
                    r_ram += ram_gb

        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            continue
        except Exception as e:
            logging.debug(
                f"Error getting memory info for process {proc.info.get('name', 'unknown')}: {str(e)}"
            )
            continue

    return round(python_ram, 2), round(
        r_ram, 2
    )  # Round to 2 decimal places for cleaner display


# Activity monitor RAM collection removed


def debug_process_tree():
    """
    Debug function to understand why we're picking up incorrect child processes.
    Shows the full process tree and explains why each process is being counted.
    """
    try:
        main_process = psutil.Process(os.getpid())

        logging.debug("=== Process Tree Debug ===")
        logging.debug(
            f"Main ROBIN process: {main_process.name()} (PID: {main_process.pid})"
        )
        logging.debug(f"Command line: {' '.join(main_process.cmdline())}")
        logging.debug(f"Parent PID: {main_process.ppid()}")

        # Get all processes in the system
        all_processes = []
        for proc in psutil.process_iter(["pid", "name", "ppid", "cmdline"]):
            try:
                all_processes.append(proc.info)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue

        # Find all processes that could be considered "children"
        potential_children = []
        for proc_info in all_processes:
            pid = proc_info["pid"]
            ppid = proc_info["ppid"]

            # Check if this process is in our process tree
            if _is_in_process_tree(pid, main_process.pid, all_processes):
                potential_children.append(proc_info)

        logging.debug(
            f"\nFound {len(potential_children)} processes in ROBIN's process tree:"
        )

        for i, proc_info in enumerate(potential_children):
            pid = proc_info["pid"]
            name = proc_info["name"]
            ppid = proc_info["ppid"]
            cmdline = " ".join(proc_info["cmdline"]) if proc_info["cmdline"] else "N/A"

            # Determine why this process is being counted
            reason = _explain_process_inclusion(proc_info, main_process.pid)

            logging.debug(f"\n{i+1}. {name} (PID: {pid}, PPID: {ppid})")
            logging.debug(f"   Command: {cmdline}")
            logging.debug(f"   Reason: {reason}")

            # Show memory usage if available
            try:
                proc = psutil.Process(pid)
                rss = proc.memory_info().rss / (1024 * 1024 * 1024)
                logging.debug(f"   Memory: {rss:.2f}GB RSS")
            except:
                logging.debug("   Memory: <access denied>")

        logging.debug("\n=== End Process Tree Debug ===")

    except Exception as e:
        logging.debug(f"Process tree debug failed: {e}")


def _is_in_process_tree(target_pid, root_pid, all_processes):
    """
    Check if a process is in the process tree starting from root_pid.
    """
    if target_pid == root_pid:
        return True

    # Find the process
    target_proc = None
    for proc_info in all_processes:
        if proc_info["pid"] == target_pid:
            target_proc = proc_info
            break

    if not target_proc:
        return False

    # Recursively check if parent is in the tree
    return _is_in_process_tree(target_proc["ppid"], root_pid, all_processes)


def _explain_process_inclusion(proc_info, main_pid):
    """
    Explain why a process is being included in our child process count.
    """
    pid = proc_info["pid"]
    ppid = proc_info["ppid"]
    name = proc_info["name"]
    cmdline = proc_info["cmdline"]

    reasons = []

    # Direct child
    if ppid == main_pid:
        reasons.append("Direct child of main ROBIN process")

    # Python process with ROBIN in command line
    if name in ["python", "python3", "python3.9", "robin"]:
        if cmdline and any("robin" in arg.lower() for arg in cmdline):
            reasons.append("Python process with ROBIN in command line")

    # Analysis tool
    if name in ["modkit", "matkit", "sturgeon", "R", "Rscript", "bedtools"]:
        reasons.append("Analysis tool spawned by ROBIN")

    # Process in tree but not direct child
    if ppid != main_pid and ppid != 1:
        reasons.append("Process in ROBIN's process tree (indirect child)")

    # Reparented process
    if ppid == 1:
        reasons.append("Reparented to init/systemd (orphaned process)")

    # Shared library or dependency
    if name in ["libc", "libpython", "libssl", "libcrypto"]:
        reasons.append("Shared library process")

    if not reasons:
        reasons.append("Unknown reason - process in tree but unclear relationship")

    return "; ".join(reasons)
