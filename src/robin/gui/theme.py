"""
Module: theme

This module defines the theme and layout for the entire application, implementing Material Design 3 (M3) principles
for a modern, accessible, and consistent user experience across all pages.

It includes:

- A context manager `frame` to create a custom page frame with navigation, header, and footer.
- Utility functions to handle dark mode and remote access toggling.
- Material Design 3 color tokens, typography scale, and spacing system.
- Modern component styling with proper elevation and surface treatments.
- Responsive design patterns and accessibility improvements.

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

# These will be set by the get_imagefile() and get_version() functions
IMAGEFILE = None
__about__ = None


import os
import psutil
import platform

# Check if we're in development mode
is_development_mode = os.environ.get("ROBIN_DEV_MODE", "").lower() in ("1", "true", "yes", "on")


def get_imagefile():
    """Get the path to the ROBIN logo image file."""
    global IMAGEFILE
    if IMAGEFILE is None:
        try:
            from robin.gui import images
            IMAGEFILE = os.path.join(
                os.path.dirname(os.path.abspath(images.__file__)), "ROBIN_logo_small.png"
            )
        except (ImportError, AttributeError):
            # Fallback path when running standalone
            IMAGEFILE = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "images", "ROBIN_logo_small.png"
            )
    return IMAGEFILE


def get_about():
    """Get the __about__ module with version information."""
    global __about__
    if __about__ is None:
        try:
            from robin import __about__
        except ImportError:
            # Fallback when running standalone - create a minimal __about__ object
            class MockAbout:
                __version__ = "standalone-test"
            __about__ = MockAbout()
    return __about__


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


def styled_table(*, columns, rows=None, pagination=20, class_size="table-xs", **kwargs):
    """Create a NiceGUI table with Material Design 3 styling and compact layout.

    Args:
        columns: columns definition passed to ui.table
        rows: initial rows
        pagination: rows per page (0 disables)
        class_size: table size class (e.g., "table-xs", "table-sm")
        **kwargs: forwarded to ui.table

    Returns:
        Tuple of (container, table) where container is the overflow wrapper column and table is the ui.table instance.
    """
    # Outer container with M3 styling
    container = ui.column().classes(
        "w-full overflow-x-auto compact-table elevation-1 rounded-lg"
    )
    with container:
        table = ui.table(
            columns=columns, rows=rows or [], pagination=pagination, **kwargs
        )
        try:
            table.classes(replace=f"table w-full {class_size} text-xs")
        except Exception:
            table.classes(f"table w-full {class_size} text-xs")
        # Use Quasar's dense mode with M3 styling for maximum compactness
        try:
            table.props("dense flat wrap-cells")
        except Exception:
            pass
    return container, table


async def check_version():
    """
    Check the current version against the remote version on GitHub.
    Shows a notification or dialog to the user about their version status.
    """
    # Skip version check in development mode
    if is_development_mode:
        return
        
    # Check if version has already been checked in this app session
    try:
        if app.storage.general.get("version_checked", False):
            return
    except RuntimeError:
        # Storage not available in this context, continue with version check
        pass

    try:
        remote_version_str = await run.io_bound(get_version_from_github)

        if not remote_version_str:
            with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
                ui.label("Version Check Error").classes(
                    "text-headline-small text-weight-bold q-mb-md"
                )
                ui.label(
                    "Could not determine remote version. Please check manually."
                ).classes("text-body-medium q-mb-md")
                ui.button("OK", on_click=dialog.close).classes(
                    "bg-primary text-white rounded-md"
                )
            dialog.open()
            return

        local_version = version.parse(get_about().__version__)
        remote_version = version.parse(remote_version_str)

        if local_version == remote_version:
            ui.notify("Your ROBIN installation is up to date!", type="positive")
        elif local_version < remote_version:
            with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
                ui.label("Update Available!").classes(
                    "text-headline-small text-weight-bold q-mb-md"
                )
                ui.label(f"Your version: {local_version}").classes("text-body-medium")
                ui.label(f"Latest version: {remote_version}").classes(
                    "text-body-medium"
                )
                ui.label(
                    "Would you like to visit the GitHub repository to update?"
                ).classes("text-body-medium q-mb-md")
                with ui.row().classes("w-full justify-center gap-3"):
                    ui.button(
                        "Continue with current version", on_click=dialog.close
                    ).classes("rounded-md")
                    ui.button(
                        "Visit GitHub",
                        on_click=lambda: ui.open("https://github.com/LooseLab/ROBIN"),
                    ).classes("bg-primary text-white rounded-md elevation-2")
            dialog.open()
        else:
            with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
                ui.label("Development Version").classes(
                    "text-headline-small text-weight-bold q-mb-md"
                )
                ui.label(
                    f"You are running a development version ({local_version})."
                ).classes("text-body-medium")
                ui.label(f"Latest release: {remote_version}").classes(
                    "text-body-medium"
                )
                ui.label(
                    "This version may be unstable and is only for testing purposes. It is not recommended for production use."
                ).classes("text-body-medium q-mb-md")
                ui.label("Please consider using the latest release instead.").classes(
                    "text-body-medium"
                )
                ui.button("OK", on_click=dialog.close).classes(
                    "bg-primary text-white rounded-md"
                )
            dialog.open()

    except requests.RequestException:
        with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
            ui.label("Connection Error").classes(
                "text-headline-small text-weight-bold q-mb-md"
            )
            ui.label("Could not check for updates.").classes("text-body-medium")
            ui.label(
                "Either you are not connected to the internet or you cannot access https://www.github.com/looselab/robin."
            ).classes("text-body-medium q-mb-md")
            ui.label("Please manually check for updates.").classes("text-body-medium")
            ui.button("OK", on_click=dialog.close).classes(
                "bg-primary text-white rounded-md"
            )
        dialog.open()
    except Exception as e:
        with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
            ui.label("Error").classes("text-headline-small text-weight-bold q-mb-md")
            ui.label(f"Error checking version: {str(e)}").classes("text-body-medium")
            ui.button("OK", on_click=dialog.close).classes(
                "bg-primary text-white rounded-md"
            )
        dialog.open()

    # Mark version as checked for this app session
    try:
        app.storage.general["version_checked"] = True
    except RuntimeError:
        # Storage not available in this context, skip setting the flag
        pass


# IMAGEFILE is now defined above with fallback handling

# Module-level variables
quitdialog = None

MENU_BREAKPOINT = 1200


class GlobalSystemMetrics:
    """Global system metrics singleton that provides CPU and RAM usage data."""
    
    _instance = None
    _timer_active = False
    
    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.cpu = 0
            cls._instance.ram = 0
        return cls._instance
    
    def start_timer(self):
        """Start the global metrics timer if not already running."""
        if not self._timer_active:
            self._timer_active = True
            ui.timer(1.0, self.update_metrics)
    
    def update_metrics(self):
        """Update CPU and RAM metrics."""
        try:
            self.cpu = round(psutil.getloadavg()[1] / os.cpu_count() * 100, 1)
            self.ram = round(psutil.virtual_memory()[2], 1)
        except (OSError, AttributeError):
            # Handle cases where psutil might fail or os.cpu_count() returns None
            self.cpu = 0
            self.ram = 0


# Global metrics instance
global_metrics = GlobalSystemMetrics()

# Read the HTML content for the header
HEADER_HTML = (Path(__file__).parent / "static" / "header.html").read_text()

# Read the CSS styles for the application
STYLE_CSS = (Path(__file__).parent / "static" / "styles.css").read_text()
M3_COMPONENTS_CSS = (Path(__file__).parent / "static" / "m3-components.css").read_text()
MOSAIC_COMPONENTS_CSS = (Path(__file__).parent / "static" / "mosaic-components.css").read_text()


@contextmanager
def frame(navtitle: str, batphone=False, smalltitle=None, center: str = None):
    """
    Context manager to create a custom page frame with Material Design 3 styling and consistent behavior across all pages.

    Args:
        navtitle (str): The title to display in the navigation header.
        batphone (bool): Whether to show the BATMAN mode title.
        smalltitle (str): The title to display on small screens.
        center (str): Center ID running the analysis.

    Yields:
        None
    """
    global quitdialog
    if batphone:
        navtitle = f"BATMAN & {navtitle}"

    # Store center in app storage if provided
    if center:
        app.storage.general["center"] = center

    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.2.0/dist/igv.min.js"></script>'
    )
    ui.add_head_html(
        HEADER_HTML + f"<style>{STYLE_CSS}</style><style>{M3_COMPONENTS_CSS}</style><style>{MOSAIC_COMPONENTS_CSS}</style>"
    )
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

    # Create disclaimer dialog that appears on first visit with M3 styling
    async def show_disclaimer():
        # Skip disclaimer in development mode
        if is_development_mode:
            return
            
        try:
            disclaimer_acknowledged = app.storage.general.get("disclaimer_acknowledged", False)
        except RuntimeError:
            # Storage not available in this context, show disclaimer
            disclaimer_acknowledged = False
        
        if not disclaimer_acknowledged:
            with ui.dialog().props(
                "persistent"
            ) as disclaimer_dialog, ui.card().classes("w-160 elevation-4 rounded-xl"):
                ui.label("DISCLAIMER").classes(
                    "text-headline-small text-weight-bold q-mb-md"
                )
                ui.label(
                    "This tool and the data generated by it are intended for research use only and should not be used for "
                    "direct diagnostic purposes. The methylation-based classifications and other analyses provided here may "
                    "be considered by neuropathologists as supplementary information in the context of comprehensive "
                    "diagnostic assessment, which should include clinical history, radiological findings, and complete "
                    "histopathological and molecular evaluation. The final interpretation and diagnosis should always be "
                    "made by qualified healthcare professionals based on all available information."
                ).classes("text-body-medium q-mb-md")

                if batphone:
                    ui.label("BATMAN Mode").classes(
                        "text-headline-small text-weight-bold q-mb-md"
                    )
                    ui.label(
                        "You are running this tool in BATMAN mode. "
                        "This is a beta version of the tool and may not be fully functional. "
                        "BATMAN means: Breakpoint Adaptive Targeting alongside Methylation Analysis on Nanopore. "
                        "This means that the target regions will be updated in real-time based on detected breakpoints. "
                        "This code only works with ReadFish at this time. "
                    ).classes("text-body-medium")

                def acknowledge():
                    try:
                        app.storage.general["disclaimer_acknowledged"] = True
                    except RuntimeError:
                        # Storage not available in this context, skip setting the flag
                        pass
                    disclaimer_dialog.close()

                ui.button("I agree", on_click=acknowledge).classes(
                    "bg-primary text-white rounded-md elevation-2"
                )
            disclaimer_dialog.open()

    ui.timer(0.5, show_disclaimer, once=True)

    # Add version check timer
    ui.timer(1.0, check_version, once=True)

    # Create a persistent dialog for quitting the app with M3 styling
    quitdialog = ui.dialog().props("persistent")

    async def quit_app():
        quitdialog.close()
        await cleanup_and_exit()

    with quitdialog, ui.card().classes("elevation-4 rounded-xl"):
        ui.label(
            "Quitting the app will stop running methylation analysis. Are you sure?"
        ).classes("text-headline-small text-weight-bold q-mb-md")
        ui.label("If you want to keep analysis running, click Cancel.").classes(
            "text-body-medium"
        )
        ui.label(
            "You can safely close this window and analysis will keep running in the background."
        ).classes("text-body-medium q-mb-md")
        with ui.row().classes("w-full justify-center gap-3"):
            ui.button("Cancel", on_click=quitdialog.close).classes("rounded-md")
            ui.button("Really Quit", icon="logout", on_click=quit_app).classes(
                "bg-error text-white rounded-md elevation-2"
            )

    # Create a header with navigation title and menu using M3 styling
    header_classes = "items-center duration-200 p-0 px-4 no-wrap elevation-1"
    if batphone:
        header_classes += " batphone"

    with ui.header(elevated=True).classes(header_classes):
        with ui.grid(columns=2).style("width: 100%"):
            with ui.row().classes(
                f"max-[{MENU_BREAKPOINT}px]:hidden items-center align-left px-4"
            ):
                ui.html(navtitle, sanitize=False).classes("text-headline-medium drop-shadow font-bold").style(
                    "font-weight: 600; font-family: var(--font-primary)"
                )
            with ui.row().classes(
                f"min-[{MENU_BREAKPOINT+1}px]:hidden items-center align-left px-4"
            ):
                ui.html(smalltitle, sanitize=False).classes("text-headline-medium drop-shadow font-bold").style(
                    "font-weight: 600; font-family: var(--font-primary)"
                )
            with ui.row().classes("ml-auto align-top"):
                with ui.row().classes("items-center m-auto gap-3"):
                    ui.label(f"Viewing: {platform.node()}").classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden text-body-medium"
                    )
                    ui.label("CPU").classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden text-body-medium"
                    )
                    cpu_activity = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )
                    ui.label("RAM").classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden text-body-medium"
                    )
                    ram_utilisation = ui.circular_progress(max=100).classes(
                        f"max-[{MENU_BREAKPOINT}px]:hidden"
                    )

                    # Start global metrics timer if not already running
                    global_metrics.start_timer()

                    # Bind the progress indicators to the global metrics
                    cpu_activity.bind_value_from(global_metrics, "cpu")
                    ram_utilisation.bind_value_from(global_metrics, "ram")

                    with ui.button(icon="menu").classes("rounded-md"):
                        with ui.menu() as menu:
                            ui.menu_item("Home", lambda: ui.navigate.to("/")).classes(
                                "text-body-medium"
                            )
                            ui.menu_item(
                                "View Samples", lambda: ui.navigate.to("/live_data")
                            ).classes("text-body-medium")
                            ui.menu_item(
                                "Activity Monitor",
                                lambda: ui.navigate.to("/robin"),
                            ).classes("text-body-medium")
                            ui.menu_item(
                                "Workflow",
                                lambda: ui.navigate.to("/workflow"),
                            ).classes("text-body-medium")
                            ui.menu_item(
                                "Documentation",
                                lambda: ui.navigate.to(
                                    "https://looselab.github.io/ROBIN/"
                                ),
                            ).classes("text-body-medium")
                            ui.separator()
                            ui.switch("Allow Remote Access").classes(
                                "ml-4 bg-transparent"
                            ).props('color="primary"').bind_value(
                                app.storage.general, "use_on_air"
                            )
                            ui.separator()
                            ui.switch("Dark Mode").classes("ml-4 bg-transparent").props(
                                'color="primary"'
                            ).bind_value(app.storage.browser, "dark_mode")
                            ui.dark_mode().bind_value(app.storage.browser, "dark_mode")
                            ui.separator()
                            ui.menu_item("Close", menu.close).classes(
                                "text-body-medium"
                            )
                            ui.button(
                                "Quit", icon="logout", on_click=quitdialog.open
                            ).classes("bg-error text-white rounded-md")
                    ui.image(get_imagefile()).style("width: 50px")

    with ui.column().classes(
        "w-full h-full max-w-full overflow-hidden"
    ) as main_content:
        pass

    # Create a footer with useful information and quit button using M3 styling
    footer_classes = "items-center elevation-1"
    if batphone:
        footer_classes += " batphone"
    with ui.footer().classes(footer_classes):
        with ui.dialog() as dialog, ui.card().classes("elevation-3 rounded-xl"):
            ui.label("Links").classes("text-headline-small text-weight-bold q-mb-md")
            ui.separator()
            ui.link("Code on GitHub", "https://github.com/looselab/robin").classes(
                "text-body-medium"
            )
            ui.link(
                "Rapid CNS2 Paper",
                "https://link.springer.com/article/10.1007/s00401-022-02415-6",
            ).classes("text-body-medium")
            ui.link(
                "Sturgeon Classifier",
                "https://www.nature.com/articles/s41586-023-06615-2",
            ).classes("text-body-medium")
            ui.link(
                "Protocol",
                "https://www.protocols.io/view/intra-operative-nanopore-sequencing-to-classify-br-c65qzg5w",
            ).classes("text-body-medium")
            ui.link("Oxford Nanopore", "https://nanoporetech.com/").classes(
                "text-body-medium"
            )
            ui.link("epi2me labs", "https://labs.epi2me.io/").classes(
                "text-body-medium"
            )
            ui.link("Looselab", "https://looselab.github.io/").classes(
                "text-body-medium"
            )
            ui.button("Close", on_click=dialog.close).classes(
                "bg-primary text-white rounded-md"
            )
        ui.image(get_imagefile()).style("width: 40px")
        ui.colors(primary="#4F9153")  # Green primary color
        ui.button("Links", on_click=dialog.open).classes("rounded-md")

        with ui.button(icon="info").classes("rounded-md"):
            with ui.menu() as menu:
                ui.label().bind_text_from(
                    app, "urls", backward=lambda n: f"Available urls: {n}"
                ).classes("text-body-medium")
                ui.label("Version: " + get_about().__version__).classes(
                    "text-body-medium"
                )
        ui.label(
            "Some aspects of this application are ©Looselab - all analyses provided for research use only."
        ).classes(
            f"max-[{MENU_BREAKPOINT}px]:hidden text-body-small text-weight-italic"
        )
        ui.label("©Looselab").classes(
            f"min-[{MENU_BREAKPOINT+1}px]:hidden text-body-small text-weight-italic"
        )
        ui.label("Not for diagnostic use.").classes(
            f"min-[{MENU_BREAKPOINT+1}px]:hidden text-body-small text-weight-italic"
        )
        ui.label(f"Center: {app.storage.general.get('center', 'Not set')}").classes(
            f"max-[{MENU_BREAKPOINT}px]:hidden text-body-small text-weight-medium"
        )

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

    # Create and show shutdown modal with M3 styling
    with ui.dialog().props("persistent") as shutdown_dialog, ui.card().classes(
        "w-96 elevation-4 rounded-xl"
    ):
        ui.label("Shutting Down").classes(
            "text-headline-small text-weight-bold q-mb-md"
        )
        ui.label("ROBIN is shutting down. Please wait while we clean up...").classes(
            "text-body-medium q-mb-md"
        )
        with ui.row().classes("w-full justify-center"):
            ui.spinner(size="lg", color="primary")

    # Close the quit dialog if it's open
    if quitdialog:
        quitdialog.close()

    # Show the shutdown dialog
    shutdown_dialog.open()
    logging.info("Shutdown dialog opened")

    # Wait a moment to ensure the dialog is visible
    await asyncio.sleep(1.0)
    logging.info("Initial wait complete")

    # Perform cleanup
    logging.info("Performing cleanup operations...")
    logging.info("Shutting down ROBIN... from theme.py")
    logging.info(
        "Here we need to do some very graceful shutdown to make sure we don't leave any threads running and we don't leave any files open."
    )

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


def create_home_page():
    """Create the home page content."""
    with frame(
        "<strong>R</strong>apid nanop<strong>O</strong>re <strong>B</strong>rain intraoperat<strong>I</strong>ve classificatio<strong>N</strong>",
        smalltitle="<strong>R.O.B.I.N</strong>",
    ):
        
        ui.label("Welcome to the Application").classes(
            "text-headline-large text-center"
        )
        
        with ui.card().classes("w-full").style("border: 2px solid var(--md-primary)"):
            with ui.row().classes("w-full flex justify-between items-center"):
                ui.label('Sample Name').classes("text-headline-medium ")
                ui.button("Button").classes("bg-primary text-white rounded-md")
            ui.separator().classes().style("border: 1px solid var(--md-primary)")
            with ui.card().classes("w-full bg-gradient-to-r from-blue-50 to-indigo-50"):
                ui.label("Run Information").classes("text-lg font-semibold mb-3 text-blue-800")
                with ui.row().classes("w-full gap-8 items-center flex-wrap"):
                    ui.label("Run").classes("text-body-medium")
                    ui.label("Model").classes("text-body-medium")
                    ui.label("Device").classes("text-body-medium")
                    ui.label("Flow Cell").classes("text-body-medium")
                    ui.label("Sample").classes("text-body-medium")


            with ui.card().classes("w-full"):
                ui.label("Classification Results").classes("text-lg font-semibold mb-3 text-blue-800")
                with ui.row().classes("w-full flex justify-between items-center"):
                    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-blue-50 to-indigo-50 border-l-4 border-blue-500"):
                        ui.label("Sturgeon Classification").classes("font-bold text-blue-800 mb-2")
                        ui.label("Class: --").classes("font-bold text-medium text-blue-600")
                        ui.label("Confidence: --%").classes("text-sm text-blue-600")
                        ui.label("Probes: --").classes("text-sm text-blue-600")
                        ui.label("Model: --").classes("text-sm text-blue-600")
                        ui.label("Features: --").classes("text-sm text-blue-600")

                    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-green-50 to-green-100 border-l-4 border-green-500"):
                        ui.label("NanoDX Classification").classes("font-bold text-green-800 mb-2")
                        ui.label("Class: --").classes("font-bold text-medium text-green-600")
                        ui.label("Confidence: --%").classes("text-sm text-green-600")
                        ui.label("Probes: --").classes("text-sm text-green-600")
                        ui.label("Model: --").classes("text-sm text-blue-600")
                        ui.label("Features: --").classes("text-sm text-blue-600")

                    
                    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-purple-50 to-purple-100 border-l-4 border-purple-500"):
                        ui.label("PanNanoDX Classification").classes("font-bold text-purple-800 mb-2")
                        ui.label("Class: --").classes("font-bold text-medium text-purple-600")
                        ui.label("Confidence: --%").classes("text-sm text-purple-600")
                        ui.label("Probes: --").classes("text-sm text-purple-600")
                        ui.label("Model: --").classes("text-sm text-blue-600")
                        ui.label("Features: --").classes("text-sm text-blue-600")

                    
                    with ui.card().classes("flex-1 elevation-4 rounded-xl bg-gradient-to-br from-orange-50 to-orange-100 border-l-4 border-orange-500"):
                        ui.label("Random Forest Classification").classes("font-bold text-orange-800 mb-2")
                        ui.label("Class: --").classes("font-bold text-medium text-orange-600")
                        ui.label("Confidence: --%").classes("text-sm text-orange-600")
                        ui.label("Probes: --").classes("text-sm text-orange-600")
                        ui.label("Model: --").classes("text-sm text-blue-600")
                        ui.label("Features: --").classes("text-sm text-blue-600")


            ui.label('text below').classes("text-body-medium")
            with ui.card_section():
                ui.image('https://picsum.photos/id/684/640/360').classes()
                ui.label('Lorem ipsum dolor sit amet, consectetur adipiscing elit, ...')


def create_standalone_page():
    """Create a simplified standalone page without storage dependencies."""
    # Add custom HTML and CSS to the head of the page
    ui.add_head_html(
        '<script src="https://cdn.jsdelivr.net/npm/igv@3.2.0/dist/igv.min.js"></script>'
    )
    ui.add_head_html(
        HEADER_HTML + f"<style>{STYLE_CSS}</style><style>{M3_COMPONENTS_CSS}</style><style>{MOSAIC_COMPONENTS_CSS}</style>"
    )
    
    # Create a simple header
    with ui.header(elevated=True).classes("items-center duration-200 p-0 px-4 no-wrap elevation-1"):
        ui.html("<strong>R.O.B.I.N</strong>", sanitize=False).classes("text-headline-medium drop-shadow font-bold").style(
            "font-weight: 600; font-family: var(--font-primary)"
        )
        ui.image(get_imagefile()).style("width: 50px").classes("ml-auto")
    
    # Create main content
    with ui.column().classes("w-full h-full max-w-full overflow-hidden p-8"):
        ui.label("Welcome to ROBIN Theme Test").classes("text-headline-large text-center")
        ui.label("This is a standalone test of the ROBIN theme system.").classes("text-body-large text-center mt-4")
        ui.label(f"Version: {get_about().__version__}").classes("text-body-medium text-center mt-2")
    
    # Create a simple footer
    with ui.footer().classes("items-center elevation-1"):
        ui.image(get_imagefile()).style("width: 40px")
        ui.label("ROBIN Theme Test - Standalone Mode").classes("text-body-small")


def create_workflow_page():
    """Display the ROBIN workflow diagram with M3 styling."""
    with frame(
        "ROBIN Workflow",
        smalltitle="Workflow",
    ):
        # Get versions for the diagram
        sturgeon_version = get_sturgeon_version()
        crossnn_version = get_crossnn_version()
        cnv_from_bam_version = get_cnv_from_bam_version()

        ui.mermaid(
            f"""
flowchart TD
    %% Style definitions with M3 color palette
    classDef minKNOW fill:#E8DEF8,stroke:#6750A4,stroke-width:2px,color:#1D192B,font-size:14px,font-weight:500
    classDef robin fill:#EADDFF,stroke:#6750A4,stroke-width:2px,color:#21005D,font-size:14px,font-weight:500
    classDef classifier fill:#FFD8E4,stroke:#7D5260,stroke-width:2px,color:#31111D,font-size:14px,font-weight:500
    classDef analysis fill:#FFD8E4,stroke:#7D5260,stroke-width:2px,color:#31111D,font-size:14px,font-weight:500
    classDef output fill:#EADDFF,stroke:#6750A4,stroke-width:2px,color:#21005D,font-size:14px,font-weight:500
    classDef headerLabel fill:#ffffff00,stroke:#ffffff00,color:#1C1B1F,font-size:18px,font-weight:700

    %% Custom header nodes
    robinLabel["<b>R.O.B.I.N v{get_about().__version__}</b>"]:::headerLabel
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
    %% Subgraph background coloring with M3 colors
    style MinKNOW fill:#E8DEF8,stroke:#6750A4,stroke-width:2px
    style ROBIN fill:#EADDFF,stroke:#6750A4,stroke-width:2px
""",
            config={
                "theme": "redux",
                "look": "neo",
                "flowchart": {"curve": "basis", "defaultRenderer": "elk"},
            },
        ).classes("w-full elevation-2 rounded-xl")


def register_theme_pages():
    """Register the theme pages. This function should be called when the module is imported."""
    @ui.page("/")
    def home_page():
        create_home_page()

    @ui.page("/workflow")
    def workflow_page():
        create_workflow_page()


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
    # Add Apple HIG + Material Design 3 fonts
    ui.add_head_html(
        """
        <link rel="preconnect" href="https://fonts.googleapis.com">
        <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
        <link href="https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
        <style>
            /* Apple HIG: System font fallbacks */
            :root {
                --font-primary: "Inter", "SF Pro Display", "SF Pro", -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
                --font-mono: "JetBrains Mono", "SF Mono", "Monaco", "Menlo", "Consolas", "Liberation Mono", "Courier New", monospace;
            }
        </style>
    """
    )

    # Register some fonts that we might need later on (guarded)
    fonts_dir = Path(__file__).parent / "fonts"
    if fonts_dir.exists():
        app.add_static_files("/fonts", str(fonts_dir))

    # For standalone execution, create UI directly without @ui.page decorators
    # This completely avoids the global scope issue with NiceGUI
    create_standalone_page()

    ui.run(storage_secret="robin")


# Register theme pages when module is imported (but not when run as main)
if __name__ not in {"__main__", "__mp_main__"}:
    register_theme_pages()

if __name__ in {"__main__", "__mp_main__"}:
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
            except Exception as e:
                logging.debug(f"   Memory: <access denied>: {e}")
                
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
