"""
Local folder picker based on NiceGUI's local_file_picker example.

Allows selecting a directory from the local filesystem where NiceGUI is running.
Uses an in-dialog browser - no native window or new GUI instance.
"""

import platform
from pathlib import Path
from typing import Optional

from nicegui import ui


class local_folder_picker(ui.dialog):
    """Folder picker dialog for selecting a directory on the server filesystem."""

    def __init__(
        self,
        directory: str,
        *,
        upper_limit: Optional[str] = ...,
        show_hidden_files: bool = False,
    ) -> None:
        """Local folder picker.

        :param directory: The directory to start in.
        :param upper_limit: The directory to stop at (None: no limit, default: same as starting directory).
            Use None to allow full filesystem navigation.
        :param show_hidden_files: Whether to show hidden files/folders.
        """
        super().__init__()

        self.path = Path(directory).expanduser()
        if not self.path.exists():
            self.path = Path.home()
        if not self.path.is_dir():
            self.path = self.path.parent

        if upper_limit is None:
            self.upper_limit = None
        else:
            self.upper_limit = Path(
                directory if upper_limit == ... else upper_limit
            ).expanduser()
        self.show_hidden_files = show_hidden_files

        with self, ui.card().classes(
            "robin-dialog-surface p-4 w-full max-w-md min-w-[18rem]"
        ):
            self.add_drives_toggle()
            self.path_label = ui.label().classes(
                "text-sm workflow-folder-picker-path mb-2"
            )
            self.list_container = ui.column().classes(
                "w-full max-w-md max-h-80 overflow-y-auto min-w-0"
            )
            with ui.row().classes("w-full justify-end gap-2 mt-2 flex-wrap"):
                ui.button("Cancel", on_click=self.close).props("flat no-caps outline")
                ui.button(
                    "Select this folder",
                    on_click=lambda: self.submit([str(self.path)]),
                    icon="check",
                ).props("color=primary no-caps")
            self.update_list()

    def add_drives_toggle(self) -> None:
        """Add drive selector on Windows (optional, may fail without pywin32)."""
        if platform.system() == "Windows":
            try:
                import win32api

                drives = win32api.GetLogicalDriveStrings().split("\000")[:-1]
                if drives:
                    self.drives_toggle = ui.toggle(
                        drives, value=drives[0], on_change=self.update_drive
                    )
            except ImportError:
                pass

    def update_drive(self) -> None:
        """Handle drive change on Windows."""
        if hasattr(self, "drives_toggle"):
            self.path = Path(self.drives_toggle.value).expanduser()
            self.update_list()

    def update_list(self) -> None:
        """Refresh the list with current directory contents."""
        self.path_label.set_text(str(self.path))

        self.list_container.clear()
        with self.list_container:
            try:
                paths = list(self.path.glob("*"))
            except (OSError, PermissionError):
                paths = []

            if not self.show_hidden_files:
                paths = [p for p in paths if not p.name.startswith(".")]
            paths.sort(key=lambda p: p.name.lower())
            paths.sort(key=lambda p: not p.is_dir())

            # Parent directory (..) - single click to go up (always show when not at root)
            if self.path != self.path.parent and (
                self.upper_limit is None or self.path != self.upper_limit
            ):
                parent_path = self.path.parent

                def go_up(_=None):
                    self.path = parent_path
                    self.update_list()

                with ui.row().classes(
                    "w-full items-center gap-2 p-2 cursor-pointer rounded "
                    "workflow-folder-picker-row"
                ).on("click", go_up):
                    ui.icon("arrow_upward").classes("text-base shrink-0")
                    ui.label("..").classes("font-medium")

            for p in paths:
                path_str = str(p)

                def make_click_handler(target):
                    def handler(_=None):
                        target_path = Path(target)
                        if target_path.is_dir():
                            self.path = target_path
                            self.update_list()

                    return handler

                if p.is_dir():
                    with ui.row().classes(
                        "w-full items-center gap-2 p-2 cursor-pointer rounded "
                        "workflow-folder-picker-row"
                    ).on("click", make_click_handler(path_str)):
                        ui.icon("folder").classes("text-base shrink-0")
                        ui.label(p.name).classes("font-medium")
                else:
                    with ui.row().classes(
                        "w-full items-center gap-2 p-2 workflow-folder-picker-file"
                    ):
                        ui.icon("insert_drive_file").classes(
                            "text-base shrink-0 opacity-70"
                        )
                        ui.label(p.name)
