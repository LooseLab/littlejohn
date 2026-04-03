# What happens when you start ROBIN

This page describes the **typical experience** when you run **`robin workflow`** from a terminal: checks, prompts, the **research disclaimer**, **GUI password**, and where the **NiceGUI** monitor appears. Behaviour matches the current codebase; always check **`robin workflow --help`** for your installed version.

For install steps, see [Installation](installation.md). For command-line flags, see the [command-line reference](../cli/workflow.md).

---

## 1. Model assets

Before anything else, ROBIN verifies that **required model files** are present (same manifest as `robin utils update-models`). If files are missing, the process exits with a message telling you to run **`robin utils update-models`** (and set **`GITHUB_TOKEN`** if your assets are on private GitHub).

---

## 2. Reference genome (if you pass `--reference` / `-r`)

If you supply a reference **FASTA**, ROBIN validates it and ensures an index (e.g. `.fai`) can be created or found. On failure, it exits with an error — fix the path or omit `--reference` only if your workflow truly does not need it.

---

## 3. Research disclaimer (`I agree`)

ROBIN prints the **research-use disclaimer** and waits for you to type **`I agree`** exactly (case-sensitive). This is required before the workflow starts.

---

## 4. Optional warning: large BAMs

If **`ROBIN_PROCESS_LARGE_BAMS`** is enabled, ROBIN prints a **warning** that this mode must **not** be used together with live sequencing runs.

---

## 5. Configuration summary

The terminal prints a summary of your run: paths, **`--center`**, workflow steps (including automatic **`preprocessing`** / **`bed_conversion`** insertion), logging, Ray/threading mode, etc.

---

## 6. Execution engine (Ray vs threading)

- **Default (`--use-ray`)** — Ray is initialised (with optional dashboard, CPU limits, presets). The **Ray Core** workflow driver runs asynchronously.
- **`--no-use-ray`** — Falls back to threaded workers instead of Ray.

---

## 7. NiceGUI workflow monitor (default: on)

With **`--with-gui`** (the default), ROBIN starts a **web-based** workflow monitor (NiceGUI) so you can follow progress in a browser.

### When the GUI actually starts

| Mode | GUI behaviour |
|------|----------------|
| **Ray workflow (default)** | The GUI is started **inside** the Ray workflow driver **only if** you pass **`--work-dir`** (`-d`). If `--work-dir` is omitted, the driver may skip launching the GUI and print that **`--work-dir` was not provided**. |
| **`--no-use-ray`** | The CLI launches the GUI when **`--with-gui`** is set, using **`--work-dir`** if provided, otherwise the **watched BAM directory** as the monitored path. |

So for the **default Ray** path, plan to pass **both** a data directory (positional `PATH`) **and** **`--work-dir`** if you want the browser UI.

### URLs

When the GUI starts, the terminal prints a base URL such as **`http://<gui-host>:<gui-port>`** (defaults: host **`0.0.0.0`**, port **`8081`**). Typical entry points include:

- Welcome / root: `/`
- Workflow monitor: `/robin`
- Sample-oriented views: paths under `/live_data` (exact routes are printed at startup)

Use **`--gui-host`** and **`--gui-port`** to change bind address and port; use **`--no-gui`** to disable the web UI entirely.

### Disable the GUI

```bash
robin workflow ... --no-gui
```

---

## 8. GUI password (terminal prompts)

Access to the web UI is protected by a **password** (requires **`argon2-cffi`** for secure handling). ROBIN prompts in the terminal; passwords are **not** echoed.

### First run (password not set yet)

If no password has been set and **stdin is a TTY** (interactive terminal), ROBIN prompts:

```text
Set GUI password:
Confirm GUI password:
```

You must enter the same password twice. The password is stored for future logins.

If **no password has been set yet** and **stdin is not a TTY** (e.g. some automated contexts), startup **fails** with a message to run ROBIN from a terminal so the password can be set interactively.

### Later runs (password already set)

If stdin is a **TTY**, you are prompted once:

```text
GUI password:
```

Enter the same password you set earlier. Wrong password → **Invalid password.** and the GUI does not start.

### Setting or changing the password without a full workflow

Use the dedicated command (see [GUI password](../cli/password.md)):

```bash
robin password set
```

This can **replace** an existing password after confirmation. It uses the same mechanism as the first-run prompts.

---

## 9. Workflow hooks (optional message)

If the GUI starts successfully, ROBIN may install **workflow hooks** for live updates. If hook installation fails, you may see a message that the GUI will show **static** information only.

---

## 10. Watching BAMs and shutting down

After the above, the runner **watches** the input directory for **`*.bam`** files (subject to `--no-process-existing`, `--no-watch`, etc.) and schedules jobs. Use **Ctrl+C** to stop; ROBIN attempts a **graceful shutdown** (workflow manager, Ray, GUI), though complex runs may take a moment to exit.

---

## Quick reference — order of prompts

1. (Automatic) Model check  
2. (If `-r`) Reference validation  
3. **Type `I agree`** — research disclaimer  
4. (If GUI enabled) **GUI password** — set twice first time, or single verify later  
5. Browser → open printed URL to monitor  

---

## See also

- [Quickstart](quickstart.md)  
- [`robin workflow`](../cli/workflow.md)  
- [GUI password command](../cli/password.md)  
- [CLI overview](../cli/index.md)  
