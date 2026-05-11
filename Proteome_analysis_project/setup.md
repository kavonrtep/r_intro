# Setup

This page covers the tools needed for the project on **Ubuntu** (tested on 24.04 LTS):

1. A coding agent — either **GitHub Copilot CLI** (requires a Copilot subscription) or **OpenCode** (works with free model providers).
2. **VS Code** as a general-purpose editor. This tool is not strictly necessary, but it provides more versatility compared to RStudio.

If you already have a working coding agent and editor, you can skip this page.

---

# Coding agent

Pick **one** of the two options below. Both run in the terminal and let you interact with an LLM that can read, edit, and run code in the current directory.

## Option A — GitHub Copilot CLI

The official agentic CLI from GitHub. Requires:

- An active **Copilot subscription** (Pro, Pro+, Business, or Enterprise — student Pro accounts also work).

```bash
curl -fsSL https://gh.io/copilot-install | bash
```

### Run and log in

```bash
copilot
```
Note: In some cases the copilot command may not be immediately available after installation. You can run it by specifying the full path, which is usually `~/.local/bin/copilot`:

```bash
~/.local/bin/copilot
```

On the first launch of `copilot`, you will be prompted for authentication. Type:

```
/login
```

inside the TUI — this opens a browser window for the standard GitHub OAuth/device-code flow. Sign in with the GitHub account that holds your Copilot subscription.

Useful slash commands once inside:
- `/login` — authenticate
- `/model` — switch model (default is Claude Sonnet 4.5)
- `/help` — list all commands
- `/exit` — quit

Docs: <https://docs.github.com/copilot/concepts/agents/about-copilot-cli>

---

## Option B — OpenCode

Open-source coding agent that supports many model providers, including several with **free tiers** (OpenRouter free models, NVIDIA build, Z.AI, etc.). Use this if you don't have a Copilot subscription. Free-tier models typically do not require authentication, so you can skip the login step.

### Install

The official one-line installer:

```bash
curl -fsSL https://opencode.ai/install | bash
```
After installation, the `opencode` command should be available in your terminal. If not, you may need to restart your terminal.

### Run `opencode`

```bash
opencode
```

Inside `opencode`, run the command `/models` to see the available models and select one with a free tier (e.g. `Big Pickle`).


Docs: <https://opencode.ai/docs/>

---

> **Note:** When using `copilot` or `opencode`, make sure you start the agent in the project directory, for example `~/github/r_intro/Proteome_analysis_project`. This way the agent can read and edit files in the current directory. If you start the agent in a different directory, exit it with the `/exit` command, `cd` into the project directory, and start the agent again.


# VS Code (optional)

VS Code is a general-purpose editor and works well alongside either coding agent (you can run the agent in VS Code's integrated terminal). Go to <https://code.visualstudio.com/download> and download the latest version for your OS (`.deb` for Ubuntu).

```bash
cd ~/Downloads
su
gdebi DOWNLOAD_LINK.deb   # replace with the actual filename
```

Recommended extensions for this project:

- **R** (REditorSupport.r) — R language support, linting, REPL integration.
- **Quarto** (quarto.quarto) — preview and render `.qmd` reports.


---

# Windows (10 / 11)

Two practical paths for Windows users:

- **Path A — WSL (Ubuntu on Windows)** — recommended.
- **Path B — Native Windows** (PowerShell).

Pick **Path A** unless you have a reason to avoid WSL: it gives you a real Ubuntu shell inside Windows, so the **Ubuntu** instructions above apply verbatim.

---

## Path A — WSL (Ubuntu on Windows)

Open **PowerShell as Administrator** (right-click the Start menu → "Terminal (Admin)") and run:

```powershell
wsl --install
```

This installs WSL 2 and an Ubuntu distribution. After it finishes:

1. Restart Windows.
2. Open "Ubuntu" from the Start menu. You will be prompted to create a Linux username and password the first time.
3. From that point on, follow the **Ubuntu** instructions at the top of this document **verbatim**, inside the Ubuntu terminal — coding agent, VS Code remote integration, everything.

To use VS Code on top of WSL:

1. Install VS Code natively on Windows (see *VS Code* under Path B below).
2. In VS Code, install the **WSL** extension (`ms-vscode-remote.remote-wsl`).
3. From the Ubuntu terminal, `cd` into your project directory and run `code .` — VS Code will open and connect to the WSL filesystem.

Microsoft's WSL docs: <https://learn.microsoft.com/windows/wsl/install>

---

## Path B — Native Windows

Use **PowerShell** (right-click the Start menu → "Terminal" or "Windows PowerShell"). The commands below assume PowerShell, not the legacy `cmd.exe`.

### Coding agent

Pick **one** of the two options below.

#### Option A — GitHub Copilot CLI

Requires an active **Copilot subscription** (Pro, Pro+, Business, or Enterprise — student Pro accounts also work) and **PowerShell 7+** (the default `Windows PowerShell 5.1` is too old). Install both via `winget`, which is built into Windows 10 / 11:

```powershell
winget install Microsoft.PowerShell   # PowerShell 7+ (needed by Copilot CLI)
winget install GitHub.Copilot
```

Open a **new** PowerShell 7 window (the Start menu entry is "PowerShell 7" or just type `pwsh`) and run:

```powershell
copilot
```

Type `/login` inside the TUI to authenticate via the browser. The slash commands (`/login`, `/model`, `/help`, `/exit`) work the same as on Linux.

Docs: <https://docs.github.com/copilot/concepts/agents/about-copilot-cli>

#### Option B — OpenCode

OpenCode is not on `winget` (yet). The simplest no-admin installer on Windows is **Scoop**. One-time setup of Scoop (skip if you already have it):

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
irm get.scoop.sh | iex
```

Then install and run OpenCode:

```powershell
scoop install opencode
opencode
```

Inside `opencode`, run `/models` to see the available models and select one with a free tier (e.g. `Big Pickle`). Free-tier models typically do not require authentication.

Docs: <https://opencode.ai/docs/>

> **Note:** As on Linux, always start `copilot` or `opencode` in the project directory. In PowerShell:
> ```powershell
> cd C:\Users\<you>\github\r_intro\Proteome_analysis_project
> copilot   # or opencode
> ```

### VS Code (optional)

Download the installer from <https://code.visualstudio.com/download> (the `.exe` for Windows) and run it. Alternatively, with the Windows Package Manager (built into Windows 10 / 11):

```powershell
winget install Microsoft.VisualStudioCode
```

Recommended extensions (same as on Linux):

- **R** (REditorSupport.r) — R language support, linting, REPL integration.
- **Quarto** (quarto.quarto) — preview and render `.qmd` reports.

