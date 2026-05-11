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


