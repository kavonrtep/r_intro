# Setup

This page covers the tools needed for the project on **Ubuntu** (tested on 24.04 LTS):

1. A coding agent — either **GitHub Copilot CLI** (requires a Copilot subscription) or **OpenCode** (works with free model providers).
2. **VS Code** as a general-purpose editor.
3. **R packages** for tree handling, plotting, and reporting.

If you already have a working coding agent and editor, skip directly to the *R packages* section at the bottom.

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
On first launch you will be prompted for authentication. Type:

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

Open-source coding agent that supports many model providers, including several with **free tiers** (OpenRouter free models, NVIDIA build, Z.AI, etc.). Use this if you don't have a Copilot subscription.

### Install

The official one-line installer:

```bash
curl -fsSL https://opencode.ai/install | bash
```

### Run and configure a provider

```bash
opencode
```

```bash
opencode auth list
```

Docs: <https://opencode.ai/docs/>

---

# VS Code

VS Code is a general-purpose editor and works well alongside either coding agent (you can run the agent in VS Code's integrated terminal). Go to <https://code.visualstudio.com/download> and download the latest version for your OS (.deb for Ubuntu).

```bash
cd ~/Downloads
sudo gdebi DOWNLOAD_LINK.deb
```

Recommended extensions for this project:

- **R** (REditorSupport.r) — R language support, linting, REPL integration.
- **Quarto** (quarto.quarto) — preview and render `.qmd` reports.

---

# R packages

Install once before starting the analysis. CRAN packages first:

```r
install.packages(c(
  "ape", "treeio", "ggplot2", "dplyr", "tidyr", "readr",
  "tibble", "lubridate", "patchwork",
  "ggnewscale", "ggimage", "viridisLite",
  "phytools"
))
```

Bioconductor packages (`ggtree`, `ggtreeExtra`, `Biostrings`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("ggtree", "ggtreeExtra", "Biostrings"))
```

Quarto is needed for the final report. On Ubuntu:

```bash
# Latest .deb from https://quarto.org/docs/get-started/
sudo gdebi quarto-*-linux-amd64.deb
```

Verify the toolchain:

```r
library(ape); library(treeio); library(ggtree); library(ggtreeExtra)
library(ggplot2); library(ggnewscale); library(dplyr); library(tidyr)
library(readr); library(lubridate); library(patchwork); library(Biostrings)
```

If all `library()` calls succeed without warnings, the project is ready to run.
