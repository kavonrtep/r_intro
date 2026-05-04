# Proteome analysis project instructions

This directory is an educational example project for using a coding agent in a small, practical data analysis workflow. Keep solutions simple, readable, and functional. Prefer straightforward scripts and explicit steps over abstraction-heavy designs.

## Project goals

- Work with proteome-related data in a clear, reproducible way.
- Store intermediate and final tabular outputs as CSV files.
- Produce simple, interpretable plots and reports.
- Favor teaching value and maintainability over clever code.

## Expected project structure

Use and preserve this layout:

- `data/` - input data, downloaded resources, and small processed helper files
- `results/` - analysis outputs, mainly CSV tables
- `reports/` - rendered reports, figures, and presentation-ready outputs
- `scripts/` - R scripts for data import, cleaning, analysis, and plotting

If a directory is missing and is needed by the workflow, create it with that exact name.

## UniProt download guidance

For UniProt-specific download details, use the local project document:

- `./uniprot_howto.md`

Treat that file as the source of truth for how to query UniProt, resolve proteomes, and download FASTA files.

## R coding style

Write R code in a teaching-friendly style:

- Prefer the **tidyverse** for data import, wrangling, and pipelines
- Prefer **ggplot2** for plots
- Prefer **Quarto** for reporting
- Keep code short and easy to follow
- Use clear variable names
- Split work into small functions only when it improves readability
- Avoid unnecessary metaprogramming, deep nesting, and overly generic helper layers
- Do not introduce complexity just to make the code look advanced

## Implementation guidance for agents

When adding or editing analysis code:

- Assume scripts are run from the project root
- Use relative paths such as `data/...`, `results/...`, and `reports/...`
- Save tabular outputs as CSV unless there is a strong reason to use another format
- Make plots publication-friendly but simple
- Prefer one script per logical task rather than large monolithic workflows
- Add brief comments only where they help a learner understand the step

## Reporting guidance

When creating reports:

- Use Quarto documents when a narrative report is needed
- Keep text concise and focused on biological or analytical interpretation
- Show code and outputs in a way that is easy for students to inspect
- Save rendered outputs under `reports/`

## What to avoid

- Overly complex object-oriented designs
- Large frameworks or unnecessary package dependencies
- Hidden side effects
- Clever but hard-to-read one-liners
- Writing code that is optimized for scale instead of clarity

# Ignores

Ignore these files and directories: notes