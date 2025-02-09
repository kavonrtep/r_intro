---
title: "RStudio Introduction Outline"
author: "Instructor Notes"
format: html
---

# RStudio Introduction Outline

This outline covers the essential topics for introducing RStudio. It is designed to ensure that all key points are discussed during the session.

## 1. Introduction
    - Install RStudio.
    - Navigate the RStudio interface.
    - Write and execute R code.
    - Create and export plots.
    - Create interactive documents using Quarto.
    - Customize the RStudio appearance.

## 2. What is RStudio?
- **Definition**: An Integrated Development Environment (IDE) for R.
- **Main Components/Windows:**
  - **Code Editor**: Where you write R scripts.
  - **Console**: Executes code and displays output.
  - **Environment/Workspace**: Shows variables and work history.
  - **Plots/Files/Packages/Help**: Manages plots, file navigation, package activation, and help documentation.
- **Comparison**: Unlike a basic setup (separate code editor and console), RStudio brings all these elements together in one interface.

## 3. Installing RStudio
- **Steps:**
  - Visit the RStudio website.
  - Download the appropriate installer for your platform.
  - Follow the installation instructions.
- **Prerequisite**: R must already be installed on your computer.

## 4. Navigating the RStudio Interface
- **Starting an R Script:**
  - Go to **File > New File > R Script** to open the code editor.
- **Overview of the Interface:**
  - Identify the four main panels: Code Editor, Console, Environment/Workspace, and the combined Plots/Files/Packages/Help window.

## 5. Customizing Appearance
- **Changing Themes:**
  - Access **RStudio Preferences > Appearance**.
  - Select an editor theme (various background, syntax highlighting, and font color options).
  - Adjust font type and size for readability.
- **Customizing Pane Layout:**
  - Rearrange panels (e.g., move the Console or Files pane) to suit your workflow.

## 6. Writing and Running R Code
- **Creating Variables:**
  - Demonstrate assigning a vector to variables (e.g., `x` and `y`).
  - Introduce the useful keyboard shortcut for the assignment operator (Option/Alt + `-`).
- **Executing Code:**
  - Run a single line with the **Run** button or using **Command/Ctrl + Enter**.
  - Show how multiple lines can be run by highlighting them before executing.
- **Observing Results:**
  - Point out how variables appear in the Global Environment.
  - Explain how the Console displays the output.

## 7. Plotting in RStudio
- **Creating a Basic Plot:**
  - Example: `plot(x, y)`
- **Modifying the Plot:**
  - Change parameters (e.g., `type="b"` for both points and lines).
- **Plot Navigation and Export:**
  - Use the arrows to review multiple plots.
  - Demonstrate the zoom feature.
  - Explain how to export plots as images or PDFs.

## 8. Working with Files and Directories
- **Managing Files:**
  - Use the Files tab to navigate directories.
- **Setting the Working Directory:**
  - Demonstrate setting the working directory via **More > Set Working Directory** or by using the command (e.g., `setwd("path/to/directory")`).
- **Reading Data:**
  - Example: `read.csv("mydata.csv")`
  - Show how to assign the read data to a variable and view it in the Environment pane.

## 9. Managing Packages
- **Installing Packages:**
  - Use `install.packages("package_name")`.
  - Highlight RStudio’s autocompletion features when typing function names.
- **Accessing Documentation:**
  - Use `?function_name` and `??partial_name` to search for help.
- **Activating Packages:**
  - Demonstrate using the Packages tab (checking the box) or via `library(package_name)`.

## 10. Plotting with Package Data
- **Example with an External Package:**
  - Remove and then install a package (e.g., "nycflights13") to show how to install packages.
  - Activate the package and explore a dataset (e.g., the `flights` data).
  - Create a histogram (or other plot) from a column in the dataset.

## 11. Creating Quarto Documents
*Note: This section updates the original “R Notebook” instructions.*
- **What is a Quarto Document?**
  - A modern document format similar to R Markdown that allows interactive code chunks.
- **Creating a New Quarto Document:**
  - Go to **File > New File > Quarto Document**.
  - Note: The template will be pre-populated with some basic information.
- **Working with Code Chunks:**
  - Insert new code chunks using the toolbar button or the shortcut (Command/Option + I).
  - Explain how to execute code chunks individually (e.g., Command/Shift + Enter).
  - Show how to preview the document so that both the code and its output are visible.
  - Demonstrate options to hide or show code chunks.
- **Benefits:**
  - Combines code, output, and narrative in one document.
  - Useful for teaching, reporting, and reproducible analysis.

## 12. Summary and Additional Resources
- **Recap:**
  - Installation and basic navigation of RStudio.
  - Writing and running R code.
  - Creating plots and managing files.
  - Installing and activating packages.
  - Creating interactive documents using Quarto.
- **Encouragement:**
  - Explore additional features in RStudio.
  - Review the RStudio documentation for advanced tips and tricks.
- **Discussion Points:**
  - Ask for questions.
  - Invite students to experiment with customizations and shortcuts.

# End of Outline
