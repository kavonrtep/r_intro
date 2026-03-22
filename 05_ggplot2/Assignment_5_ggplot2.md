# Assignment 5: Data Visualisation with ggplot2 — Outbreak Epidemiology

## Background

Visualising epidemiological data is one of the most impactful applications of data science in public health. In this assignment you will build a series of publication-quality figures using three real outbreak datasets originally distributed through the [`outbreaks`](https://CRAN.R-project.org/package=outbreaks) R package:

| File | Dataset | Description |
|---|---|---|
| `ebola_kikwit_1995.csv` | Ebola virus disease, Kikwit, DRC | Daily epidemic curve: onset and death counts |
| `fluH7N9_china_2013.csv` | Avian influenza A(H7N9), China | Individual case linelist: 136 patients |
| `covid19_cases_wide.csv` | COVID-19, WHO reports 2020 | Cumulative cases per country (wide format) |

All files are in `05_ggplot2/data/`. The assignment practices data import, `scale_x_date()`, epidemic curve plotting, multi-line time series, faceting, violin + jitter plots, log scales, and `ggsave()`.

---

## Setup

```r
library(ggplot2)
library(tidyr)    # for pivot_longer()
```

---

## Task 1: Load and Inspect All Three Datasets

Load each CSV file and inspect its structure.

- Use `read.csv()` with `stringsAsFactors = FALSE`.
- For date columns, convert character strings to `Date` objects using `as.Date()`.
- Print the dimensions and first few rows of each dataset.
- How many rows and columns does each dataset have? What does one row represent in each?

<details>
<summary>💡 Hint</summary>

Use `read.csv("05_ggplot2/data/ebola_kikwit_1995.csv", stringsAsFactors = FALSE)` and repeat for the other two files. After loading, convert date columns with `as.Date()` — for example `ebola$date <- as.Date(ebola$date)`. Inspect with `dim()` and `head()`.

The flu dataset has several date columns (`date_of_onset`, `date_of_hospitalisation`, `date_of_outcome`) — convert all of them.

</details>

---

## Task 2: Ebola Epidemic Curve — Bar Chart with `scale_x_date()`

Plot the daily number of new **onset** cases as a bar chart over time.

- Map `date` to x and `onset` to y.
- Use `geom_col()` (or `geom_bar(stat = "identity")`).
- Format the x-axis to show months using `scale_x_date(date_labels = "%b %Y", date_breaks = "1 month")`.
- Rotate x-axis labels 45° with `theme(axis.text.x = element_text(angle = 45, hjust = 1))`.
- Add a meaningful title and axis labels.

<details>
<summary>💡 Hint</summary>

`scale_x_date()` works like `scale_x_continuous()` but understands `Date` objects. The `date_labels` argument uses `strftime` format codes: `%b` = abbreviated month name, `%Y` = 4-digit year. The `date_breaks` argument accepts strings like `"1 month"` or `"2 weeks"`.

</details>

---

## Task 3: Ebola — Overlay Onset and Deaths

Extend the previous plot by overlaying both **onset** and **death** counts on the same axes.

- Reshape the `ebola` data from wide to long format using `pivot_longer()`.
- Map `type` (onset/death) to `fill`.
- Use `position = "dodge"` to place bars side by side, or keep them overlapping with different `alpha`.
- Use `scale_fill_manual()` to assign meaningful colours (e.g., red for onset, black for death).

<details>
<summary>💡 Hint</summary>

`pivot_longer()` (from the `tidyr` package) converts a **wide** data frame into a **long** one.
Right now `ebola` has two separate columns `onset` and `death` — one row per day.
After reshaping, each day becomes **two rows**: one for onset, one for death — and a new
column records which type it is. This is the format ggplot2 needs to map the series to `fill`.

```
# BEFORE (wide):
#   date       onset  death
#   1995-04-01   3      1

# AFTER (long):
#   date       type    count
#   1995-04-01 onset   3
#   1995-04-01 death   1
```

- `cols` — which columns to collapse (here: `onset` and `death`)
- `names_to` — name of the new column that will hold the old column names (`"type"`)
- `values_to` — name of the new column that will hold the values (`"count"`)

```r
ebola_long <- pivot_longer(ebola,
                           cols      = c(onset, death),
                           names_to  = "type",
                           values_to = "count")
```

Then build your `ggplot()` using `ebola_long`, mapping `count` to y and `type` to `fill`.

</details>

---

## Task 4: Flu H7N9 — Cases by Province (Reordered Bar Chart)

Count the number of H7N9 cases per Chinese province and visualise as a horizontal bar chart, ordered from most to fewest cases.

- Use `table()` to tally cases per `province`, then convert to a data frame with `as.data.frame()`.
- Reorder provinces by count using `reorder()`.
- Use `geom_col()` + `coord_flip()` for a horizontal layout.
- Colour bars by province or use a single colour; add value labels with `geom_text()` if you like.

<details>
<summary>💡 Hint</summary>

`as.data.frame(table(flu$province))` gives a two-column data frame — rename the columns with `colnames()` so they have meaningful names before passing to `ggplot()`.

Inside `aes()`, wrap the province variable in `reorder(province, count)` to sort bars by count. After `coord_flip()`, the longest bar will appear at the top.

</details>

---

## Task 5: Flu H7N9 — Outcome by Gender (Violin + Jitter)

Compare the time from onset to hospitalisation between patients who died and those who recovered, split by gender.

- Calculate `days_to_hospital` as the difference in days between `date_of_hospitalisation` and `date_of_onset` — remember to handle `NA` values.
- Map `outcome` to x and `days_to_hospital` to y.
- Use `geom_violin(fill = "lightblue")` as the base layer.
- Add `geom_jitter(width = 0.15, alpha = 0.6)` on top to show individual points.
- Facet by `gender` using `facet_wrap(~ gender)`.

<details>
<summary>💡 Hint</summary>

Subtracting two `Date` columns gives a `difftime` object — wrap it in `as.numeric()` to get plain days.

Rows where either `days_to_hospital` or `outcome` is `NA` will cause problems — filter them out with `flu[!is.na(flu$days_to_hospital) & !is.na(flu$outcome), ]` before plotting.

Add `theme(legend.position = "none")` to suppress the redundant fill legend when outcome is already on the x-axis.

</details>

---

## Task 6: COVID-19 — Reshape Wide to Long and Plot Multi-Country Lines

The COVID-19 dataset is in **wide format** (one column per country). Reshape it to long format and plot cumulative cases over time for all countries.

- Use `pivot_longer()` to convert country columns to rows.
- Map `date` to x, `cases` to y, and `country` to colour.
- Use `geom_line(linewidth = 1)`.
- Add `scale_x_date(date_labels = "%b", date_breaks = "2 weeks")` and rotate labels.
- Use `scale_color_brewer(palette = "Set1")` or `scale_color_manual()` for distinct colours.

<details>
<summary>💡 Hint</summary>

The COVID-19 dataset is in **wide format**: one column per country, each holding that country's
cumulative case count. ggplot2 cannot directly draw nine separate lines from this layout — it
needs the data in **long format**: one row per date × country combination, with a single `cases`
column and a `country` column identifying which country each row belongs to.

```
# BEFORE (wide) — 1 row per date, 9 country columns:
#   date       China  South_Korea  Japan  ...
#   2020-01-20   278            1      1

# AFTER (long) — 9 rows per date, 1 value column:
#   date       country      cases
#   2020-01-20 China          278
#   2020-01-20 South_Korea      1
#   2020-01-20 Japan            1
#   ...
```

`cols = -date` means "pivot all columns **except** `date`". The `-` sign excludes that column
from the reshape; it stays as the grouping identifier.

```r
covid_long <- pivot_longer(covid_wide,
                           cols      = -date,   # pivot everything except the date column
                           names_to  = "country",
                           values_to = "cases")
```

Then build your `ggplot()` using `covid_long`.

</details>

---

## Task 7: COVID-19 — Log10 Scale and Faceting

Plotting on a linear scale obscures countries with fewer cases. Apply a log10 transformation and then facet by country.

- Add `scale_y_log10()` to the previous plot. What changes?
- Create a second version using `facet_wrap(~ country, scales = "free_y")` **without** the log scale — does free scaling reveal anything that the log scale hides?
- Add `geom_vline(xintercept = as.Date("2020-03-11"), linetype = "dashed")` to mark the WHO pandemic declaration date.

<details>
<summary>💡 Hint</summary>

`scale_y_log10()` drops in as a replacement for the default linear y scale — add it as an extra layer. Note that `log10(0)` is undefined, so if any country has 0 cases you may see warnings; adding `+ 1` before the scale (`y = cases + 1`) avoids this.

For `geom_vline()`, pass `xintercept = as.Date("2020-03-11")` — the date must be a `Date` object, not a plain string, because the x-axis is a date scale.

To exclude China from the faceted version, subset `covid_long` with `[covid_long$country != "China", ]`.

</details>

---

## Task 8: Export Your Best Figure

Choose one of the plots you created and export it as both a PNG (for web/slides) and a PDF (for print/publication).

- Use `ggsave()` with explicit `width`, `height`, and `dpi` arguments.
- Save to `work_dir/` (create it with `dir.create("work_dir", showWarnings = FALSE)` if needed).
- Try `width = 10, height = 6, dpi = 300` for the PNG and `width = 18, height = 12, units = "cm"` for the PDF.

<details>
<summary>💡 Hint</summary>

`ggsave()` saves the last plot printed by default, or you can pass a plot object explicitly with `ggsave("file.png", plot = p, ...)`. PDFs are vector graphics and do not need a `dpi` argument. Use `units = "cm"` when you want to specify physical dimensions precisely.

</details>

---

## Bonus Task: Annotated Epidemic Curve

Return to the Ebola epidemic curve from Task 2 and add annotations marking key events:

- Draw a vertical line (`geom_vline()`) on `1995-05-06` — the date the outbreak was officially recognised.
- Add a text annotation (`annotate("text", ...)`) labelling that line.
- Add a shaded rectangle (`annotate("rect", ...)`) covering the peak weeks.
- Use `geom_smooth(method = "loess", se = TRUE)` to overlay a smoothed trend on the onset bars.

<details>
<summary>💡 Hint</summary>

`annotate()` adds fixed elements not linked to the data. Key arguments:
- `"rect"` requires `xmin`, `xmax`, `ymin`, `ymax` — use `as.Date()` for the x values and `Inf` for `ymax` to extend to the top of the panel.
- `"text"` requires `x`, `y`, and `label`. Use `hjust = 0` to left-align the text from the given x position.

Layer order matters: put `annotate("rect", ...)` before `geom_col()` so the shading appears behind the bars.

</details>
