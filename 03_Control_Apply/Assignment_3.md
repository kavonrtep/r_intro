## Assignment: CpG Island Detection — Beyond the For Loop

### Overview

In the in-class exercise, you detected CpG islands by scanning a DNA sequence with
a `for` loop. That approach works, but R offers more elegant and efficient
alternatives. In this assignment, you will solve the **same CpG island detection
problem** using three different approaches — each introducing a key R concept:

- Replacing a for loop with `sapply()`
- Using **vectorized** functions for maximum performance
- Writing a **reusable function** with default arguments
- Comparing execution time with `system.time()`

### Setup

Load the DNA sequence and the stringi library (same as in the exercise):

```R
library(stringi)
dna_sequence <- readRDS("data/dna_sequence.rds")
window_size <- 100
n_windows <- nchar(dna_sequence) - window_size + 1
```

---

### Task 1: Replace the for loop with `sapply()`

**Concept:** The `sapply()` function applies a function to each element of a
vector and returns a simplified result. Instead of managing an index variable and
pre-allocating a result vector, you provide a vector and a function:

```R
# For loop pattern:
result <- numeric(length(x))
for (i in seq_along(x)) {
  result[i] <- some_function(x[i])
}

# Equivalent sapply pattern:
result <- sapply(x, some_function)
```

The function can be an **anonymous function** defined inline:
`sapply(x, function(i) { ... })`.

**Your task:**

1. Create a vector of all starting positions: `positions <- 1:n_windows`
2. Use `sapply()` with an anonymous function that takes a position `i`, extracts
   the window with `substr()`, and returns the CpG count with
   `stri_count_fixed()`
3. Plot the result — it should look identical to the for-loop version

<details>
<summary>💡 Hint</summary>

```R
positions <- 1:n_windows
cpg_sapply <- sapply(positions, function(i) {
  window <- substr(dna_sequence, i, i + window_size - 1)
  stri_count_fixed(window, "CG")
})

plot(cpg_sapply, type = "l",
     xlab = "Position (bp)", ylab = "CpG count (per 100 bp)",
     main = "CpG Profile (sapply)", col = "darkblue")
abline(h = mean(cpg_sapply), col = "red", lty = 2)
```

</details>

---

### Task 2: Fully vectorized approach

**Concept:** Many R functions are **vectorized** — they can operate on entire
vectors element-wise in a single call, without any loop or apply. This is
typically the fastest approach in R because the iteration happens in optimized
C code under the hood.

Key insight — `substr()` accepts vectors for `start` and `stop`:

```R
# This extracts three substrings at once:
substr("ABCDEFGH", c(1, 3, 5), c(3, 5, 7))
# [1] "ABC" "CDE" "EFG"
```

Similarly, `stri_count_fixed()` is vectorized — it counts the pattern in each
element of a character vector.

**Your task:**

1. Create a vector of start positions and a vector of end positions for all
   windows
2. Use a **single** `substr()` call to extract all 9901 windows at once (the
   result is a character vector)
3. Use a **single** `stri_count_fixed()` call to count `"CG"` in every window
4. Plot the result

<details>
<summary>💡 Hint</summary>

```R
starts <- 1:n_windows
ends <- starts + window_size - 1
all_windows <- substr(dna_sequence, starts, ends)
cpg_vectorized <- stri_count_fixed(all_windows, "CG")

plot(cpg_vectorized, type = "l",
     xlab = "Position (bp)", ylab = "CpG count (per 100 bp)",
     main = "CpG Profile (vectorized)", col = "darkblue")
abline(h = mean(cpg_vectorized), col = "red", lty = 2)
```

</details>

---

### Task 3: Write a reusable function

**Concept:** When you solve a problem that could apply to different inputs, wrap
it in a **function**. Good functions have clear parameters with sensible defaults
so they can be called with minimal effort:

```R
my_function <- function(data, param = default_value) {
  # compute something
  return(result)
}
```

**Your task:**

1. Write a function `cpg_profile(sequence, window_size = 100)` that:
   - Takes a DNA string and a window size as arguments
   - Uses the vectorized approach from Task 2 inside
   - Returns a numeric vector of CpG counts per window
2. Test your function on `dna_sequence` with window sizes **50**, **100**, and
   **200**
3. Create a 3-panel plot using `par(mfrow = c(3, 1))` showing the CpG profile
   at each resolution. How does window size affect the smoothness of the profile?
   Reset the layout with `par(mfrow = c(1, 1))` when done.

<details>
<summary>💡 Hint</summary>

```R
cpg_profile <- function(sequence, window_size = 100) {
  n <- nchar(sequence) - window_size + 1
  starts <- 1:n
  ends <- starts + window_size - 1
  windows <- substr(sequence, starts, ends)
  stri_count_fixed(windows, "CG")
}

par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
for (ws in c(50, 100, 200)) {
  profile <- cpg_profile(dna_sequence, ws)
  plot(profile, type = "l",
       main = paste("Window size:", ws, "bp"),
       xlab = "Position (bp)", ylab = "CpG count",
       col = "darkblue")
  abline(h = mean(profile), col = "red", lty = 2)
}
par(mfrow = c(1, 1))
```

</details>

---

### Task 4: Compare performance

**Concept:** R provides `system.time()` to measure how long a block of code takes
to execute. The `elapsed` component reports wall-clock time in seconds:

```R
system.time({
  # code to measure
})
```

**Your task:**

1. Wrap each of your three approaches (for loop from the exercise, `sapply` from
   Task 1, vectorized from Task 2) inside `system.time({ ... })` and note the
   `elapsed` time
2. Which approach is fastest? By roughly how much?
3. Think about why: `sapply` still calls `substr` and `stri_count_fixed` once
   per window (like the for loop), while the vectorized version calls each
   function only once for all windows

<details>
<summary>💡 Hint</summary>

```R
# For loop
time_loop <- system.time({
  cpg_loop <- numeric(n_windows)
  for (i in 1:n_windows) {
    window <- substr(dna_sequence, i, i + window_size - 1)
    cpg_loop[i] <- stri_count_fixed(window, "CG")
  }
})

# sapply
time_sapply <- system.time({
  cpg_sapply <- sapply(1:n_windows, function(i) {
    stri_count_fixed(substr(dna_sequence, i, i + window_size - 1), "CG")
  })
})

# Vectorized
time_vec <- system.time({
  starts <- 1:n_windows
  windows <- substr(dna_sequence, starts, starts + window_size - 1)
  cpg_vec <- stri_count_fixed(windows, "CG")
})

cat("For loop: ", time_loop["elapsed"], "s\n")
cat("sapply:   ", time_sapply["elapsed"], "s\n")
cat("Vectorized:", time_vec["elapsed"], "s\n")
```

</details>

---

## Submission

- **Script File:** Save your completed assignment as an R script file (e.g.,
  `cpg_island_assignment.R`)
- **Error-Free:** Ensure your script runs without errors and produces the
  expected plots
- **Submission:** Send your R script file as an attachment via email
