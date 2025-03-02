## Follow-Up Assignment: Exploring Chick Growth Data 

### Overview 
In this assignment, you will work with the **ChickWeight**  dataset. This dataset records the weight of chicks measured over time while they are fed different diets. You will practice loading data, subsetting data frames, calculating simple statistics, and creating basic plots.

---


### Task 1: Data Exploration 
 
1. **Load and Inspect the Data:**  
  - **Step:**  Ensure the **ChickWeight**  dataset is available in R (it is built-in).
   - **Instructions:**  
    - Use `head(ChickWeight)` to display the first few rows.
    - Use `str(ChickWeight)` and `summary(ChickWeight)` to view the structure and summary statistics of the dataset.
   - **Hint:**  You can type `?ChickWeight` in the console to read more about the dataset.

---

### Task 2: Subsetting Data for a Specific Diet and Time 
 
1. **Subset for Diet 1:**  
  - **Step:**  Create a new data frame that includes only the observations where the diet is 1.
  - **Instructions:**  
    - Use logical indexing or the `subset()` function to filter the data.
    - Save the result in a variable named `diet1_data`.
<details>
<summary>💡 Hint</summary>
  You might write something like

```R
diet1_data <- ChickWeight[ChickWeight$Diet == 1, ]
```
or
```R
diet1_data <- subset(ChickWeight, Diet == 1)
```

</details>

2. **Subset for the Final Time Point:**  
  - **Step:**  From `diet1_data`, extract only the rows where the time is 21 (the last measurement day).
 
  - **Instructions:**  
    - Save this subset as `diet1_time21`.
    - Use `head(diet1_time21)` to check your new data frame.

<details>
<summary>💡 Hint</summary> 
Use a similar approach as before, for example:

```R
diet1_time21 <- diet1_data[diet1_data$Time == 21, ]
```
</details>


---


### Task 3: Calculating Average Weight 
 
1. **Compute the Mean Weight:**  
  - **Step:**  For the subset `diet1_time21`, calculate the average (mean) chick weight.
 
  - **Instructions:**  
    - Use the `mean()` function on the weight column.
    - Print the result.
 
<details>
<summary>💡 Hint</summary> 

```R
average_weight <- mean(diet1_time21$weight)
print(average_weight)
```
This will give you the average weight of chicks on Diet 1 at time 21.
</details>

---


### Task 4: Basic Plotting 
 
1. **Scatter Plot of Chick Growth Over Time:**  
  - **Step:**  Create a scatter plot to visualize how chick weight changes over time using the entire dataset.
 
  - **Instructions:**  
    - Plot `Time` (x-axis) versus `weight` (y-axis) using the `plot()` function.
    - Label the axes as "Time (days)" and "Chick Weight".
    - Add a title "Chick Growth Over Time".
 
<details>
<summary>💡 Hint</summary> 

```R
plot(ChickWeight$Time, ChickWeight$weight,
     xlab = "Time (days)",
     ylab = "Chick Weight",
     main = "Chick Growth Over Time")
```
</details>

---


### Task 5: Bonus Challenge (Optional) 
 
1. **Separate Plots for Each Diet:**  
  - **Step:**  Create individual scatter plots for each diet type (Diet 1, Diet 2, etc.).
   - **Instructions:**  
    - Use the `subset()` function to create a separate data frame for each diet.
    - Create a scatter plot for each subset.
 
<details>
<summary>💡 Hint</summary>
For example, for Diet 2:

```R
diet2_data <- subset(ChickWeight, Diet == 2)
plot(diet2_data$Time, diet2_data$weight,
     xlab = "Time (days)",
     ylab = "Chick Weight",
     main = "Chick Growth Over Time for Diet 2")
```
</details>

---


## Submission 
 
- **Script File:**  Save your completed assignment in an R script file (e.g., `chick_weight_assignment.R`). 
- **Error-Free:**  Ensure your script runs without errors and that each task produces the expected output.
- Send me a emails with the script attached 
