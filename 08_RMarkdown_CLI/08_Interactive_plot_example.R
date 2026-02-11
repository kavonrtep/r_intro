# Install and load necessary libraries
# Run these commands only once
# install.packages("ggplot2")
# install.packages("plotly")
# install.packages("htmlwidgets")
# install.packages("palmerpenguins")

library(ggplot2)
library(plotly)
library(htmlwidgets)
library(palmerpenguins)

# Create scatter plot with meaningful text aesthetics (for tooltips)
gg <- ggplot(data = penguins, aes(x = bill_length_mm, y = bill_depth_mm,
                                  color = species, shape = species,
                                  text = paste("Species: ", species,
                                               "<br>Island: ", island,
                                               "<br>Bill Length: ", bill_length_mm, "mm",
                                               "<br>Bill Depth: ", bill_depth_mm, "mm"))) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "Palmer Penguins Bill Dimensions",
       x = "Bill Length (mm)",
       y = "Bill Depth (mm)") +
  theme_minimal()

# non interactive plot
gg


# Convert ggplot to interactive plot using ggplotly, specifying tooltip
gg_interactive <- ggplotly(gg, tooltip = "text")

# View interactive plot
gg_interactive

# Export interactive plot as HTML
saveWidget(gg_interactive, file = "penguin_scatter_interactive.html")