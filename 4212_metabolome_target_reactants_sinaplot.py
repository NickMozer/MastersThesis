# Sinaplot of distribution of the all fold changes of response metabolites that are a known reactant of target enzymes

# Import libraries
import pandas as pd
import os
import plotnine as p9

# Initialize date and output path
output_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(output_dir)

# Load fc, rsd tables
fc_table = pd.read_excel("./219_metabolome_fc_table_low_rsds.xlsx", "Sheet1")

# Reshape the data for boxplot
fc_values = fc_table.iloc[:, 1:].values.flatten()

# Create a DataFrame from the flattened fold change values
df_fc = pd.DataFrame({'fold_change': fc_values, 'category': 'Fold Change'})

# Configuration for plot size, dot colors, Y-axis ticks, and grid
figure_width = 5 # inches
figure_height = 10 # inches
p9.options.figure_size = (figure_width, figure_height)

dot_facecolor = "orange"
dot_edgecolor = "black"
dot_stroke_thickness = 0.5

# Define custom y-axis breaks
y_tick_breaks = [-6, -4, -2, 0, 2, 4, 6, 8, 10]

# Create the sinaplot
sinaplot = (
    p9.ggplot(df_fc, p9.aes(x='category', y='fold_change'))
    + p9.geom_sina(
        alpha=0.7,
        size=3,
        fill=dot_facecolor,
        color=dot_edgecolor,
        stroke=dot_stroke_thickness
    )
    + p9.stat_summary(fun_data="mean_sdl", geom='pointrange', color='gray', size=0.8)
    + p9.labels.ylab("Log2 Fold Change")
    + p9.labels.xlab("Metabolite targets")
    + p9.labels.ggtitle("")
    # Y-axis tick adjustment
    + p9.scale_y_continuous(breaks=y_tick_breaks)
    # Remove background grid
    + p9.theme_minimal()
    + p9.theme(
        panel_grid_major=p9.element_blank(), # Removes major grid lines
        panel_grid_minor=p9.element_blank(), # Removes minor grid lines
        # Set the background of the plotting panel to white
        panel_background=p9.element_rect(fill="white"),
        # Set the overall plot background to white
        plot_background=p9.element_rect(fill="white"),
        # Hide the x-axis tick labels (the category name)
        axis_line_y=p9.element_line(color="black", size=0.5),
        axis_line_x=p9.element_blank(),
        axis_text_x=p9.element_blank(),
        axis_text_y=p9.element_text(size=12, margin={"r": 10, "units": "pt"}),
        axis_ticks_major_y=p9.element_line(color="black", size=0.5)
    )
)
# Save the plot
sinaplot.save("sinaplot_all_fcs.png", dpi=300)

# Display the plot
print(sinaplot)