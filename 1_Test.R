# Install necessary package if not installed
if (!require(treemap)) install.packages("treemap", dependencies = TRUE)
library(treemap)

# Example data: Replace with your actual vaccine terms and study numbers
vaccine_data <- data.frame(
  Vaccine = c("Vaccine A", "Vaccine B", "Vaccine C", "Vaccine D", "Vaccine E",
              "Vaccine F", "Vaccine G", "Vaccine H", "Vaccine I", "Vaccine J",
              "Vaccine K", "Vaccine L", "Vaccine M", "Vaccine N", "Vaccine O",
              "Vaccine P", "Vaccine Q", "Vaccine R", "Vaccine S", "Vaccine T"),
  Studies = c(500, 300, 250, 200, 150, 400, 350, 320, 280, 260,
              450, 420, 410, 390, 370, 360, 340, 330, 310, 290),
  Subgroup = c("Group 1", "Group 1", "Group 1", "Group 1", "Group 1",
               "Group 2", "Group 2", "Group 2", "Group 2", "Group 2",
               "Group 3", "Group 3", "Group 3", "Group 3", "Group 3",
               "Group 4", "Group 4", "Group 4", "Group 4", "Group 4")
)

# Create treemap
treemap(vaccine_data,
        index = c("Subgroup", "Vaccine"),  # Grouping first by Subgroup, then Vaccine
        vSize = "Studies",  # Variable to determine rectangle size
        title = "Vaccine Studies Treemap",
        palette = "Set3",  # Color palette
        border.col = "white")  # White borders for readability
