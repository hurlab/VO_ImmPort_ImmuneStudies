################################################################################
#####################Set Current Working Directory ###########################
################################################################################

# Ensure rstudioapi is installed
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}
library(rstudioapi)

# Set the working directory to the script’s location
current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
base_dir <- dirname(current_path)

################################################################################
#####################  Packages  ###############################################
################################################################################

# Install pacman for package management if not already installed
if (!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman")
}
library(pacman)

# Required CRAN packages
cran_packages <- c(
  "ontologyIndex", "ontologySimilarity", "readxl", "dplyr", "tidyr", 
  "data.table", "treemapify", "ggplot2"
)

# Automatically install and load CRAN packages
p_load(char = cran_packages, install = TRUE)

################################################################################
#####################  Load Ontology from OBO File #############################
################################################################################

# Define the OBO file path
obo_file <- file.path(base_dir, "vo.obo")

# Load ontology using `ontologyIndex`
ontology <- get_OBO(obo_file, extract_tags = "everything")

# Ensure root ontology term is present
if (!"VO:0000001" %in% ontology$id) {
  stop("VO:0000001 is not present in the ontology structure. Please check the OBO file.")
}

################################################################################
#####################  Load Study Data from Excel ##############################
################################################################################

# Define the Excel file path
excel_file <- file.path(base_dir, "UniqueStudy_VOIDs.xlsx")

# Load study data
figure_data <- read_excel(excel_file, sheet = "Data for Figure") %>%
  rename(VO_ID = 1, Vaccine = 2, Study_num = 3, Pathogen_Type = 4, Disease = 5)

immport_data <- read_excel(excel_file, sheet = "ImmPortSTUDY_VO_Simple") %>%
  rename(VO_ID = 1, STUDY_ACCESSION = 2, VO_Label = 3)

# Convert VO_ID format from "VO_0000738" to "VO:0000738"
convert_vo_format <- function(vo_id) {
  return(gsub("VO_", "VO:", vo_id))
}

figure_data$VO_ID <- sapply(figure_data$VO_ID, convert_vo_format)
immport_data$VO_ID <- sapply(immport_data$VO_ID, convert_vo_format)

# Merge data
vo_map <- figure_data %>%
  select(VO_ID, Pathogen_Type, Disease) %>%
  filter(complete.cases(.)) %>%
  distinct()

merged_data <- immport_data %>%
  left_join(vo_map, by = "VO_ID") %>%
  filter(complete.cases(.)) %>%
  distinct()

# Print data summary
cat("Original immport data rows:", nrow(immport_data), "\n")
cat("Merged data rows with complete data:", nrow(merged_data), "\n")

# Create study count table
explicit_counts <- merged_data %>%
  group_by(VO_ID, VO_Label, Pathogen_Type, Disease) %>%
  summarise(
    STUDY_ACCESS_IDs = paste(sort(unique(STUDY_ACCESSION)), collapse=", "),
    explicit_study_num = n_distinct(STUDY_ACCESSION),
    .groups = "drop"
  )

################################################################################
#####################  Propagate Study Assignments Up the Tree #################
################################################################################

# Function to get all ancestors of a given ontology term
get_all_ancestors <- function(ontology, term) {
  if (!term %in% ontology$id) return(NULL)
  ancestors <- get_ancestors(ontology, term)
  
  # Filter only those that start with "VO:"
  valid_ancestors <- ancestors[grepl("^VO:", ancestors)]
  
  return(valid_ancestors)
}

# Get all unique ontology terms from the study data
unique_vo_terms <- unique(explicit_counts$VO_ID)

# Create a mapping of ontology terms to all their ancestor terms
ontology_mapping <- lapply(unique_vo_terms, function(term) {
  data.frame(Ontology_Term = term, Ancestor_Term = get_all_ancestors(ontology, term))
}) %>%
  bind_rows()

# Merge study data with all ancestor terms
expanded_study_data <- explicit_counts %>%
  rename(Explicit_Ontology_Term = VO_ID) %>%
  left_join(ontology_mapping, by = c("Explicit_Ontology_Term" = "Ontology_Term"))

################################################################################
#####################  Count Studies for Each Ontology Term ####################
################################################################################

# Count unique study accessions for each ontology term (both explicit and propagated)
study_counts <- expanded_study_data %>%
  group_by(Ancestor_Term) %>%
  summarise(
    Study_Count = n_distinct(STUDY_ACCESS_IDs),
    .groups = "drop"
  )

# Print summary
cat("\nFinal Ontology Study Assignment Summary:\n")
cat("Total unique ontology terms:", nrow(study_counts), "\n")



################################################################################
#####################  Extract Level 2 and Level 3 Terms from Ontology #########
################################################################################

# Get direct children of VO:0000001 (Level 2 terms)
level2_terms <- ontology$children[["VO:0000001"]]

# Ensure valid Level 2 terms exist
if (is.null(level2_terms) || length(level2_terms) == 0) {
  stop("No Level 2 terms found under VO:0000001. Check the ontology structure.")
}

# Get Level 3 terms (children of Level 2 terms)
level3_terms <- unlist(lapply(level2_terms, function(term) {
  ontology$children[[term]]
}))

# Remove any NULL values
level3_terms <- level3_terms[!is.na(level3_terms)]

# Ensure Level 3 terms are found
if (length(level3_terms) == 0) {
  stop("No Level 3 terms found. Check the ontology structure.")
}

# Function to find the closest Level 2 ancestor for a given VO_ID
find_closest_level2 <- function(term) {
  if (!term %in% ontology$id) return(NA)  # Skip if term is not in ontology
  
  # Get all ancestors of this term
  ancestors <- get_ancestors(ontology, term)
  
  # Find the closest ancestor that belongs to Level 2
  closest_level2 <- ancestors[ancestors %in% level2_terms]
  
  if (length(closest_level2) > 0) {
    return(tail(closest_level2, 1))  # Get the closest (last in the list)
  } else {
    return(NA)
  }
}

# Function to find the closest Level 3 ancestor for a given VO_ID
find_closest_level3 <- function(term) {
  if (!term %in% ontology$id) return(NA)  # Skip if term is not in ontology
  
  # Get all ancestors of this term
  ancestors <- get_ancestors(ontology, term)
  
  # Find the closest ancestor that belongs to Level 3
  closest_level3 <- ancestors[ancestors %in% level3_terms]
  
  if (length(closest_level3) > 0) {
    return(tail(closest_level3, 1))  # Get the closest (last in the list)
  } else {
    return(NA)
  }
}

################################################################################
#####################  Create Simplified Study Data (Top 3 Levels) #############
################################################################################

# Apply functions to find closest Level 2 and Level 3 ancestors for each term
simplified_mapping <- data.frame(
  VO_ID = unique_vo_terms,
  Level2_Ancestor = sapply(unique_vo_terms, find_closest_level2),
  Level3_Ancestor = sapply(unique_vo_terms, find_closest_level3),
  stringsAsFactors = FALSE
)

# Merge study data with simplified mapping
simplified_study_data <- explicit_counts %>%
  left_join(simplified_mapping, by = "VO_ID") %>%
  filter(!is.na(Level2_Ancestor) | !is.na(Level3_Ancestor))  # Keep terms with Level 2 or 3 mapping

# Aggregate study counts at the Level 2 ancestor level
simplified_study_counts_level2 <- simplified_study_data %>%
  group_by(Level2_Ancestor) %>%
  summarise(
    Study_Count = sum(explicit_study_num, na.rm = TRUE),
    .groups = "drop"
  )

# Aggregate study counts at the Level 3 ancestor level
simplified_study_counts_level3 <- simplified_study_data %>%
  group_by(Level3_Ancestor) %>%
  summarise(
    Study_Count = sum(explicit_study_num, na.rm = TRUE),
    .groups = "drop"
  )

# Print summaries
cat("\nSimplified Study Assignment Summary (Top 2 Levels):\n")
cat("Total Level 2 terms:", nrow(simplified_study_counts_level2), "\n")

cat("\nSimplified Study Assignment Summary (Top 3 Levels):\n")
cat("Total Level 3 terms:", nrow(simplified_study_counts_level3), "\n")







################################################################################
#####################  Generate Treemap for Each Pathogen Type #################
################################################################################

# Prepare the dataset for Level 2 treemap plotting
treemap_data_level2 <- simplified_study_data %>%
  select(VO_ID, Level2_Ancestor, Pathogen_Type, explicit_study_num) %>%
  rename(Study_num = explicit_study_num) %>%
  filter(!is.na(Level2_Ancestor))  # Keep only terms mapped to Level 2

# Prepare the dataset for Level 3 treemap plotting
treemap_data_level3 <- simplified_study_data %>%
  select(VO_ID, Level3_Ancestor, Pathogen_Type, explicit_study_num) %>%
  rename(Study_num = explicit_study_num) %>%
  filter(!is.na(Level3_Ancestor))  # Keep only terms mapped to Level 3

# Function to generate nested treemap for each pathogen type
plot_nested_treemap <- function(treemap_data, level) {
  pathogen_types <- unique(treemap_data$Pathogen_Type)
  
  for (pathogen in pathogen_types) {
    data_subset <- treemap_data %>% filter(Pathogen_Type == pathogen)
    
    if (nrow(data_subset) > 0) {
      p <- ggplot(data_subset, aes(area = Study_num, fill = !!sym(level), subgroup = !!sym(level), label = VO_ID)) +
        geom_treemap() +
        geom_treemap_subgroup_border() +  # Borders for subgroup
        geom_treemap_subgroup_text(place = "centre", grow = TRUE, alpha = 0.5, color = "black") +  # Labels for subgroup
        geom_treemap_text(color = "black", place = "centre", size = 5) +  # Labels for original terms
        theme(legend.position = "none") +
        labs(title = paste("Vaccine Study Counts -", pathogen, "(", level, ")"))
      
      print(p)  # Display the plot
    }
  }
}

# Generate hierarchical treemaps for Level 2 and Level 3
cat("\nGenerating Level 2 Treemaps...\n")
plot_nested_treemap(treemap_data_level2, "Level2_Ancestor")

cat("\nGenerating Level 3 Treemaps...\n")
plot_nested_treemap(treemap_data_level3, "Level3_Ancestor")


