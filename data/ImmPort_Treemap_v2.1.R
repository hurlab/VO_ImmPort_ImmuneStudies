################################################################################
#####################  Set Current Working Directory ###########################
################################################################################

# Ensure rstudioapi is installed
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}
library(rstudioapi)

# Set the working directory to the script's location
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
  "rdflib", "readxl", "ggplot2", "treemapify", "dplyr", "tidyr", "data.table",
  "igraph"
)

# Automatically install and load CRAN packages
p_load(char = cran_packages, install = TRUE)

################################################################################
#####################  Load and Parse OWL File ##################################
################################################################################

# Define the OWL file path
owl_file <- file.path(base_dir, "vo.owl")

# Load OWL file as RDF
rdf <- rdf_parse(owl_file)

# Function to extract ontology relationships
get_filtered_ontology_relations <- function(rdf) {
  query <- '
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    SELECT ?child ?parent WHERE { ?child rdfs:subClassOf ?parent . }'
  
  triples <- rdf_query(rdf, query)
  
  if (is.null(triples) || nrow(triples) == 0) {
    stop("No subclass relationships found. Ensure the OWL file is structured correctly.")
  }
  
  colnames(triples) <- c("Child", "Parent")
  triples <- triples %>%
    mutate(
      Child = gsub(".*/", "", Child),
      Parent = gsub(".*/", "", Parent)
    ) %>%
    filter(grepl("^VO_", Child) & grepl("^VO_", Parent)) %>%
    distinct(Child, Parent)
  
  return(triples)
}

################################################################################
#####################  Build Hierarchy Network (Find Level 2 Ancestors) ########
################################################################################

# Extract ontology relations
ontology_relations <- get_filtered_ontology_relations(rdf)

# Create directed graph
g <- graph_from_data_frame(ontology_relations, directed = TRUE)

# Ensure VO_0000001 is present
if (!"VO_0000001" %in% V(g)$name) {
  stop("VO_0000001 is not in the ontology. Check the OWL file.")
}

# Get all direct children of VO_0000001 (Level 2 terms)
level2_terms <- neighbors(g, "VO_0000001", mode = "in")$name

# Function to find the closest Level 2 ancestor for a given VO_ID
find_closest_level2 <- function(vo_id) {
  if (!vo_id %in% V(g)$name) return(NA)  # Skip if term is not in ontology
  
  # Get all shortest paths from the given term to Level 2 terms
  paths <- shortest_paths(g, from = vo_id, to = level2_terms, mode = "out")$vpath
  
  # Find the shortest path to any Level 2 term
  valid_paths <- Filter(length, paths)  # Remove empty paths
  if (length(valid_paths) == 0) return(NA)  # No path found
  
  shortest_path <- valid_paths[[which.min(lengths(valid_paths))]]  # Select the shortest one
  closest_ancestor <- tail(shortest_path, n = 1)$name  # Get last node (Level 2 term)
  
  return(closest_ancestor)
}

################################################################################
#####################  Load Excel Data Files  ##################################
################################################################################

# Define the Excel file path
excel_file <- file.path(base_dir, "UniqueStudy_VOIDs.xlsx")

# Load study data
figure_data <- read_excel(excel_file, sheet = "Data for Figure") %>%
  rename(VO_ID = 1, Vaccine = 2, Study_num = 3, Pathogen_Type = 4, Disease = 5)

immport_data <- read_excel(excel_file, sheet = "ImmPortSTUDY_VO_Simple") %>%
  rename(VO_ID = 1, STUDY_ACCESSION = 2, VO_Label = 3)

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
#####################  Map Each VO_ID to Its Level 2 Ancestor ##################
################################################################################

# Get VO terms from the Excel file
vo_terms_from_excel <- unique(figure_data$VO_ID)

# Apply function to find closest Level 2 ancestor
vo_hierarchy_mapping <- data.frame(
  VO_ID = vo_terms_from_excel,
  Level2_Ancestor = sapply(vo_terms_from_excel, find_closest_level2),
  stringsAsFactors = FALSE
)

# Merge hierarchy mapping into study data
explicit_counts <- explicit_counts %>%
  left_join(vo_hierarchy_mapping, by = "VO_ID")

# Print mapping summary
cat("\nHierarchy Mapping Summary:\n")
cat("Total VO terms in Excel:", length(vo_terms_from_excel), "\n")
cat("Successfully mapped to Level 2 terms:", sum(!is.na(vo_hierarchy_mapping$Level2_Ancestor)), "\n")

################################################################################
#####################  Calculate Implicit Study Counts  ########################
################################################################################

calculate_implicit_counts <- function(explicit_counts) {
  all_counts <- explicit_counts %>%
    mutate(implicit_study_num = explicit_study_num)
  
  # Aggregate implicit study counts by Level 2 ancestor
  implicit_counts <- explicit_counts %>%
    group_by(Level2_Ancestor) %>%
    summarise(
      implicit_study_num = sum(explicit_study_num, na.rm = TRUE),
      STUDY_ACCESS_IDs = paste(unique(STUDY_ACCESS_IDs), collapse=", "),
      .groups = "drop"
    ) %>%
    rename(VO_ID = Level2_Ancestor)
  
  # Merge implicit counts back into main dataset
  all_counts <- all_counts %>%
    left_join(implicit_counts, by = "VO_ID", suffix = c("_explicit", "_implicit"))
  
  return(all_counts)
}

# Compute implicit counts
final_counts <- calculate_implicit_counts(explicit_counts)

# Print summary
cat("\nFinal Study Count Summary:\n")
cat("Total explicit study terms:", nrow(explicit_counts), "\n")
cat("Total implicit study terms:", nrow(final_counts) - nrow(explicit_counts), "\n")










































################################################################################
#####################  Generate Treemaps for Each Pathogen Type  ###############
################################################################################

plot_pathogen_treemaps <- function(final_counts) {
  # Get unique pathogen types (excluding NA and Multiple)
  pathogen_types <- unique(final_counts$Pathogen_Type)
  pathogen_types <- pathogen_types[!is.na(pathogen_types) & pathogen_types != "Multiple"]
  
  # Create a treemap for each pathogen type
  for (pathogen in pathogen_types) {
    # Filter data for this pathogen
    pathogen_data <- final_counts %>%
      filter(Pathogen_Type == pathogen) %>%
      # Use implicit count for visualization
      mutate(study_count = implicit_study_num) %>%
      # Filter out zero counts
      filter(study_count > 0)
    
    # Skip if no data
    if (nrow(pathogen_data) == 0) {
      cat("No data for pathogen type:", pathogen, "\n")
      next
    }
    
    cat("Creating treemap for", pathogen, "with", nrow(pathogen_data), "entries\n")
    
    # Create treemap
    p <- ggplot(pathogen_data, 
                aes(area = study_count, 
                    fill = VO_Label, 
                    label = paste0(VO_Label, "\n(", study_count, ")"))) +
      geom_treemap() +
      geom_treemap_text(color = "black", place = "centre", 
                        grow = TRUE, min.size = 8, max.size = 14) +