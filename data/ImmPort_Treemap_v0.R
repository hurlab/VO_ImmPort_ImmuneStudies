################################################################################
#####################  Set Current Working Directory ###########################
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
  "rdflib", "igraph", "readxl", "ggplot2", "treemap", "dplyr", "tidyr", "xml2"
)

# Automatically install and load CRAN packages
p_load(char = cran_packages, install = TRUE)

################################################################################
#####################  Load and Parse OWL File ##################################
################################################################################

# Define the OWL file path (update accordingly)
owl_file <- file.path(base_dir, "your_ontology_file.owl")

# Load OWL file as RDF
rdf <- rdf_parse(owl_file)

# Extract Parent-Child Relationships with Explicit Namespace Declaration
get_ontology_relations <- function(rdf) {
  query <- '
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    SELECT ?child ?parent WHERE {
      ?child rdfs:subClassOf ?parent .
    }'
  
  triples <- rdf_query(rdf, query)
  
  if (is.null(triples) || nrow(triples) == 0) {
    stop("No subclass relationships found. Ensure the OWL file is structured correctly.")
  }
  
  return(triples)
}

ontology_relations <- get_ontology_relations(rdf)
colnames(ontology_relations) <- c("Child", "Parent")

# Convert URIs to readable VO IDs
ontology_relations <- ontology_relations %>%
  mutate(
    Child = gsub(".*/", "", Child),
    Parent = gsub(".*/", "", Parent)
  )

################################################################################
#####################  Load Vaccine Data  #######################################
################################################################################

# Define the vaccine data file path
file_path <- file.path(base_dir, "20250303_immune_exposure_studies.xlsx")
sheet_name <- "Data for Figure"

# Load the vaccine dataset
vaccine_data <- read_excel(file_path, sheet = sheet_name)

# Ensure correct column names
colnames(vaccine_data) <- c("VO_ID", "Vaccine", "Study_num", "Pathogen_Type", "Disease")

################################################################################
#####################  Propagate Study_num PER Pathogen Type ###################
################################################################################

# Convert ontology into an igraph object for hierarchical processing
graph <- graph_from_data_frame(ontology_relations, directed = TRUE)

# Function to propagate study counts from children to parents within the same pathogen type
propagate_study_numbers_by_pathogen <- function(vaccine_data, ontology_relations) {
  # Initialize result dataframe
  adjusted_data <- vaccine_data
  
  # Process each pathogen type independently
  pathogen_types <- unique(vaccine_data$Pathogen_Type)
  
  for (pathogen in pathogen_types) {
    cat("Processing pathogen type:", pathogen, "\n")
    
    # Filter vaccine data for this pathogen type
    subset_data <- vaccine_data %>% filter(Pathogen_Type == pathogen)
    
    # Filter ontology relationships for only relevant vaccines
    relevant_ontology <- ontology_relations %>%
      filter(Child %in% subset_data$VO_ID | Parent %in% subset_data$VO_ID)
    
    # Initialize study counts as a named vector
    study_map <- setNames(subset_data$Study_num, subset_data$VO_ID)
    
    # Iterate over parent-child relationships within the same pathogen type
    for (i in seq_len(nrow(relevant_ontology))) {
      parent <- relevant_ontology$Parent[i]
      child <- relevant_ontology$Child[i]
      
      if (!is.na(study_map[child])) {
        study_map[parent] <- sum(study_map[parent], study_map[child], na.rm = TRUE)
      }
    }
    
    # Convert back to a data frame
    adjusted_study_counts <- data.frame(VO_ID = names(study_map), Study_num_Adjusted = study_map, row.names = NULL)
    
    # Merge with original data
    adjusted_data <- adjusted_data %>%
      left_join(adjusted_study_counts, by = "VO_ID")
  }
  
  return(adjusted_data)
}

vaccine_data_adjusted <- propagate_study_numbers_by_pathogen(vaccine_data, ontology_relations)

################################################################################
#####################  Generate Treemap for Each Pathogen Type #################
################################################################################

# Create treemaps for each Pathogen Type
plot_treemap <- function(vaccine_data_adjusted) {
  pathogen_types <- unique(vaccine_data_adjusted$Pathogen_Type)
  
  for (pathogen in pathogen_types) {
    data_subset <- vaccine_data_adjusted %>% filter(Pathogen_Type == pathogen)
    
    # Remove NA values
    data_subset <- data_subset %>% drop_na(Study_num_Adjusted)
    
    if (nrow(data_subset) > 0) {
      treemap(
        data_subset,
        index = "Vaccine",
        vSize = "Study_num_Adjusted",
        title = paste("Vaccine Study Counts -", pathogen),
        palette = "Blues",
        fontsize.labels = 12,
        fontface.labels = "bold",
        border.col = "white"
      )
    }
  }
}

# Generate treemaps
plot_treemap(vaccine_data_adjusted)
