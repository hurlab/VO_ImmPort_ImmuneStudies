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
  "rdflib", "readxl", "ggplot2", "treemapify", "dplyr", "tidyr", "xml2", "data.table"
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

# Function to detect correct ontology relationship type
detect_subclass_property <- function(rdf) {
  properties <- rdf_query(rdf, 'SELECT DISTINCT ?p WHERE { ?s ?p ?o }')
  
  if ("http://www.w3.org/2000/01/rdf-schema#subClassOf" %in% properties$p) {
    return("rdfs:subClassOf")
  } else if ("http://www.w3.org/2002/07/owl#subClassOf" %in% properties$p) {
    return("owl:subClassOf")
  } else {
    stop("No subclass relationships found. Check the OWL file structure.")
  }
}

# Detect the correct property
subclass_property <- detect_subclass_property(rdf)

# Function to extract only relevant parent-child relationships
get_relevant_ontology_relations <- function(rdf, subclass_property, vaccine_data) {
  query <- sprintf('
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    SELECT ?child ?parent WHERE { ?child %s ?parent . }', subclass_property)
  
  triples <- rdf_query(rdf, query)
  
  if (is.null(triples) || nrow(triples) == 0) {
    stop("No subclass relationships found. Ensure the OWL file is structured correctly.")
  }
  
  # Format VO terms
  colnames(triples) <- c("Child", "Parent")
  triples <- triples %>%
    mutate(
      Child = gsub(".*/", "", Child),
      Parent = gsub(".*/", "", Parent)
    )
  
  # Filter only the terms that exist in the vaccine dataset
  relevant_terms <- unique(vaccine_data$VO_ID)
  triples <- triples %>%
    filter(Child %in% relevant_terms | Parent %in% relevant_terms)
  
  return(triples)
}

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
#####################  Extract Only Relevant Ontology Terms ####################
################################################################################

# Extract parent-child relationships **ONLY for VO terms in the vaccine dataset**
ontology_relations <- get_relevant_ontology_relations(rdf, subclass_property, vaccine_data)

################################################################################
#####################  Propagate Study_num Within Given Terms ##################
################################################################################

propagate_study_numbers <- function(vaccine_data, ontology_relations) {
  vaccine_data <- as.data.table(vaccine_data)
  ontology_relations <- as.data.table(ontology_relations)
  
  # Merge **only relevant parent-child relationships**
  for (i in 1:3) {  # Limit to 3 levels to keep it fast
    vaccine_data <- merge(vaccine_data, ontology_relations, by.x = "VO_ID", by.y = "Child", all.x = TRUE)
    
    # Replace VO_ID with Parent if applicable
    vaccine_data <- vaccine_data %>%
      mutate(VO_ID = ifelse(is.na(Parent), VO_ID, Parent)) %>%
      select(-Parent)
  }
  
  return(vaccine_data)
}

# Apply study number propagation
vaccine_data_adjusted <- propagate_study_numbers(vaccine_data, ontology_relations)

################################################################################
#####################  Aggregate Study Numbers Without Loops ###################
################################################################################

aggregate_study_numbers <- function(vaccine_data_adjusted) {
  vaccine_data_adjusted %>%
    group_by(VO_ID, Pathogen_Type) %>%
    summarise(Study_num = sum(Study_num, na.rm = TRUE), .groups = "drop")
}

# Apply aggregation
vaccine_data_final <- aggregate_study_numbers(vaccine_data_adjusted)

################################################################################
#####################  Generate Treemap for Each Pathogen Type #################
################################################################################

plot_nested_treemap <- function(vaccine_data_final) {
  pathogen_types <- unique(vaccine_data_final$Pathogen_Type)
  
  for (pathogen in pathogen_types) {
    data_subset <- vaccine_data_final %>% filter(Pathogen_Type == pathogen)
    
    if (nrow(data_subset) > 0) {
      ggplot(data_subset, aes(area = Study_num, fill = VO_ID, label = VO_ID)) +
        geom_treemap() +
        geom_treemap_text(color = "black", place = "centre", size = 5) +
        theme(legend.position = "none") +
        labs(title = paste("Vaccine Study Counts -", pathogen))
    }
  }
}

# Generate hierarchical treemap
plot_nested_treemap(vaccine_data_final)
