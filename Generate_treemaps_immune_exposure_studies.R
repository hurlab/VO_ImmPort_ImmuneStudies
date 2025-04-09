################################################################################
##################### Set Current Working Directory ###########################
################################################################################

if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")
library(rstudioapi)

current_path <- rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path))
print(current_path)
base_dir <- dirname(current_path)

output_dir <- file.path(base_dir, "output")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
data_dir <- file.path(base_dir, "data")

################################################################################
##################### Packages #################################################
################################################################################

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
library(pacman)

p_load(char = c(
  "ontologyIndex", "ontologySimilarity", "readxl", "dplyr",
  "tidyr", "data.table", "ggplot2", "treemapify"
))

################################################################################
##################### Load Ontology from OBO File ##############################
################################################################################

obo_file <- file.path(data_dir, "vo.obo")
ontology <- get_OBO(obo_file, extract_tags = "everything")
if (!"VO:0000001" %in% ontology$id) stop("VO:0000001 not found in ontology.")

################################################################################
##################### Load Study Data from Excel ###############################
################################################################################

excel_file <- file.path(data_dir, "UniqueStudy_VOIDs.xlsx")

figure_data <- read_excel(excel_file, sheet = "Data for Figure") %>%
  rename(VO_ID = 1, Vaccine = 2, Study_num = 3, Pathogen_Type = 4, Disease = 5)

immport_data <- read_excel(excel_file, sheet = "ImmPortSTUDY_VO_Simple") %>%
  rename(VO_ID = 1, STUDY_ACCESSION = 2, VO_Label = 3)

convert_vo_format <- function(vo_id) gsub("VO_", "VO:", vo_id)
figure_data$VO_ID <- sapply(figure_data$VO_ID, convert_vo_format)
immport_data$VO_ID <- sapply(immport_data$VO_ID, convert_vo_format)

vo_map <- figure_data %>%
  select(VO_ID, Pathogen_Type, Disease) %>%
  filter(complete.cases(.)) %>%
  distinct()

merged_data <- immport_data %>%
  left_join(vo_map, by = "VO_ID") %>%
  filter(complete.cases(.)) %>%
  distinct()

explicit_counts <- merged_data %>%
  group_by(VO_ID, VO_Label, Pathogen_Type, Disease) %>%
  summarise(
    STUDY_ACCESS_IDs = paste(sort(unique(STUDY_ACCESSION)), collapse = ", "),
    explicit_study_num = n_distinct(STUDY_ACCESSION),
    .groups = "drop"
  )

################################################################################
##################### Labeling Helper Functions ################################
################################################################################

get_vo_label <- function(vo_id) {
  label_index <- which(ontology$id == vo_id)
  if (length(label_index) > 0) ontology$name[label_index] else vo_id
}

truncate_label <- function(label, vo_id, count, max_length = 60) {
  base_info <- paste0(" (", vo_id, "; n=", count, ")")
  max_name_length <- max_length - nchar(base_info)
  if (nchar(label) > max_name_length) {
    short_label <- substr(label, 1, max_name_length - 1)
    label <- paste0(short_label, "~")
  }
  paste0(label, base_info)
}

################################################################################
##################### Treemap Plotting Function ################################
################################################################################

plot_treemap_versions_pdf <- function(data, count_column, output_dir, prefix = "Treemap") {
  pathogen_types <- unique(data$Pathogen_Type)
  
  for (pathogen in pathogen_types) {
    subset_data <- data %>%
      filter(Pathogen_Type == pathogen & !!sym(count_column) > 0) %>%
      mutate(
        Label_Simple = VO_Label,
        Label_Detailed = mapply(truncate_label, VO_Label, VO_ID, .data[[count_column]])
      )
    
    if (nrow(subset_data) > 0) {
      safe_pathogen <- gsub("[^A-Za-z0-9]", "_", pathogen)
      
      # Simple Label
      pdf(file = file.path(output_dir, paste0(prefix, "_", safe_pathogen, "_Simple.pdf")), width = 10, height = 8)
      p1 <- ggplot(subset_data, aes(area = !!sym(count_column), fill = VO_Label, label = Label_Simple)) +
        geom_treemap() +
        geom_treemap_text(place = "centre", grow = TRUE, reflow = TRUE, color = "black") +
        theme(legend.position = "none") +
        labs(title = paste("Study Counts -", pathogen, "(Label: VO_Label)"))
      print(p1)
      dev.off()
      
      # Detailed Label
      pdf(file = file.path(output_dir, paste0(prefix, "_", safe_pathogen, "_Detailed.pdf")), width = 12, height = 10)
      p2 <- ggplot(subset_data, aes(area = !!sym(count_column), fill = VO_Label, label = Label_Detailed)) +
        geom_treemap() +
        geom_treemap_text(place = "centre", grow = TRUE, reflow = TRUE, color = "black") +
        theme(legend.position = "none") +
        labs(title = paste("Study Counts -", pathogen, "(Detailed Label)"))
      print(p2)
      dev.off()
    }
  }
}

################################################################################
##################### Simplify Study Data for Aggregation ######################
################################################################################

vo_terms_file <- file.path(data_dir, "VO_term_updated_250403.txt")
vo_term_list <- read.delim(vo_terms_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
rownames(vo_term_list) <- gsub("VO_", "VO:", rownames(vo_term_list))
allowed_vo_ids <- rownames(vo_term_list)

level2_IDs <- ontology$children[["VO:0000001"]]
level3_IDs <- unlist(lapply(level2_IDs, function(term) ontology$children[[term]]))
level2_IDs_filtered <- intersect(level2_IDs, allowed_vo_ids)
level3_IDs_filtered <- intersect(level3_IDs, allowed_vo_ids)

create_level_aggregated <- function(filtered_level_IDs, study_data, ontology) {
  study_long <- study_data %>%
    select(VO_ID, Pathogen_Type, STUDY_ACCESS_IDs) %>%
    separate_rows(STUDY_ACCESS_IDs, sep = ", ") %>%
    filter(!is.na(STUDY_ACCESS_IDs))
  
  ancestor_mapping <- list()
  
  for (i in seq_len(nrow(study_long))) {
    vo_id <- study_long$VO_ID[i]
    study_id <- study_long$STUDY_ACCESS_IDs[i]
    pathogen_type <- study_long$Pathogen_Type[i]
    
    if (!vo_id %in% ontology$id) next
    ancestors <- get_ancestors(ontology, vo_id)
    ancestors <- ancestors[ancestors %in% filtered_level_IDs]
    
    for (ancestor in ancestors) {
      key <- paste(ancestor, pathogen_type, sep = "|")
      if (is.null(ancestor_mapping[[key]])) {
        ancestor_mapping[[key]] <- list(studies = character(), pathogen = pathogen_type)
      }
      ancestor_mapping[[key]]$studies <- c(ancestor_mapping[[key]]$studies, study_id)
    }
  }
  
  result_df <- do.call(rbind, lapply(names(ancestor_mapping), function(key) {
    parts <- strsplit(key, "\\|")[[1]]
    vo_id <- parts[1]
    pathogen <- parts[2]
    study_ids <- unique(ancestor_mapping[[key]]$studies)
    data.frame(
      VO_ID = vo_id,
      VO_Label = get_vo_label(vo_id),
      Pathogen_Type = pathogen,
      STUDY_ACCESS_IDs = paste(study_ids, collapse = ", "),
      implicit_study_num = length(study_ids),
      stringsAsFactors = FALSE
    )
  }))
  
  return(result_df)
}

################################################################################
##################### Generate & Plot Treemaps #################################
################################################################################

cat("\nAggregating implicit study data...\n")
level2_implicit_counts <- create_level_aggregated(level2_IDs_filtered, explicit_counts, ontology)
level3_implicit_counts <- create_level_aggregated(level3_IDs_filtered, explicit_counts, ontology)

cat("Level 2 Aggregated Terms:", nrow(level2_implicit_counts), "\n")
cat("Level 3 Aggregated Terms:", nrow(level3_implicit_counts), "\n")

cat("\nGenerating Treemaps for Explicit Study Counts...\n")
plot_treemap_versions_pdf(explicit_counts, "explicit_study_num", output_dir, prefix = "Treemap_Explicit")

cat("\nGenerating Treemaps for Implicit Study Counts...\n")
plot_treemap_versions_pdf(level2_implicit_counts, "implicit_study_num", output_dir, prefix = "Treemap_Implicit_Level2")
plot_treemap_versions_pdf(level3_implicit_counts, "implicit_study_num", output_dir, prefix = "Treemap_Implicit_Level3")



################################################################################
##################### Disease Summary Treemap ##################################
################################################################################

plot_disease_summary_treemap <- function(data, output_dir, prefix = "Treemap_Disease_Summary") {
  disease_summary <- data %>%
    group_by(Disease) %>%
    summarise(
      Study_num = n_distinct(unlist(strsplit(STUDY_ACCESS_IDs, ",\\s*"))),
      VO_Count = n_distinct(VO_ID),
      .groups = "drop"
    ) %>%
    mutate(
      Label = paste0(Disease, " (", VO_Count, " VOs with ", Study_num, " studies)")
    )
  
  pdf(file = file.path(output_dir, paste0(prefix, ".pdf")), width = 12, height = 10)
  p <- ggplot(disease_summary, aes(area = Study_num, fill = Disease, label = Label)) +
    geom_treemap() +
    geom_treemap_text(color = "black", place = "centre", grow = TRUE, reflow = TRUE) +
    theme(legend.position = "none") +
    labs(title = "Treemap of Diseases by Total Study Count and VO Coverage")
  print(p)
  dev.off()
}

cat("\nDisease summary treemap...\n")
plot_disease_summary_treemap(explicit_counts, output_dir)


