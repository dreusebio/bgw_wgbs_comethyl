#First Conda Activate Comethyl before running this script
#Make sure to change to the directory of interest before runninng this script 

#This script will run enrichment analysis on the DMRs that were identified in the previous steps
#Make sure that your DMRs_annotated.xlsx file is in the working directory and that those DMRs have been annotated with gene symbols.
# The DMRs should have a column named 'geneSymbol' and a column named '
# These will be grouped together and separated by directionality
# Load required packages
#run with the command below
#Rscript /quobyte/lasallegrp/George/ECHO_NBDS/DM_cyto/cytosine_reports_combined/Enrichment.R
#the first block chunk will do it seperate by hyper andhypo. Then next is for all.

library(enrichR)
library(dplyr)
library(readxl)
library(openxlsx)
library(glue)
library(ggplot2)
library(viridis)
library(forcats)

# Set up Enrichr
setEnrichrSite("Enrichr")
enrichr_libraries <- c(
  "GO_Biological_Process_2025",
  "GO_Cellular_Component_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2021_Human",
  "Reactome_2022",
  "SynGO_2022"
)
#  "Panther_2016",
# Create output folder
dir.create("EnrichR", showWarnings = FALSE)

#all gene names considered
#I added this part 11/10/2025 to do with reading in the files
# --- Find the DMR annotated file (one main file, flexible naming) ---
dmr_file <- list.files(pattern = "^DMRs_annotated.*\\.xlsx$", ignore.case = TRUE)

if (length(dmr_file) == 0) {
  stop("❌ No file found starting with 'DMRs_annotated' in this directory.")
}

if (length(dmr_file) > 1) {
  warning("⚠️ Multiple files found. Using the first one: ", dmr_file[1])
  dmr_file <- dmr_file[1]
}

cat("📂 Reading file:", dmr_file, "\n")

# --- Read the file ---
dmrs <- read_excel(dmr_file) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Add file name for reference
dmrs$SourceFile <- sub("\\.xlsx$", "", basename(dmr_file))

cat("✅ DMRs loaded successfully — rows:", nrow(dmrs), "columns:", ncol(dmrs), "\n")

# Check structure
head(dmrs)

# Filter out NA gene symbols
genes <- unique(na.omit(dmrs$geneSymbol))

# Optionally, limit to protein-coding genes or certain direction
# genes <- unique(na.omit(dmr_data$geneSymbol[dmr_data$direction == "+"]))

# Check number of genes
cat("Number of unique genes:", length(genes), "\n")

# Define Enrichr libraries
enrichr_libraries <- c(
  "GO_Biological_Process_2025",
  "GO_Cellular_Component_2025",
  "GO_Molecular_Function_2025",
  "KEGG_2021_Human",
  "Reactome_2022",
  "SynGO_2022"
)
  #"Panther_2016",
# Run enrichment
enrichr_results <- enrichr(genes, enrichr_libraries)

# Save results to Excel file
wb <- createWorkbook()
for (i in seq_along(enrichr_results)) {
  df <- enrichr_results[[i]]
  sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", names(enrichr_results)[i]), 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, x = df)
}
saveWorkbook(wb, file = "EnrichR/DMR_EnrichR_Results.xlsx", overwrite = TRUE)


# Define plotting function: always plot top 25 terms
plot_and_save <- function(df, db_name) {
  if (is.null(df) || nrow(df) == 0) {
    cat("No data to plot for", db_name, "\n")
    return()
  }

  # Sort and select top 25 (by Adjusted.P.value, fallback to P.value)
  if ("Adjusted.P.value" %in% colnames(df)) {
    df <- df %>% arrange(Adjusted.P.value)
  } else {
    df <- df %>% arrange(P.value)
  }
  df_top25 <- df %>% slice_head(n = 25)

  # Plot and save
  pdf(glue("EnrichR/DMR_enrichr_plot_{gsub('[^A-Za-z0-9]', '_', db_name)}.pdf"), width = 15, height = 7)
  print(
    plotEnrich(df_top25, showTerms = 25, numChar = 75, y = "Count", orderBy = "P.value") +
      ggtitle(glue("Top 25 Enrichment Terms - {db_name}"))
  )
  dev.off()
}

# Generate plots for each database
for (db in names(enrichr_results)) {
  df <- enrichr_results[[db]]
  plot_and_save(df, db)
}


# plot significant terms for adjusted p-value and raw p-value if adjusted p-value is not available
# Plotting function with fallback
plot_and_save <- function(df, db_name) {
  if (is.null(df) || nrow(df) == 0) {
    cat("No data to plot for", db_name, "\n")
    return()
  }

  # First try adjusted p-value
  df_filtered <- df %>% filter(Adjusted.P.value <= 0.05)

  # Fallback to raw p-value if too few terms
  if (nrow(df_filtered) < 5) {
    cat("Too few adjusted terms for", db_name, "- using raw P.value instead\n")
    df_filtered <- df %>% filter(P.value <= 0.05)
  }

  if (nrow(df_filtered) == 0) {
    cat("No significant terms to plot for", db_name, "\n")
    return()
  }

  # Save PDF
  pdf(glue("EnrichR/DMR_enrichr_plot_{gsub('[^A-Za-z0-9]', '_', db_name)}.pdf"), height = 7, width = 15)
  print(
    plotEnrich(df_filtered, showTerms = 25, numChar = 75, y = "Count", orderBy = "P.value") +
      ggtitle(glue("DMR Enrichment - {db_name}"))
  )
  dev.off()
}

# Run plotting
for (db in names(enrichr_results)) {
  df <- enrichr_results[[db]]
  plot_and_save(df, db)
}

# plot only significant terms for adjusted p-value
# Plot top 25 terms for each database
for (db in names(enrichr_results)) {
  df <- enrichr_results[[db]]
  df <- df %>%
    filter(Adjusted.P.value <= 0.05) %>%
    arrange(Adjusted.P.value) %>%
    head(25)

  if (nrow(df) > 0) {
    p <- ggplot(df, aes(x = reorder(Term, -log10(Adjusted.P.value)), 
                        y = -log10(Adjusted.P.value), 
                        size = Combined.Score, 
                        fill = Adjusted.P.value)) +
      geom_point(shape = 21) +
      coord_flip() +
      scale_fill_viridis(option = "D", direction = -1, name = "Adj. P-value") +
      labs(
        title = paste("Top 25 Enriched Terms in", db),
        x = "Term",
        y = "-log10(Adjusted P-value)"
      ) +
      theme_minimal(base_size = 12)

    ggsave(filename = glue("EnrichR/DMR_EnrichR_{db}_dotplot.pdf"), plot = p, width = 10, height = 8)
  }
}

# # Read DMR file
# dmrs <- readxl::read_excel("DMRs_annotated.xlsx")

library(readxl)
library(dplyr)
#I added this part 11/10/2025
# --- Find the DMR annotated file (one main file, flexible naming) ---
dmr_file <- list.files(pattern = "^DMRs_annotated.*\\.xlsx$", ignore.case = TRUE)

if (length(dmr_file) == 0) {
  stop("❌ No file found starting with 'DMRs_annotated' in this directory.")
}

if (length(dmr_file) > 1) {
  warning("⚠️ Multiple files found. Using the first one: ", dmr_file[1])
  dmr_file <- dmr_file[1]
}

cat("📂 Reading file:", dmr_file, "\n")

# --- Read the file ---
dmrs <- read_excel(dmr_file) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Add file name for reference
dmrs$SourceFile <- sub("\\.xlsx$", "", basename(dmr_file))

cat("✅ DMRs loaded successfully — rows:", nrow(dmrs), "columns:", ncol(dmrs), "\n")


# Filter valid gene symbols and separate by direction
dmrs <- dmrs %>% filter(!is.na(geneSymbol), geneSymbol != "")
dmr_groups <- split(dmrs, dmrs$direction)


run_enrichr_analysis <- function(genes, label) {
  message("Running Enrichr for: ", label)
  
  genes_unique <- unique(genes)
  if (length(genes_unique) == 0) {
    warning("No genes found for ", label)
    return(NULL)
  }
  
  results <- enrichr(genes_unique, enrichr_libraries)
  wb <- createWorkbook()
  dir.create("EnrichR", showWarnings = FALSE)

  all_results <- list()
  
  for (lib in names(results)) {
    df <- results[[lib]]
    df <- df %>% mutate(Database = lib, Direction = label)
    
    # Extract gene count from "Overlap" column
    if ("Overlap" %in% colnames(df)) {
      df <- df %>%
        mutate(Count = as.numeric(sub("/.*", "", Overlap)))
    }
    
    all_results[[lib]] <- df
    addWorksheet(wb, lib)
    writeData(wb, lib, df)
    
    # Select top 25 with P.value < 0.05 and valid Count
    top_df <- df %>%
      filter(!is.na(Count), !is.na(P.value), P.value < 0.05) %>%
      arrange(desc(Count)) %>%
      head(25)
    
    if (nrow(top_df) > 0) {
      p <- ggplot(top_df, aes(x = reorder(Term, P.value), y = Count, fill = P.value)) +
        geom_col() +
        coord_flip() +
        scale_fill_gradient(low = "red", high = "blue", name = "P-value", trans = "log10") +
        labs(
          title = paste0("Top 25 Enriched Terms - ", label, " - ", lib),
          x = "Enriched Term",
          y = "Gene Count"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
          axis.text.y = element_text(size = 10),
          plot.title = element_text(face = "bold", size = 14),
          legend.position = "right"
        )

      ggsave(
        filename = file.path("EnrichR", paste0(label, "_", gsub("[^A-Za-z0-9]", "_", lib), "_trial.pdf")),
        plot = p, width = 12, height = 7
      )
    }
  }
  
  saveWorkbook(wb, file.path("EnrichR", paste0(label, "_Enrichr.xlsx")), overwrite = TRUE)
  return(bind_rows(all_results))
}

# this code makes bargraphs based on the p values, not used anymore.
# # Function to run enrichment and save plots
# run_enrichr_analysis <- function(genes, label) {
#   message("Running Enrichr for: ", label)
  
#   genes_unique <- unique(genes)
#   if (length(genes_unique) == 0) {
#     warning("No genes found for ", label)
#     return(NULL)
#   }
  
#   results <- enrichr(genes_unique, enrichr_libraries)
#   wb <- createWorkbook()
  
#   for (lib in names(results)) {
#     df <- results[[lib]]
#     df <- df %>% mutate(Database = lib, Direction = label)
#     addWorksheet(wb, lib)
#     writeData(wb, lib, df)
    
#     # Always plot top 25 (even if not significant)
#     top_df <- df %>%
#       arrange(P.value) %>%
#       head(25)
    
#     if (nrow(top_df) > 0) {
#       p <- ggplot(top_df, aes(x = reorder(Term, -P.value), y = -log10(P.value))) +
#         geom_col(fill = "steelblue") +
#         coord_flip() +
#         labs(title = paste0("Top 25 Terms - ", label, " - ", lib),
#              x = "", y = "-log10(P-value)") +
#         theme_minimal(base_size = 10)
#       ggsave(filename = file.path("EnrichR", paste0(label, "_", gsub("[^A-Za-z0-9]", "_", lib), ".pdf")),
#              plot = p, width = 10, height = 6)
#     }
#   }
  
#   saveWorkbook(wb, file.path("EnrichR", paste0(label, "_Enrichr.xlsx")), overwrite = TRUE)
#   return(bind_rows(results))
# }



# Run enrichment for each direction
# Run enrichment by direction
# all_enrichment <- lapply(names(dmr_groups), function(dir) {
#   run_enrichr_analysis(dmr_groups[[dir]]$geneSymbol, dir)
# }) %>% bind_rows()


#This worked.
#THis will be the script to run enrichment analysis on the DMRs that were identified in the previous steps
# These are grouped togeether and not separated by directionality
# # ===========================
# # Enrichment Analysis from DMRichR Output
# # ===========================

# # Load required libraries
# library(enrichR)
# library(readxl)
# library(openxlsx)
# library(dplyr)
# library(ggplot2)
# library(glue)
# library(viridis)

# # Set Enrichr site
# listEnrichrSites()
# setEnrichrSite("Enrichr") # Human genes

# # Create output directory
# dir.create("EnrichR", showWarnings = FALSE)
