library(data.table)

library(ggfortify)
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop("Package 'openxlsx' not installed. Install with:\nRscript -e 'install.packages(\"openxlsx\", repos=\"https://cloud.r-project.org\")'")
}

library(ggfortify)
library(openxlsx)

# Function to create PCA plot from smoothed methylation data
SmoothMethPCAPlot <- function(file, sample_info, name){
  x <- read.delim(file, header = TRUE, check.names = FALSE)
  x <- x[,-1:-16]
  t_x <- data.table::transpose(x)
  colnames(t_x) <- rownames(x)
  rownames(t_x) <- colnames(x)
  sample_info <- read.xlsx(sample_info, rowNames = TRUE)
  y <- merge(sample_info, t_x, by = 'row.names', all = FALSE)
  PCA <- prcomp(y[,-1:-3], center = TRUE, scale. = TRUE)
  plot <- autoplot(PCA, data = y, colour = "Diagnosis", frame = TRUE, frame.type = 'norm') +
                   theme_classic() + 
    theme(plot.background = element_rect(fill = "white"),
          axis.text = element_text(color = "black", size = 12, family = "sans"),
          axis.title.x = element_text(color = "black", size = 12, family = "sans"),
          axis.title.y = element_text(color = "black", size = 12, family = "sans"),
          axis.line = element_line(color = "black", linewidth = 0.2),
          legend.position = "right")
  pdf(file = name, width = 7, height = 6)
  print(plot)
  dev.off()
  return(plot)
}

# # Run PCA plot for environmental risk
dmr_text_file <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/02_outputs/females_subtypes1_and_2_methylLearn/ASD_subtype1_vs_Control1/topPercent_15_V3_004/DMR_individual_smoothed_methylation_filtered.txt"
sample_info <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/00_bismark_data/ECHO_subtypes_females/ECHO_Idiopathic_ASD_subtype1_vs_Control_subtype1/V3_004/sample_info.xlsx"

ECHO_samples_CHDS_DMRs_combined <- SmoothMethPCAPlot(
  file = dmr_text_file,
  sample_info = sample_info,
  name = "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/01_scripts/01_R_scripts/ECHO_Females_Subtype1_Filtered_PCA_R.pdf"
)

# Run PCA plot for 0to6
dmr_text_file <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/02_outputs/females_subtypes1_and_2_methylLearn/ASD_subtype2_vs_Control1/topPercent_15_V4_0038/DMR_individual_smoothed_methylation_filtered.txt"
sample_info <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/00_bismark_data/ECHO_subtypes_females/ECHO_Idiopathic_ASD_subtype2_vs_Control_subtype1/V4_0038/sample_info.xlsx"

ECHO_samples_CHDS_DMRs_combined <- SmoothMethPCAPlot(
  file = dmr_text_file,
  sample_info = sample_info,
  name = "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/01_scripts/01_R_scripts/ECHO_Females_Subtype2_Filtered_PCA_R.pdf"
)

# Run PCA plot for 6to12
dmr_text_file <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/02_outputs/females_subtypes1_and_2_methylLearn/ASD_subtype2_vs_Control1/topPercent_15_V4_0038/DMR_individual_smoothed_methylation_filtered.txt"
sample_info <- "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/00_bismark_data/ECHO_subtypes_females/ECHO_Idiopathic_ASD_subtype2_vs_Control_subtype1/V4_0038/sample_info.xlsx"

ECHO_samples_CHDS_DMRs_combined <- SmoothMethPCAPlot(
  file = dmr_text_file,
  sample_info = sample_info,
  name = "/quobyte/lasallegrp/Giovani/ECHO_Machine_Learning/01_scripts/01_R_scripts/ECHO_Females_Subtype2_Filtered_PCA_R.pdf"
)