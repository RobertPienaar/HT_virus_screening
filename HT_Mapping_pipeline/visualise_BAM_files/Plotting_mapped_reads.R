# Heatmap BAM for BSF virus screening

rm(list=ls())

# The script

Input_directory_address <- c("path/to/working/directory")

name_of_output_heatmap <- c(paste("Heatmap_reads_virus_screening_", Sys.Date(), ".svg", sep = ""))


# Set the working directory for the script, you should have your input files here and this will also be the directory where the output will be stored.
getwd()
setwd(Input_directory_address)
getwd()


# Loading libraries 
library("Rsamtools")
library("ComplexHeatmap")
library("dplyr")

library("GenomicRanges")
library("GenomicAlignments")

# Path to BAM files

# Virus being screened ----
Virus_name <- "HiTV2"
path_to_input_bam_directory <- "path/to/bam_files"

#bam_files <- list.files(path = "path/to/bam_files", pattern = "\\.bam$", full.names = TRUE, recursive = TRUE)
bam_files <- list.files(path = path_to_input_bam_directory, pattern = "\\.bam$", full.names = TRUE, recursive = TRUE)


#----------------------------------------------------------------------------------------------------------------------------------

#Note: if the warning "[W::hts_idx_load3] The index file is older than the data file:...*.bai" appears, let the script complete successfully before troubleshooting this, it is not necessarily an issue.

{

# Dimensions for output images
# 
wdth <- 2500
hght <- 2773

# Dimensions 131 samples
# For PNG
wdth_png <- 2500
hght_png <- 8500

# For SVG, the values here will be much smaller than for the PNG
wdth_svg=10
hght_svg=30

#######################



# Initialize list to store coverage values and mapped read counts for each sample
coverage_list <- list()
mapped_read_counts_list <- list()

# Loop through each BAM file
for (bam_file in bam_files) {
  # Read BAM file
  # Using readGAlignments to obtain the max length of the reference sequence automatically
  param_auto <- ScanBamParam(what = c("rname", "strand", "pos", "qwidth", "mapq"))
  aln_auto <- readGAlignments(bam_file, param = param_auto)
  
  # Obtaining the sequence name and max length automatically, and setting the coverage range
  Ref_seq_fasta_header <- aln_auto@seqinfo@seqnames
  max_length <- aln_auto@seqinfo@seqlengths
  cov_range <- 1:max_length
  
  # You can specify regions if needed, or leave it to use the max length and reference sequence name provided by "aln"
  regions <- GRanges(Ref_seq_fasta_header, IRanges(1,max_length))
  
 # Now the the reference sequence name and max length has been set, the script will reimport the bam file coverage and put it into a list container
  
  # Perform the pileup to calculate coverage
  # The `cov` object will contain the coverage information
  cov <- pileup(bam_file, scanBamParam=ScanBamParam(which=regions))
  
 # Assuming 'cov' is your pileup result
  # Create a data frame with all positions
  all_positions <- data.frame(pos = cov_range, cov_count = 0)
  
  
  # Rename 'count' column in 'cov'
  cov_df <- cov %>%
    rename(cov_count = count) # Change 'count' to 'cov_count'
  
  
  # Aggregate coverage by position
  cov_df_aggregated <- cov_df %>%
    group_by(pos) %>%
    summarise(cov_count = sum(cov_count, na.rm = TRUE)) %>%
    ungroup()  # Ensure we remove the grouping for further operations
  
  
  # Perform the left join, mutate the two resulting coverage count columns to replace the NAs with 0s and remove column "cov_count.x"
  cov_complete <- all_positions %>%
    left_join(cov_df_aggregated, by = "pos") %>%
    # Replace NA values in cov_count.y with 0
    mutate(cov_count = ifelse(is.na(cov_count.y), 0, cov_count.y)) %>%
    # Remove cov_count.x and cov_count.y, keep the new cov_count and other columns
    select(-cov_count.x, -cov_count.y)
  
  # Check the first few rows to confirm that the remaining coloumns are named correctly
  head(cov_complete)
  
  # Extracting and transposing the coverage values
  transposed_coverage <- t(cov_complete$cov_count)
  

  # Function to count mapped reads in a BAM file
  count_mapped_reads <- function(bam_file_path) {
    # Use scanBam to read flags from the BAM file
    bam_data <- scanBam(bam_file_path, param = ScanBamParam(what = "flag"))
    
    # Extract the flags and ensure they are integers
    flags <- unlist(lapply(bam_data, function(x) as.integer(x$flag)))
    
    # Count reads that are not marked as unmapped (0x4 flag is not set)
    sum(bitwAnd(flags, 0x4) == 0)
  }
  
 
    # Count mapped reads for the current BAM file
  mapped_reads <- count_mapped_reads(bam_file)
  
  # Store mapped read count and coverage values in the list using the sample name as key
  sample_name <- basename(bam_file)
  coverage_list[[sample_name]] <- as.integer(transposed_coverage)
  mapped_read_counts_list[[sample_name]] <- mapped_reads
  

}

# Converting the coverage list into a matrix for plotting
{
# Combine coverage values into a matrix
readMappingMatrix <- do.call(rbind, coverage_list)
  

####

# Get the current row names
current_row_names <- row.names(readMappingMatrix)

# Remove "_mappedsortedfiltered.bam" from the row names
cleaned_row_names <- gsub("_mappedsortedfiltered.bam", "", current_row_names)

# Set the cleaned row names
row.names(readMappingMatrix) <- cleaned_row_names

####
ncol(readMappingMatrix)



# Create the bottom_annotation without specifying annotation_name_side
# Calculate the total number of positions
total_positions <- ncol(readMappingMatrix)

# Initialize the vector with empty strings (or NA, based on preference)
position_annotations <- rep("", total_positions)

# Set the marker text every 1000 bases, starting from 1
marker_indices <- seq(1, total_positions, by = 1000)
position_annotations[marker_indices] <- as.character(marker_indices - 1)  # Adjust to start from 0 and then add 1

# Correct the first position manually if starting from 1 is preferred
position_annotations[1] <- "1"

# Create the bottom_annotation with text annotation
bottom_annotation <- HeatmapAnnotation(Position = anno_text(position_annotations, which = "column", just = "center", rot = 0))
}

#-------------------------------------------------------------

# Convert mapped_read_counts_list to a vector if it's not already
mapped_counts_vector <- unname(unlist(mapped_read_counts_list))
# Note: Ensure the order of `mapped_counts_vector` matches the order of rows in your heatmap matrix
row_annotation <- rowAnnotation(mapped_reads = anno_text(mapped_counts_vector, which = "row"))

row_annotation <- rowAnnotation(
  Counts = anno_text(mapped_counts_vector, which = "row", just = "left"),
  width = unit(2, "cm") # Adjust width as necessary
)

#--------------------------------------------------------------
# If you want to quickly determine if enough reads are observed to determine if the virus is likely there or not ----

above_10_coverage_matrix <- ifelse(readMappingMatrix > 10, 10, readMappingMatrix)

# Plotting image if needed to adjust resolution 

ht  <-  Heatmap(above_10_coverage_matrix, use_raster = TRUE,
                name = "Read Density", 
                row_names_side = "left", 
                column_names_side = "bottom",  # This controls where column names are placed relative to the heatmap
                show_row_names = TRUE, 
                show_column_names = TRUE,  # Assuming you want to show column names
                cluster_rows = FALSE, 
                cluster_columns = FALSE, 
                column_title = "Position (nt)", 
                row_title = "Samples",
                col = colorRampPalette(c("#FBFDED", "#009999", "#FF6600"))(10),
                heatmap_legend_param = list(title = "Read Hits", at = c(min(readMappingMatrix, na.rm = TRUE), max(readMappingMatrix, na.rm = TRUE)), labels = c("Low", "High")),
                bottom_annotation = bottom_annotation,  # Add your bottom_annotation here
                right_annotation = row_annotation  # Add the row annotation on the right
        )

draw(ht)


# Saving image as PNG
{
  png(file=paste(Input_directory_address,paste(Virus_name,"_Heatmap_virus_screening_above_10_reads_", Sys.Date(), ".png", sep = ""),sep="/"), width=wdth_png,height=hght_png,units="px",res=300)
  draw(ht)
  # Decorate the heatmap body with horizontal lines between rows
  ht_lines <-  decorate_heatmap_body("Read Density", {
    for (i in 1:(nrow(above_10_coverage_matrix) - 1)) {
      grid.lines(x = c(0, 1), y = rep(1 - i / nrow(above_10_coverage_matrix), 2), gp = gpar(lty = 1, lwd = 0.2))
    }
  })
  dev.off()
} 

# Saving image as SVG
{
  svg(file=paste(Input_directory_address,paste(Virus_name,"_Heatmap_virus_screening_above_10_reads_", Sys.Date(), ".svg", sep = ""),sep="/"), width=wdth_svg,height=hght_svg)
  draw(ht)
  # Decorate the heatmap body with horizontal lines between rows
  ht_lines <-  decorate_heatmap_body("Read Density", {
    for (i in 1:(nrow(above_10_coverage_matrix) - 1)) {
      grid.lines(x = c(0, 1), y = rep(1 - i / nrow(above_10_coverage_matrix), 2), gp = gpar(lty = 1, lwd = 0.2))
    }
  })
  dev.off()
} 

#-----------------------------------------------------------------------------------------------------------------


# If you want to have a general scope of the coverage across the reference genome ----

ht_max <- Heatmap(readMappingMatrix, use_raster = TRUE,
                  name = "Read Density", 
                  row_names_side = "left", 
                  column_names_side = "bottom",  # This controls where column names are placed relative to the heatmap
                  show_row_names = TRUE, 
                  show_column_names = TRUE,  # Assuming you want to show column names
                  cluster_rows = FALSE, 
                  cluster_columns = FALSE, 
                  column_title = "Position (nt)", 
                  row_title = "Samples",
                  col = colorRampPalette(c("white", "blue", "red"))(100),
                  heatmap_legend_param = list(title = "Read Hits", at = c(min(readMappingMatrix, na.rm = TRUE), max(readMappingMatrix, na.rm = TRUE)), labels = c("Low", "High")),
                  bottom_annotation = bottom_annotation,  # Add your bottom_annotation here
                  right_annotation = row_annotation  # Add the row annotation on the right
)


# Saving image as PNG
{
  png(file=paste(Input_directory_address,paste(Virus_name,"_Heatmap_virus_screening_max_", Sys.Date(), ".png", sep = ""),sep="/"), width=wdth_png,height=hght_png,units="px",res=300)
  draw(ht_max)
  # Decorate the heatmap body with horizontal lines between rows
  ht_lines <-  decorate_heatmap_body("Read Density", {
    for (i in 1:(nrow(above_10_coverage_matrix) - 1)) {
      grid.lines(x = c(0, 1), y = rep(1 - i / nrow(above_10_coverage_matrix), 2), gp = gpar(lty = 1, lwd = 0.2))
    }
  })
  dev.off()
} 

# Saving image as SVG
{
  svg(file=paste(Input_directory_address,paste(Virus_name,"_Heatmap_virus_screening_max_", Sys.Date(), ".svg", sep = ""),sep="/"), width=wdth_svg,height=hght_svg)
  draw(ht_max)
  # Decorate the heatmap body with horizontal lines between rows
  ht_lines <-  decorate_heatmap_body("Read Density", {
    for (i in 1:(nrow(above_10_coverage_matrix) - 1)) {
      grid.lines(x = c(0, 1), y = rep(1 - i / nrow(above_10_coverage_matrix), 2), gp = gpar(lty = 1, lwd = 0.2))
    }
  })
  dev.off()
} 
}
#-------------------------------------------------------------------------------------------


# Save session information
library("devtools")
{ 
  # Tell R to print into textfile and name of text file.
  sink(paste("session_info_for_Heatmap_virus_screening_ver6_", Sys.Date(), ".txt", sep = ""))  
  # generate session info
  sf <- session_info()
  print(sf)
  sink()  #stop sinking, =sink(NULL)
}
