##################################################################################
######                             Main Script                              ######
##################################################################################

###################       Load all required libraries     ########################

# Attempt to install and load the 'vegan' package

tryCatch({
    # webr::shim_install()
    # webr::install("vegan")
    library("vegan")
}, error = function(e) {
    cat("An error occurred:", e$message, "\n")
})

###################       Read all required input files      ####################

# Load the tab-delimited file containing the values to be be checked (rownames in the first column)
otu_table <- read.table(file_name, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")
if (dim(otu_table)[1] == 0) {
    stop("OTU table is empty or not properly formatted.")
} else {
    print(sprintf("OTU table loaded with %d rows and %d columns.", nrow(otu_table), ncol(otu_table)))
}

# Making sure tha taxonomy column gets lower-case
col_names <- colnames(otu_table)
tax_ind <- which(sapply(tolower(col_names),
                        function(x) "taxonomy" %in% x, USE.NAMES = FALSE))
if (length(tax_ind) != 0) {
  col_names[tax_ind] <- tolower(col_names[tax_ind])
  colnames(otu_table) <- col_names
}
rm(col_names)
rm(tax_ind)

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

####################       Normalize OTU Table          ###################


# Save taxonomy information in vector
taxonomy <- as.vector(otu_table$taxonomy)

# Delete column with taxonomy information in dataframe
otu_table$taxonomy <- NULL



if (level == 0) {
  # Calculate the minimum sum of all columns/samples
  min_sum <- min(colSums(otu_table))
} else {
  # The minimum size is set to a fixed reference level
  min_sum <- normCutoff
}


if (method == 0) {
  # Divide each value by the sum of the sample and multiply by the minimal sample sum
  norm_otu_table <- t(min_sum * t(otu_table) / colSums(otu_table))
} else {
  # Rarefy the OTU table to an equal sequencing depth
  norm_otu_table <- Rarefy(t(otu_table),depth = min_sum)
  norm_otu_table <- t(as.data.frame(norm_otu_table$otu.tab.rff))
}

# Calculate relative abundances for all OTUs over all samples
# Divide each value by the sum of the sample and multiply by 100
rel_otu_table <- t(100 * t(otu_table) / colSums(otu_table))


# Re-insert the taxonomy information in normalized counts table
norm_otu_table_tax <- cbind(norm_otu_table,taxonomy)

# Re-insert the taxonomy information in relative abundance table
rel_otu_table_tax <- cbind(rel_otu_table,taxonomy)

# Debugging: Check the contents of the tables
print(head(norm_otu_table))
print(head(rel_otu_table))
print("Summaries of processed tables:")
print(summary(norm_otu_table))
print(summary(rel_otu_table))

################################################################################
# Generate a two-sided pdf with a rarefaction curve for all samples and a curve
pdf(file = "RarefactionCurve.pdf")

# Plot the rarefaction curve for all samples
rarefactionCurve <- rarecurve(data.frame(t(otu_table)),
                              step = 20,
                              col = "black",
                              lty = "solid",
                              label = F,
                              xlab = "Number of Reads",
                              ylab = "Number of Species",
                              main = "Rarefaction Curves of All Samples")

# Generate empty vectors for the analysis of the rarefaction curve
slope = vector()
SampleID = vector()
angle <- c()

# Iterate through all samples
for (i in seq_along(rarefactionCurve)) {
  # If the sequencing depth is greater than 100, the slope of the line that passes between the last and last-100 count is calculated
  richness <- ifelse(length(rarefactionCurve[[i]]) > 6, 
                     (rarefactionCurve[[i]][length(rarefactionCurve[[i]])] - rarefactionCurve[[i]][length(rarefactionCurve[[i]])-5])/(attr(rarefactionCurve[[i]], "Subsample")[length(rarefactionCurve[[i]])]-attr(rarefactionCurve[[i]], "Subsample")[length(rarefactionCurve[[i]])-5]) , 1000)
  angle[i] <- ifelse(richness!=1000, atan(richness)*180/pi, NA)
  slope <- c(slope,richness)
  SampleID <- c(SampleID,as.character(names(otu_table)[i]))
}

# Generate the output table for rarefaction curve
curvedf <- cbind(SampleID, slope, angle)
ordered_vector <- order(as.numeric(curvedf[,2]), decreasing = TRUE)
curvedf <- curvedf[order(as.numeric(curvedf[,2]), decreasing = TRUE),]

# Generates a graph with all samples
# Underestimated cases are shown in red
for (i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[ordered_vector[i]]], "Subsample")
  lines(N, rarefactionCurve[[ordered_vector[i]]],col="red")
}

# Determine the plotting width and height
Nmax <- sapply(rarefactionCurve[ordered_vector[1:labelCutoff]], function(x) max(attr(x, "Subsample")))
Smax <- sapply(rarefactionCurve[ordered_vector[1:labelCutoff]], max)

# Creates an empty plot for rarefaction curves of underestimated cases
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of Reads",
     ylab = "Number of Species", type = "n", main=paste(labelCutoff,"- most undersampled cases"))

for (i in 1:labelCutoff) {
  N <- attr(rarefactionCurve[[ordered_vector[i]]], "Subsample")
  lines(N, rarefactionCurve[[ordered_vector[i]]],col="red")
  text(max(attr(rarefactionCurve[[ordered_vector[i]]],"Subsample")), max(rarefactionCurve[[ordered_vector[i]]]), curvedf[i,1],cex=0.6)
}

dev.off()

#################################################################################
######                        Write Output Files                           ######
#################################################################################

# Write the normalized table in a file and copy in directories alpha-diversity and beta-diversity if existing
write.table(norm_otu_table, "OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(norm_otu_table, "../2.Alpha-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))
suppressWarnings (try(write.table(norm_otu_table, "../3.Beta-Diversity/OTUs_Table-norm.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

# Write the normalized table with taxonomy in a file
write.table(norm_otu_table_tax, "OTUs_Table-norm-tax.tab", sep = "\t",col.names = NA, quote = FALSE)

# Write the normalized relative abundance table in a file and copy in directory Serial-Group-Comparisons if existing
write.table(rel_otu_table, "OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table, "../5.Serial-Group-Comparisons/OTUs_Table-norm-rel.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

# Write the normalized relative abundance with taxonomy table in a file and copy in directory Taxonomic-Binning if existing
write.table(rel_otu_table_tax, "OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE)
suppressWarnings (try(write.table(rel_otu_table_tax, "../4.Taxonomic-Binning/OTUs_Table-norm-rel-tax.tab", sep ="\t",col.names = NA, quote = FALSE), silent =TRUE))

# Write the rarefaction table
write.table(curvedf, "RarefactionCurve.tab", sep ="\t", quote = FALSE, row.names = FALSE)

#################################################################################
######                           End of Script                             ######
#################################################################################
