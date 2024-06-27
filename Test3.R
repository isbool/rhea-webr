library(plotrix)
library(PerformanceAnalytics)
library(reshape)
library(ggplot2)
library(gridExtra)
library(gtable)
library(Matrix)
library(cowplot)
library(ggpubr)
#####################################################################################################################
####                                        Functions to be used  in main Script.                            ########
#####################################################################################################################

# Function to calculate relative abundance
rel.abundance <- function(data)
{
  total = sum(data)
  rel.data <- 100 * data / total
  return(rel.data)
}

# Replace abundance with zero if value is below given cutoff
abundance.fix <- function(data)
{
  data[data < abundance_cutoff] <- 0
  return(data)
}

# Replace zero value with NA 
fill_zero.NA <- function(data,ReplaceZero)
{
  if (ReplaceZero == "NO") {
    return(data)
  } else if (ReplaceZero == "YES") {
    data[data == 0] <- NA
    return(data)
  } else {
    return(data)
  }
}

# Return maxima and minima values of the given input value
max.fun <- function(data)
{
  data.max <- max(as.numeric(as.character(data)), na.rm = TRUE)
  return (data.max)
}

# Calculate the prevalence of the given input table
pre.fun.na <- function(data)
{
  prevalence <- nnzero(data, na.counted = FALSE)
  return(prevalence)
}

# Return the maximum median value for each group
max.med <- function(data)
{
   max.median <-max(aggregate (data ~ independent_variable,FUN = median,simplify = TRUE)[,2])
  return(max.median)
}

# Return the prevalence ratio
max.pre <- function(data)
{
 
  # Return the number of samples (excluding NA) for each group separately
  found <- aggregate (data ~ independent_variable,FUN = pre.fun.na,simplify = TRUE,na.action = na.pass)[,2]
  
  # Return the total number of samples (including missing values) for each group separately
  all <- aggregate (data ~ independent_variable,FUN = length,simplify = TRUE,na.action = na.pass)[,2]
  
  # Calculate the ratio for each group and return the maximum ratio out of the groups
  max.ratio <- max(found / all)
  return(max.ratio)
}

# Set the theme to change text for plotting (ggplot - Gtable)
mytheme <- gridExtra::ttheme_default(
  
  # Adjust settings for the text inside table
  core = list(fg_params = list(cex = 0.8)),
  
  # Adjust the test for column and row header
  colhead = list(fg_params = list(cex = 0.9)),
  rowhead = list(fg_params = list(cex = 1.0))
)

###################       Read all required input files      ####################

# Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
original_table <- read.table (file = input_filename, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")

#####################################################################################################################
####                                      Pre-processing of OTUs Table                                       ########
#####################################################################################################################

# Convert independent variable into factor to avoid errors
original_table[,independant_variable_name] <- as.factor(original_table[,independant_variable_name]) 

# Store independent variable columns from original table 
independent_variable <- original_table[[independant_variable_name]]
ifelse(group_order != "", independent_variable <- factor(independent_variable,levels=group_order), independent_variable)

# Store metadata variable columns from original table 
my_meta_data <- original_table[1:taxonomic_variables_start - 1]

# Store relative abundance values of all OTUs
my_otu_data <- original_table[taxonomic_variables_start:dim(original_table)[2]]

# Transform data by zeroing very low abundances (based on given abundance cutoff - abundance_cutoff)
my_otu_mod =  as.data.frame(apply(my_otu_data,2,abundance.fix))

# Transform data by replacing all zero values with missing values
# Column consisting of "NA" or "0" only are removed (see below)
my_otu_mod_noz = as.data.frame(apply(my_otu_mod,2,fill_zero.NA,ReplaceZero))

# Remove column if entire OTU column contain zeros or missing values
my_otu_mod_noz <- my_otu_mod_noz[,!apply(my_otu_mod_noz , 2 , function(x) all(is.na(x) | (x == 0))), drop = FALSE]
                                                       
# Transform data by removing any OTU with median relative abundance below cutoff
t_otu_mod_noz = as.data.frame(t(my_otu_mod_noz))

# Calculate median for each OTUs
t_otu_mod_noz$max.median <- apply(t_otu_mod_noz,1,max.med)

# Select OTUs above median cutoff (med.cutoff)
selected_max <- t_otu_mod_noz[t_otu_mod_noz$max.median > max_median_cutoff,]

# Remove calculated median column "max.median"
selected_max$max.median <- NULL

# Make a separate object as data frame (columns are OTUs and rows are samples)
otu_mod_noz_max <- as.data.frame(t(selected_max))

# Transpose data (columns are samples and rows are OTUs)
t_otu_mod_noz_max = as.data.frame(t(otu_mod_noz_max))

# Transform data by removing all OTUs with prevalence below the given cutoff
t_otu_mod_noz_max$pre <- apply(t_otu_mod_noz_max,1,max.pre)
selected_pre <- t_otu_mod_noz_max[t_otu_mod_noz_max$pre > prevalence_cutoff,]
selected_pre$pre <- NULL

# Transform and filter OTU table
otu_mod_noz_max_pre <- as.data.frame(t(selected_pre))

# Merge the metadata and OTU data in one table 
# This table will be used as input for the analysis
input_table <- cbind(my_meta_data,otu_mod_noz_max_pre) 


