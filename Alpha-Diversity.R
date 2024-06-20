##################################################################################
######                        Diversity Functions                           ###### 
##################################################################################

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count=sum(x[x>0.5]^0)
  return(count)
}

# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x)
{
  # Count only the OTUs that are present more than the set proportion
  total=sum(x)
  count=sum(x[x/total>eff.cutoff]^0)
  return(count)
}

# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  # Given a fixed Normalization reads depth
  total=sum(x)
  count=sum(x[norm.cutoff*x/total>0.5]^0)
  return(count)
}


# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total=sum(x)
  se=-sum(x[x>0]/total*log(x[x>0]/total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total=sum(x)
  se=round(exp(-sum(x[x>0]/total*log(x[x>0]/total))),digits =2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total=sum(x)
  si=sum((x[x>0]/total)^2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total=sum(x)
  si=round(1/sum((x[x>0]/total)^2),digits =2)
  return(si)
}

##################################################################################
######                             Main Script                              ###### 
##################################################################################

# Read a normalized OTU-table without taxonomy  
otu_table <- read.table (file_name, 
                       check.names = FALSE, 
                       header=TRUE, 
                       dec=".", 
                       sep = "\t",
                       row.names = 1)

# Clean table from empty lines
otu_table <- otu_table[!apply(is.na(otu_table) | otu_table=="",1,all),]

# Order and transpose OTU-table
my_otu_table <- otu_table[,order(names(otu_table))] 
my_otu_table <-data.frame(t(my_otu_table))

# Apply diversity functions to table
otus_div_stats<-data.frame(my_otu_table[,0])
otus_div_stats$Richness<-apply(my_otu_table,1,Species.richness)
otus_div_stats$Normalized.Richness<-apply(my_otu_table,1,Norm.Species.richness)
otus_div_stats$Effective.Richness<-apply(my_otu_table,1,Eff.Species.richness)
otus_div_stats$Shannon.Index<-apply(my_otu_table,1,Shannon.entropy)
otus_div_stats$Shannon.Effective<-apply(my_otu_table,1,Shannon.effective)
otus_div_stats$Simpson.Index<-apply(my_otu_table,1,Simpson.concentration)
otus_div_stats$Simpson.Effective<-apply(my_otu_table,1,Simpson.effective)
otus_div_stats$Evenness <- otus_div_stats$Shannon.Index/log(otus_div_stats$Richness,2)


# Write the results in a file and copy in directory "Serial-Group-Comparisons" if existing
write.table(otus_div_stats, "alpha-diversity.tab", sep="\t", col.names=NA, quote=FALSE)
suppressWarnings (try(write.table(otus_div_stats[c(1,2,3,5,7)], "../5.Serial-Group-Comparisons/alpha-diversity.tab", sep = "\t",col.names = NA, quote = FALSE), silent =TRUE))

##################################################################################
######                          End of Script                               ###### 
##################################################################################
