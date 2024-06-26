##################################################################################
######                             Main Script                              ######
##################################################################################

###################       Load all required libraries     ########################
  install.packages('phangorn', repos = c('https://isbool.r-universe.dev', 'https://repo.r-wasm.org'))

  library(ade4)
  library(vegan)
  library(GUniFrac)
  library(ape)
  library(statmod)
  library(matrixStats)
  library(inline)
  library(foreach)
  library(ggplot2)
  library(fBasics)
  library(statip)
  library(cluster)
  library(clusterSim)
  library(phangorn)

# Define a list of packages to install and their respective repositories
tryCatch({
  #webr::shim_install()
  #install.packages('ade4', repos = c('https://adeverse.r-universe.dev', 'https://repo.r-wasm.org'))
  #install.packages('vegan', repos = c('https://vegandevs.r-universe.dev', 'https://repo.r-wasm.org'))

  # To install GUniFrac (modfied) requies mounted image with the modifed GUniFrac & modifed modeest + ggrepel(?)
  # install.packages('GUniFrac', repos = 'https://repo.r-wasm.org')
  # webr::install("ape")
  # webr::install("statmod")
  # webr::install("matrixStats")
  # webr::install("inline")
  # webr::install("foreach")
  # webr::install("ggplot2")
  # webr::install("fBasics")
  # webr::install("statip")
  
  #install.packages('cluster', repos = c('https://mmaechler.r-universe.dev', 'https://repo.r-wasm.org'))
  #install.packages('clusterSim', repos = c('https://a-dudek-ue.r-universe.dev', 'https://repo.r-wasm.org'))
  

  # library(ade4)
  # library(vegan)

  # library(GUniFrac)
  # library(ape)
  # library(statmod)
  # library(matrixStats)
  # library(inline)
  # library(foreach)
  # library(ggplot2)
  # library(fBasics)
  # library(statip)

  # library(cluster)
  # library(clusterSim)
  # library(phangorn)

}, error = function(e) {
  # Catch and print any errors that occur during installation or loading
  cat("An error occurred:", e$message, "\n")
})


cluster.stats <- function (d = NULL, clustering, alt.clustering = NULL,
                           noisecluster=FALSE,
                              silhouette = TRUE, G2 = FALSE, G3 = FALSE,
                              wgap=TRUE, sepindex=TRUE, sepprob=0.1,
                              sepwithnoise=TRUE,
                              compareonly = FALSE,
                              aggregateonly = FALSE)
{
    if (!is.null(d)) 
        d <- as.dist(d)
    cn <- max(clustering)
    clusteringf <- as.factor(clustering)
    clusteringl <- levels(clusteringf)
    cnn <- length(clusteringl)
    if (cn != cnn) {
        warning("clustering renumbered because maximum != number of clusters")
        for (i in 1:cnn) clustering[clusteringf == clusteringl[i]] <- i
        cn <- cnn
    }
    n <- length(clustering)
    noisen <- 0
    cwn <- cn
    if (noisecluster){
      noisen <- sum(clustering==cn)
      cwn <- cn-1
    }
    diameter <- average.distance <- median.distance <- separation <- average.toother <- cluster.size <- within.dist <- between.dist <- numeric(0)
    for (i in 1:cn) cluster.size[i] <- sum(clustering == i)
    pk1 <- cluster.size/n
    pk10 <- pk1[pk1 > 0]
    h1 <- -sum(pk10 * log(pk10))
    corrected.rand <- vi <- NULL
    if (!is.null(alt.clustering)) {
        choose2 <- function(v) {
            out <- numeric(0)
            for (i in 1:length(v)) out[i] <- ifelse(v[i] >= 2, 
                choose(v[i], 2), 0)
            out
        }
        cn2 <- max(alt.clustering)
        clusteringf <- as.factor(alt.clustering)
        clusteringl <- levels(clusteringf)
        cnn2 <- length(clusteringl)
        if (cn2 != cnn2) {
            warning("alt.clustering renumbered because maximum != number of clusters")
            for (i in 1:cnn2) alt.clustering[clusteringf == clusteringl[i]] <- i
            cn2 <- cnn2
        }
        nij <- table(clustering, alt.clustering)
        dsum <- sum(choose2(nij))
        cs2 <- numeric(0)
        for (i in 1:cn2) cs2[i] <- sum(alt.clustering == i)
        sum1 <- sum(choose2(cluster.size))
        sum2 <- sum(choose2(cs2))
        pk2 <- cs2/n
        pk12 <- nij/n
        corrected.rand <- (dsum - sum1 * sum2/choose2(n))/((sum1 + 
            sum2)/2 - sum1 * sum2/choose2(n))
        pk20 <- pk2[pk2 > 0]
        h2 <- -sum(pk20 * log(pk20))
        icc <- 0
        for (i in 1:cn) for (j in 1:cn2) if (pk12[i, j] > 0) 
            icc <- icc + pk12[i, j] * log(pk12[i, j]/(pk1[i] * 
                pk2[j]))
        vi <- h1 + h2 - 2 * icc
    }
    if (compareonly) {
        out <- list(corrected.rand = corrected.rand, vi = vi)
    }
    else {
#        if (silhouette) 
#            require(cluster)
        dmat <- as.matrix(d)
        within.cluster.ss <- 0
        overall.ss <- nonnoise.ss <- sum(d^2)/n
        if (noisecluster)
          nonnoise.ss <- sum(as.dist(dmat[clustering<=cwn,
                                          clustering<=cwn])^2)/
                                            sum(clustering<=cwn)
        ave.between.matrix <-
          separation.matrix <- matrix(0, ncol = cn, nrow = cn)
        for (i in 1:cn) {
            cluster.size[i] <- sum(clustering == i)
            di <- as.dist(dmat[clustering == i, clustering == 
                i])
            if (i<=cwn){
              within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size[i]
              within.dist <- c(within.dist, di)
            }
            if (length(di) > 0) 
                diameter[i] <- max(di)
            else diameter[i] <- NA
            average.distance[i] <- mean(di)
            median.distance[i] <- median(di)
            bv <- numeric(0)
            for (j in 1:cn) {
                if (j != i) {
                  sij <- dmat[clustering == i, clustering == 
                    j]
                  bv <- c(bv, sij)
                  if (i < j) {
                    separation.matrix[i, j] <- separation.matrix[j, 
                      i] <- min(sij)
                    ave.between.matrix[i, j] <- ave.between.matrix[j, i] <-
                        mean(sij)
                    if (i<=cwn & j<=cwn)
                      between.dist <- c(between.dist, sij)
                  }
                }
            }
            separation[i] <- min(bv)
            average.toother[i] <- mean(bv)
        }
        average.between <- mean(between.dist)
#        average.within <- mean(within.dist)
        average.within <- weighted.mean(average.distance,cluster.size,na.rm=TRUE)
        nwithin <- length(within.dist)
        nbetween <- length(between.dist)
        between.cluster.ss <- nonnoise.ss - within.cluster.ss
        ch <- between.cluster.ss * (n - noisen - cwn)/(within.cluster.ss * 
            (cwn - 1))
        clus.avg.widths <- avg.width <- NULL
        if (silhouette) {
            sii <- silhouette(clustering, dmatrix = dmat)
            sc <- summary(sii)
            clus.avg.widths <- sc$clus.avg.widths
            if (noisecluster)
              avg.width <- mean(sii[clustering<=cwn,3])
            else
              avg.width <- sc$avg.width
        }
        g2 <- g3 <- cn2 <- cwidegap <- widestgap <- sindex <- NULL
        if (G2) {
            splus <- sminus <- 0
            for (i in 1:nwithin) {
                splus <- splus + sum(within.dist[i] < between.dist)
                sminus <- sminus + sum(within.dist[i] > between.dist)
            }
            g2 <- (splus - sminus)/(splus + sminus)
        }
        if (G3) {
            sdist <- sort(c(within.dist, between.dist))
            sr <- nwithin + nbetween
            dmin <- sum(sdist[1:nwithin])
            dmax <- sum(sdist[(sr - nwithin + 1):sr])
            g3 <- (sum(within.dist) - dmin)/(dmax - dmin)
        }
        pearsongamma <- cor(c(within.dist, between.dist), c(rep(0, 
            nwithin), rep(1, nbetween)))
        dunn <- min(separation[1:cwn])/max(diameter[1:cwn],na.rm=TRUE)
        acwn <- ave.between.matrix[1:cwn,1:cwn]
        dunn2 <- min(acwn[upper.tri(acwn)])/
          max(average.distance[1:cwn],na.rm=TRUE)
        if (wgap){
          cwidegap <- rep(0,cwn)
          for (i in 1:cwn)
            if (sum(clustering==i)>1)
              cwidegap[i] <- max(hclust(as.dist(dmat[clustering==i,
                                                   clustering==i]),
                                      method="single")$height)
          widestgap <- max(cwidegap)
        }
        if (sepindex){
          psep <- rep(NA,n)
          if (sepwithnoise | !noisecluster){
            for (i in 1:n)
              psep[i] <- min(dmat[i,clustering!=clustering[i]])
            minsep <- floor(n*sepprob)
          }
          else{
            dmatnn <- dmat[clustering<=cwn,clustering<=cwn]
            clusteringnn <- clustering[clustering<=cwn]
            for (i in 1:(n-noisen))
              psep[i] <- min(dmatnn[i,clusteringnn!=clusteringnn[i]])
            minsep <- floor((n-noisen)*sepprob)
          }
          sindex <- mean(sort(psep)[1:minsep])
        }
        if (!aggregateonly)
          out <- list(n = n, cluster.number = cn, cluster.size = cluster.size,
                      min.cluster.size = min(cluster.size[1:cwn]),
                      noisen=noisen,
            diameter = diameter, average.distance = average.distance, 
            median.distance = median.distance, separation = separation, 
            average.toother = average.toother, separation.matrix = separation.matrix,
                      ave.between.matrix=ave.between.matrix,
            average.between = average.between, average.within = average.within, 
            n.between = nbetween, n.within = nwithin,
                      max.diameter = max(diameter[1:cwn],na.rm=TRUE),
                      min.separation = sepwithnoise*min(separation)+
                      (!sepwithnoise)*min(separation[1:cwn]),
                      within.cluster.ss = within.cluster.ss, 
            clus.avg.silwidths = clus.avg.widths, avg.silwidth = avg.width, 
            g2 = g2, g3 = g3, pearsongamma = pearsongamma, dunn = dunn,
                      dunn2=dunn2,
            entropy = h1, wb.ratio = average.within/average.between, 
            ch = ch, cwidegap=cwidegap, widestgap=widestgap,
                    sindex=sindex,
                    corrected.rand = corrected.rand, vi = vi)
        else
          out <- list(n = n, cluster.number = cn,
                      min.cluster.size = min(cluster.size[1:cwn]),
                      noisen=noisen,
            average.between = average.between, average.within = average.within, 
                      max.diameter = max(diameter[1:cwn],na.rm=TRUE),
                      min.separation = sepwithnoise*min(separation)+
                      (!sepwithnoise)*min(separation[1:cwn]), 
            ave.within.cluster.ss = within.cluster.ss/(n-noisen), 
            avg.silwidth = avg.width, 
            g2 = g2, g3 = g3, pearsongamma = pearsongamma, dunn = dunn,
                      dunn2=dunn2,
            entropy = h1, wb.ratio = average.within/average.between, 
            ch = ch, widestgap=widestgap,
                    sindex=sindex,
                    corrected.rand = corrected.rand, vi = vi)
    }
    out
}

###################       Read all required input files      ####################
# Load the mapping file containing individual sample information (sample names in the first column)
meta_file <- read.table(file = input_meta, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")

# Clean table from empty lines
meta_file <- data.frame(meta_file[!apply(is.na(meta_file) | meta_file == "", 1, all),, drop = FALSE])

# Order the mapping file by sample names (ascending)
meta_file <- data.frame(meta_file[order(row.names(meta_file)),, drop = FALSE])

# Save the position of the target group name in the mapping file
meta_file_pos <- which(colnames(meta_file) == group_name)

# Select metadata group based on the pre-set group name
all_groups <- as.factor(meta_file[, meta_file_pos])

#------------------------------------------------------------------------

# Load the tab-delimited file containing the values to be analyzed (samples names in the first column)
otu_file <- read.table(file = input_otu, check.names = FALSE, header = TRUE, dec = ".", sep = "\t", row.names = 1, comment.char = "")

# Clean table from empty lines
otu_file <- otu_file[!apply(is.na(otu_file) | otu_file == "", 1, all),]

# keep only those rows that appear in the mapping file
otu_file <- otu_file[, rownames(meta_file)]

# OTU-table and mapping file should have the same order and number of sample names
# Order the OTU-table by sample names (ascending)
otu_file <- otu_file[, order(names(otu_file))]

# Transpose OTU-table and convert format to a data frame
otu_file <- data.frame(t(otu_file), check.names = FALSE)

#------------------------------------------------------------------------

# Load the phylogenetic tree calculated from the OTU sequences 
tree_file <- read.tree(input_tree)

# Remove single quotes from the tips of the tree
tree_file$tip.label <- gsub("'", "", tree_file$tip.label)

# Root the OTU tree at midpoint 
rooted_tree <- midpoint(tree_file)


####################       Calculate beta-diversity          ###################

# Create the directory where all output files are saved (is named after the target group name set above for comparisons)
dir.create(group_name)

# Calculate the UniFrac distance matrix for comparing microbial communities
unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0, 0.5, 1.0))$unifracs

# Weight on abundant lineages so the distance is not dominated by highly abundant lineages with 0.5 having the best power
unifract_dist <- unifracs[,, "d_0.5"]

################ Generate tree #######################

# Save the UniFrac output as distance object
all_dist_matrix <- as.dist(unifract_dist)

# Apply a hierarchical cluster analysis on the distance matrix based on the Ward's method
all_fit <- hclust(all_dist_matrix, method = "ward.D2")

# Generates a tree from the hierarchically generated object
tree <- as.phylo(all_fit)
my_tree_file_name <- paste(group_name, "/phylogram.pdf", sep = "")
plot_color <- rainbow(length(levels(all_groups)))[all_groups]

# Save the generated phylogram in a pdf file
pdf(my_tree_file_name)

# The tree is visualized as a Phylogram color-coded by the selected group name
plot(tree, type = "phylogram", use.edge.length = TRUE, tip.color = (plot_color), label.offset = 0.01)
print.phylo(tree)
axisPhylo()
tiplabels(pch = 16, col = plot_color)
dev.off()

#################            Build NMDS plot           ########################

# Generated figures are saved in a pdf file 
file_name <- paste(group_name, "beta-diversity.pdf", sep = "_")
pdf(paste(group_name, "/", file_name, sep = ""))

# Calculate the significance of variance to compare multivariate sample means (including two or more dependent variables)
# Omit cases where there isn't data for the sample (NA)
all_groups_comp <- all_groups[!is.na(all_groups)]
unifract_dist_comp <- unifract_dist[!is.na(all_groups), !is.na(all_groups)]
adonis <- adonis2(as.dist(unifract_dist_comp) ~ all_groups_comp)
permdisp <- permutest(betadisper(as.dist(unifract_dist_comp), as.factor(all_groups_comp), type = "median"))
all_groups_comp <- factor(all_groups_comp, levels(all_groups_comp)[unique(all_groups_comp)])

if (nrow(unifract_dist_comp) > 2) {
  # Calculate and display the MDS plot (Multidimensional Scaling plot)
  s.class(
    cmdscale(unifract_dist_comp, k = 2), col = unique(plot_color), cpoint =
      2, fac = all_groups_comp, sub = paste("MDS plot of Microbial Profiles\nPERMDISP     p=", permdisp[["tab"]][["Pr(>F)"]][1], "\n",
                                            "PERMANOVA  p=", adonis[1, 5], sep = "")
  )
  if (label_samples == 1) {
    lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
    ifelse(label_id != "", lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
    text(cmdscale(unifract_dist_comp, k = 2), labels = lab_samples, cex = 0.7, adj = c(-.1, -.8))
  }

  # Calculate and display the NMDS plot (Non-metric Multidimensional Scaling plot)
  meta <- metaMDS(unifract_dist_comp, k = 2)
  s.class(
    meta$points, col = unique(plot_color), cpoint = 2, fac = all_groups_comp,
    sub = paste("metaNMDS plot of Microbial Profiles\nPERMDISP     p=", permdisp[["tab"]][["Pr(>F)"]][1], "\n",
                "PERMANOVA  p=", adonis[1, 5], sep = "")
  )
  if (label_samples == 1) {
    lab_samples <- row.names(meta$points)
    ifelse(label_id != "", lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), ""), lab_samples)
    text(meta$points, labels = lab_samples, cex = 0.7, adj = c(-.1, -.8))
  }
}

#close the pdf file
dev.off()

###############          NMDS for pairwise analysis        ###################

# This plot is only generated if there are more than two groups included in the comparison
# Calculate the pairwise significance of variance for group pairs
# Get all groups contained in the mapping file
unique_groups <- levels(all_groups_comp)
if (dim(table(unique_groups)) > 2) {

  # Initialise vector and lists
  pVal = NULL
  permdisppval = NULL
  pairedMatrixList <- list(NULL)
  pair_1_list <- NULL
  pair_2_list <- NULL

  for (i in 1:length(combn(unique_groups, 2)[1, ])) {

    # Combine all possible pairs of groups
    pair_1 <- combn(unique_groups, 2)[1, i]
    pair_2 <- combn(unique_groups, 2)[2, i]

    # Save pairs information in a vector
    pair_1_list[i] <- pair_1
    pair_2_list[i] <- pair_2

    # Generate a subset of all samples within the mapping file related to one of the two groups
    inc_groups <-
      rownames(subset(meta_file, meta_file[, meta_file_pos] == pair_1
                      |
                        meta_file[, meta_file_pos] == pair_2))

    # Convert UniFrac distance matrix to data frame
    paired_dist <- as.data.frame(unifract_dist_comp)

    # Save all row names of the mapping file
    row_names <- rownames(paired_dist)

    # Add row names to the distance matrix
    paired_dist <- cbind(row_names, paired_dist)

    # Generate distance matrix with samples of the compared groups (column-wise)
    paired_dist <- paired_dist[sapply(paired_dist[, 1], function(x) all(x %in% inc_groups)),]

    # Remove first column with unnecessary group information
    paired_dist[, 1] <- NULL
    paired_dist <- rbind(row_names, paired_dist)

    # Generate distance matrix with samples of the compared group (row-wise)
    paired_dist <- paired_dist[, sapply(paired_dist[1,], function(x) all(x %in% inc_groups))]

    # Remove first row with unnecessary group information 
    paired_dist <- paired_dist[-1,]

    # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
    paired_matrix <- as.matrix(paired_dist)
    class(paired_matrix) <- "numeric"

    # Save paired matrix in list
    pairedMatrixList[[i]] <- paired_matrix

    # Applies multivariate analysis to a pair out of the selected groups
    adonis <- adonis2(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                        all_groups_comp == pair_2])

    permdisp <- permutest(betadisper(as.dist(paired_matrix), as.factor(all_groups_comp[all_groups_comp == pair_1 |
                                                                                        all_groups_comp == pair_2]), type = "median"), pairwise = T)

    # List p-values
    pVal[i] <- adonis[1, 5]
    permdisppval[i] <- permdisp$pairwise[2]

  }

  # Adjust p-values for multiple testing according to Benjamini-Hochberg method
  pVal_BH <- round(p.adjust(pVal, method = "BH", n = length(pVal)), 4)
  permdisppval_BH <- round(p.adjust(permdisppval, method = "BH", n = length(permdisppval)), 4)


  # Generated NMDS plots are stored in one pdf file called "pairwise-beta-diversity-nMDS.pdf"
  file_name <- paste(group_name, "pairwise-beta-diversity-NMDS.pdf", sep = "_")
  pdf(paste(group_name, "/", file_name, sep = ""))

  for (i in 1:length(combn(unique_groups, 2)[1, ])) {
    if (nrow(pairedMatrixList[[i]]) > 2) {
      meta <- metaMDS(pairedMatrixList[[i]], k = 2)
      s.class(
        meta$points,
        col = rainbow(length(levels(all_groups_comp))), cpoint = 2,
        fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                          all_groups_comp == pair_2_list[i]]),
        sub = paste("NMDS plot of Microbial Profiles\n ", pair_1_list[i], " - ", pair_2_list[i], "\n PERMDISP     p=", permdisppval[[i]], ",", "  p.adj=", permdisppval_BH[i], "\n",
                    " PERMANOVA  p=", pVal[i], ",", " p.adj=", pVal_BH[i], sep = "")
      )
    }
  }
  dev.off()

  # Generated MDS plots are stored in one pdf file called "pairwise-beta-diversity-MDS.pdf"
  file_name <- paste(group_name, "pairwise-beta-diversity-MDS.pdf", sep = "_")
  pdf(paste(group_name, "/", file_name, sep = ""))

  for (i in 1:length(combn(unique_groups, 2)[1, ])) {
    if (nrow(pairedMatrixList[[i]]) > 2) {
      # Calculate and display the MDS plot (Multidimensional Scaling plot)
      s.class(
        cmdscale(pairedMatrixList[[i]], k = 2), col = rainbow(length(levels(all_groups_comp))), cpoint =
          2, fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                               all_groups_comp == pair_2_list[i]]),
        sub = paste("MDS plot of Microbial Profiles\n ", pair_1_list[i], " - ", pair_2_list[i], "\n PERMDISP     p=", permdisppval[[i]], ",", "  p.adj=", permdisppval_BH[i], "\n",
                    " PERMANOVA  p=", pVal[i], ",", " p.adj=", pVal_BH[i], sep = "")
      )
    }
  }
  dev.off()

}

######                        Determine number of clusters                           ######

ch_nclusters = NULL
sil_nclusters = NULL
dunn_nclusters = NULL
db_nclusters = NULL

if (dim(otu_file)[1] - 1 <= kmers_limit) {
  kmers_limit = dim(otu_file)[1] - 1
}
for (k in 1:kmers_limit) {
  if (k == 1) {
    ch_nclusters[k] = NA
    sil_nclusters[k] = NA
    dunn_nclusters = NA
    db_nclusters = NA
  } else {
    # Partitioning the data into k clusters (max k is number of samples within the dataset)
    data_cluster = as.vector(pam(as.dist(unifract_dist_comp), k, diss = TRUE)$clustering)

    # Calculate Calinski-Harabasz and silhouette Index 
    index = cluster.stats(as.dist(unifract_dist_comp), data_cluster)
    index_db = index.DB(x = otu_file, cl = data_cluster, d = as.dist(unifract_dist_comp), centrotypes = "medoids")
    ch_nclusters[k] <- index[["ch"]]
    sil_nclusters[k] <- index[["avg.silwidth"]]
    dunn_nclusters[k] <- index[["dunn2"]]
    db_nclusters[k] <- index_db[["DB"]]
    print(k)
  }
}

# Generated plot showing the optimal number of clusters
for (i in 1:2) {
  if (i == 1) { pdf("De-novo-clustering.pdf") }
  if (i == 2) { pdf(paste0(group_name, "/", "De-novo-clustering.pdf")) }

  plot(ch_nclusters, type = "h", xlab = "k clusters", ylab = "CH index", main = "Optimal number of clusters (CH index)")
  title(sub = "*The higher the value the better", adj = 0, cex.sub = 0.9)
  plot(sil_nclusters, type = "h", xlab = "k clusters", ylab = "Average silhouette width", main = "Optimal number of clusters (Silhouette index)")
  title(sub = "*The higher the value the better", adj = 0, cex.sub = 0.9)
  plot(dunn_nclusters, type = "h", xlab = "k clusters", ylab = "Dunn Index", main = "Optimal number of clusters (Dunn Index)")
  title(sub = "*The higher the value the better", adj = 0, cex.sub = 0.9)
  plot(db_nclusters, type = "h", xlab = "k clusters", ylab = "Davies-Bouldin Index", main = "Optimal number of clusters (Davies-Bouldin Index)")
  title(sub = "*The lower the value the better", adj = 0, cex.sub = 0.9)

  dev.off()
}

#################################################################################
######                        Write Output Files                           ######
#################################################################################

# Write the distance matrix table in a file
file_name <- paste(group_name, "distance-matrix-gunif.tab", sep = "_")
write.table(unifract_dist_comp, paste(group_name, "/", file_name, sep = ""), sep = "\t", col.names = NA, quote = FALSE)
write.table(unifract_dist_comp, "distance-matrix-gunif.tab", sep = "\t", col.names = NA, quote = FALSE)
write.tree(tree, "samples-Tree.nwk", tree.names = FALSE)



#################################################################################
######                           End of Script                             ######
#################################################################################
