# Initialize the webR environment
webr::shim_install()

tryCatch({
  # Install and load packages with specified repositories
  install.packages('ade4', repos = c('https://adeverse.r-universe.dev', 'https://repo.r-wasm.org'))
  library(ade4)
  
  install.packages("rlang", repos = c("https://r-lib.r-universe.dev", "https://repo.r-wasm.org"))
  library(rlang)
  
  install.packages("plyr", repos = c("https://hadley.r-universe.dev", "https://repo.r-wasm.org"))
  library(plyr)
  
  webr::install("curl", repos = "https://timelyportfolio.github.io/webr_repo/")
  
  install.packages("BiocManager", repos = c("https://bioconductor.r-universe.dev", "https://repo.r-wasm.org"))
  library(BiocManager)
  BiocManager::install("rhdf5")
  BiocManager::install("zlibbioc")
  BiocManager::install("GenomeInfoDbData")

  install.packages("phyloseq", repos = c("https://bioc.r-universe.dev", "https://repo.r-wasm.org"))
  library(phyloseq)
  
  install.packages('GUniFrac', repos = 'https://repo.r-wasm.org')
  library(GUniFrac)

  install.packages('cluster', repos = c('https://mmaechler.r-universe.dev', 'https://repo.r-wasm.org'))
  library(cluster)

  install.packages('vegan', repos = c('https://vegandevs.r-universe.dev', 'https://repo.r-wasm.org'))
  library(vegan)

  # install.packages('clusterSim', repos = c('https://a-dudek-ue.r-universe.dev', 'https://repo.r-wasm.org'))
  # library(clusterSim)

  install.packages('phangorn', repos = c('https://isbool.r-universe.dev', 'https://repo.r-wasm.org'))
  library(phangorn)
  
  install.packages('mclust', repos = c('https://luca-scr.r-universe.dev', 'https://repo.r-wasm.org'))
  library(mclust)
  
  install.packages('fpc', repos = c('https://chrhennig.r-universe.dev', 'https://repo.r-wasm.org'))
  library(fpc)

}, error = function(e) {
  # Catch and print any errors that occur during installation or loading
  cat("An error occurred:", e$message, "\n")
})
