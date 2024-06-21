# Initialize necessary setup for webR environment
webr::shim_install()

# Function to install and load a package
install_and_load <- function(pkg_name, pkg_repo) {
  # Attempt to install the package
  if (!require(pkg_name, character.only = TRUE)) {
    message(paste("Installing package:", pkg_name))
    install.packages(pkg_name, repos = pkg_repo)
    
    # After installation, attempt to load the package
    if (!require(pkg_name, character.only = TRUE)) {
      stop(paste("Failed to load package:", pkg_name))
    }
  } else {
    message(paste("Package already loaded:", pkg_name))
  }
}

# List of packages and their repositories
packages <- list(
  ade4 = c('https://adeverse.r-universe.dev', 'https://cloud.r-project.org'),
  GUniFrac = 'https://cloud.r-project.org',
  cluster = c('https://mmaechler.r-universe.dev', 'https://cloud.r-project.org'),
  vegan = c('https://vegandevs.r-universe.dev', 'https://cloud.r-project.org'),
  clusterSim = c('https://a-dudek-ue.r-universe.dev', 'https://cloud.r-project.org'),
  phangorn = c('https://isbool.r-universe.dev', 'https://cloud.r-project.org'),
  fpc = c('https://chrhennig.r-universe.dev', 'https://cloud.r-project.org')
)

# Loop through the packages and apply the installation/loading function
tryCatch({
  for (pkg_name in names(packages)) {
    install_and_load(pkg_name, packages[[pkg_name]])
  }
}, error = function(e) {
  cat("An error occurred:", e$message, "\n")
})
