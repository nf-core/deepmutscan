# List of required packages
required_packages <- c("Biostrings", "dplyr", "ggplot2", "scales", "zoo", "reshape2", "tidyr", "tidyverse")

# Function to install and load missing packages
install_if_missing <- function(pkg) {
  # Check if package is already installed
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))

    # Ensure BiocManager is installed for Bioconductor packages
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")  # Set a CRAN mirror
    }

    # Install Bioconductor packages via BiocManager
    if (pkg %in% c("Biostrings", "GenomicRanges", "IRanges")) {  # Add more Bioconductor packages if needed
      BiocManager::install(pkg)
    } else {
      # Install CRAN packages
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  } else {
    message(paste("Package", pkg, "is already installed"))
  }
}

# Install and load all required packages
invisible(lapply(required_packages, install_if_missing))

message("All required packages are installed and loaded.")



















##Find current versions of dependencies:
# List of your packages
packages <- c("dplyr", "ggplot2", "scales", "zoo", "reshape2", "tidyr", "tidyverse", "Biostrings")

# Function to get the current version of each package
get_package_version <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    return(packageVersion(pkg))
  } else {
    return(NA)
  }
}

# Check and print versions
package_versions <- sapply(packages, get_package_version)
print(package_versions)



# Check the version of libcurl
system("curl --version")

# Check the version of OpenSSL
system("openssl version")

# Check the version of libxml2
system("xml2-config --version")

# Check the version of zlib
system("zlib-flate -version")

system("ldconfig -p | grep zlib")

