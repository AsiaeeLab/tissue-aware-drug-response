#' =============================================================================
#' Drug Response Prediction - Path Configuration
#' =============================================================================
#'
#' This script provides centralized path configuration using the here package
#' for reliable project root detection.
#'
#' Loads paths from ~/Paths/drugresponseprediction.json
#'
#' Usage:
#'   library(here)
#'   source(here("methods", "00_paths.R"))
#'
#' =============================================================================

# Load required packages
library(rjson)

# Get home directory
home <- Sys.getenv("HOME", unset = NA)
if (is.na(home)) stop("Cannot find 'HOME' from environment variables.")

# Detect project root using here package (if available)
project_root <- tryCatch(
  {
    if (requireNamespace("here", quietly = TRUE)) {
      here::here()
    } else {
      getwd()
    }
  },
  error = function(e) getwd()
)

# Project name for JSON lookup
project_name <- "drugresponseprediction"

# Search for config file in standard locations
candidate_configs <- c(
  file.path(home, "Paths", paste0(project_name, ".json")),
  file.path(home, "Documents", "Paths", paste0(project_name, ".json"))
)

jinfo <- candidate_configs[file.exists(candidate_configs)][1]
if (is.na(jinfo) || !nzchar(jinfo)) {
  stop(
    "Cannot locate paths JSON for this project.\n",
    "Looked for:\n  - ",
    paste(candidate_configs, collapse = "\n  - "),
    "\n\nCreate one with paths: raw, clean, scratch"
  )
}

# Load paths from JSON
temp <- fromJSON(file = jinfo)
paths <- temp$paths

# Validate required keys
required_keys <- c("raw", "clean", "scratch")
missing_keys <- setdiff(required_keys, names(paths))
if (length(missing_keys) > 0) {
  stop(
    "Paths JSON is missing required keys: ",
    paste(missing_keys, collapse = ", "),
    "\nFile: ", jinfo
  )
}

# Add results path (inside repo) if not specified
if (is.null(paths$results)) {
  paths$results <- file.path(project_root, "results")
}

# Print confirmation
message("=== Drug Response Prediction ===")
message("Config loaded from: ", jinfo)
message("  raw:     ", paths$raw)
message("  clean:   ", paths$clean)
message("  scratch: ", paths$scratch)
message("  results: ", paths$results)

# Clean up temporary variables
rm(candidate_configs, home, jinfo, missing_keys, project_name, project_root,
   required_keys, temp)
