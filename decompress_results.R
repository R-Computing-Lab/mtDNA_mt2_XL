# Decompress all .tar.gz files across all Result* directories


repo_root   <- getwd()
result_dirs <- list.dirs(repo_root, recursive = FALSE, full.names = TRUE)
result_dirs <- result_dirs[grepl("/Result", basename(result_dirs))]

if (length(result_dirs) == 0) {
  stop("No Result* directories found in: ", repo_root, "\n",
       "Make sure you run this script from the repo root.")
}

cat("Found", length(result_dirs), "Result* directories\n")
cat("Skipping any subfolder starting with 'Analysis'\n\n")

total_done  <- 0
total_error <- 0

for (result_dir in result_dirs) {
  cat("==>", basename(result_dir), "\n")

  archives <- list.files(result_dir, pattern = "\.tar\.gz$", full.names = TRUE)

  if (length(archives) == 0) {
    cat("  No archives found, skipping\n\n")
    next
  }

  for (archive in archives) {
    cat("  Decompressing:", basename(archive), "\n")

    result <- tryCatch(
      { untar(archive, exdir = result_dir); "ok" },
      error = function(e) paste("ERROR:", e$message)
    )

    if (result == "ok") {
      cat("    Done\n")
      total_done <- total_done + 1
    } else {
      cat("   ", result, "\n")
      total_error <- total_error + 1
    }
  }
  cat("\n")
}

cat("Finished. Decompressed:", total_done, " Errors:", total_error, "\n")
