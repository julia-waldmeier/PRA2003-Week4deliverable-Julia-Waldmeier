# PRA2003 - Week 4 (originally the Week 1 question): asymmetry as a
# function of momentum
#
# For each of the 6 WT/mutant strain pairs, checks whether the asymmetry
# between the WT/normal form and the mutant/resistant form depends on
# particle momentum (p = sqrt(px^2 + py^2 + pz^2)).
#
# Uses all 10 data files (the full ~5M-event sample) - the counting-only
# scripts (analyzeAllSubsamples.R etc.) never tracked momentum, so this
# re-reads the raw files. Processed one file at a time.
#
# Method: every occurrence of one of the 12 target strains has its
# momentum computed and is sorted into one of nBins momentum bins
# (boundaries computed from Set 1's quantiles, then reused for every file
# so all 10 files share the same bin definitions). Within each bin, the
# SUB-SAMPLING method is used (same as the rest of this report, sections
# 4/5/6) rather than a Poisson-counting formula:
#   - each of the 10 files gets its OWN asymmetry for that bin,
#     A_file = (count_WT - count_mutant) / (count_WT + count_mutant)
#   - the reported asymmetry is the MEAN of the 10 per-file values
#   - the reported uncertainty is the standard deviation of the 10
#     per-file values (raw SD, no division by sqrt(10), matching the
#     convention used throughout this report)
# A file with zero particles of that strain in that bin (count_WT +
# count_mutant == 0) is excluded from that bin's mean/SD, since its
# asymmetry would be an undefined 0/0.
#
# Writes two files:
#   momentum_asymmetry_results.csv    - the final per-strain, per-bin result
#   momentum_perfile_results.csv      - the full per-file breakdown (cached
#                                        so future changes, e.g. a different
#                                        number of bins, don't require
#                                        re-reading the raw data again)

# ---- Settings ----
nBins <- 5  # change this to re-bin with a different number of momentum bins

filenames <- paste0("output-Set", 1:10, ".txt")
for (f in filenames) {
  if (!file.exists(f)) {
    stop("Cannot find '", f, "' in the current working directory (", getwd(), ").")
  }
}
nFiles <- length(filenames)

pairPositiveID <- c(211, 321, 2212, 3122, 3312, 3334)
pairNegativeID <- c(-211, -321, -2212, -3122, -3312, -3334)
pairName <- c("E. coli", "Bacillus subtilis", "Pseudomonas aeruginosa",
              "Streptococcus pneumoniae", "Mycobacterium tuberculosis",
              "Salmonella")
targetIDs <- c(pairPositiveID, pairNegativeID)
nPairs <- length(pairPositiveID)

chunkSize <- 1000000

# ---- Reads one file, returning the momentum and ID of every occurrence
# of a target strain in it ----
extractMomentumAndID <- function(filepath) {
  con <- file(filepath, "r")
  on.exit(close(con))

  momentumList <- list()
  idList <- list()
  chunkIndex <- 0

  repeat {
    chunk <- readLines(con, n = chunkSize)
    if (length(chunk) == 0) {
      break
    }
    chunkIndex <- chunkIndex + 1

    nFields <- lengths(gregexpr(" ", chunk, fixed = TRUE, useBytes = TRUE)) + 1
    isData <- (nFields == 4)
    dataLines <- chunk[isData]

    if (length(dataLines) > 0) {
      ids <- as.integer(sub("^.* ", "", dataLines, perl = TRUE, useBytes = TRUE))
      keep <- ids %in% targetIDs

      if (any(keep)) {
        keptLines <- dataLines[keep]
        parsed <- read.table(text = keptLines, colClasses = "numeric")
        p <- sqrt(parsed[[1]]^2 + parsed[[2]]^2 + parsed[[3]]^2)

        momentumList[[chunkIndex]] <- p
        idList[[chunkIndex]] <- ids[keep]
      }
    }
  }

  list(momentum = unlist(momentumList), id = unlist(idList))
}

# ---- countPerFile[file, pair, bin, 1] = WT count, [,,,2] = mutant count ----
countPerFile <- array(0, dim = c(nFiles, nPairs, nBins, 2))
breaks <- NULL

for (fileIndex in 1:nFiles) {
  filepath <- filenames[fileIndex]
  cat("Processing", filepath, "...\n")
  fileData <- extractMomentumAndID(filepath)

  # Set 1 also defines the bin boundaries (quantiles), reused for every file
  if (fileIndex == 1) {
    breaks <- quantile(fileData$momentum, probs = seq(0, 1, length.out = nBins + 1))
    breaks[1] <- -Inf
    breaks[length(breaks)] <- Inf
  }

  bin <- as.integer(cut(fileData$momentum, breaks = breaks, include.lowest = TRUE))

  for (i in 1:nPairs) {
    for (b in 1:nBins) {
      inBin <- (bin == b)
      countPerFile[fileIndex, i, b, 1] <-
        sum(fileData$id == pairPositiveID[i] & inBin)
      countPerFile[fileIndex, i, b, 2] <-
        sum(fileData$id == pairNegativeID[i] & inBin)
    }
  }
}

# ---- Cache the full per-file breakdown, so future changes don't need to
# re-read the raw data ----
perFileRows <- data.frame()
for (fileIndex in 1:nFiles) {
  for (i in 1:nPairs) {
    for (b in 1:nBins) {
      perFileRows <- rbind(perFileRows, data.frame(
        Set = fileIndex,
        Strain = pairName[i],
        MomentumBin = b,
        CountWT = countPerFile[fileIndex, i, b, 1],
        CountMutant = countPerFile[fileIndex, i, b, 2]
      ))
    }
  }
}
options(scipen = 999)
write.csv(perFileRows, "momentum_perfile_results.csv", row.names = FALSE)
cat("\nWrote momentum_perfile_results.csv (", nrow(perFileRows), "rows).\n")

# ---- Combine the 10 files per pair per bin, using the sub-sampling method
# (mean + SD of the 10 per-file asymmetries) ----
results <- data.frame()

for (i in 1:nPairs) {
  for (b in 1:nBins) {
    perFileAsymmetry <- numeric(0)

    for (fileIndex in 1:nFiles) {
      countWT <- countPerFile[fileIndex, i, b, 1]
      countMutant <- countPerFile[fileIndex, i, b, 2]
      total <- countWT + countMutant
      if (total > 0) {
        perFileAsymmetry <- c(perFileAsymmetry, (countWT - countMutant) / total)
      }
      # files with 0 total for this strain/bin are skipped (undefined 0/0)
    }

    nValidFiles <- length(perFileAsymmetry)
    if (nValidFiles >= 2) {
      asymmetry <- mean(perFileAsymmetry)
      uncertainty <- sd(perFileAsymmetry)
    } else if (nValidFiles == 1) {
      asymmetry <- perFileAsymmetry[1]
      uncertainty <- NA
    } else {
      asymmetry <- NA
      uncertainty <- NA
    }

    totalWT <- sum(countPerFile[, i, b, 1])
    totalMutant <- sum(countPerFile[, i, b, 2])

    results <- rbind(results, data.frame(
      Strain = pairName[i],
      MomentumBin = b,
      MomentumRangeLow = round(breaks[b], 4),
      MomentumRangeHigh = round(breaks[b + 1], 4),
      CountWT = totalWT,
      CountMutant = totalMutant,
      FilesUsed = nValidFiles,
      Asymmetry = round(asymmetry, 5),
      Uncertainty = round(uncertainty, 5)
    ))
  }
}

print(results)

write.csv(results, "momentum_asymmetry_results.csv", row.names = FALSE)
cat("\nWrote momentum_asymmetry_results.csv (", nrow(results), "rows).\n")
