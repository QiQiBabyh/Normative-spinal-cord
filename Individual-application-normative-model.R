rm(list = ls())

# ============================================================
# User-defined paths
# ============================================================
# Directory containing the source R scripts
# Change to your own path

datapath <- "E:/Lifespan_SCT_results/Github/Source-codes/"

# Preprocessed clinical variables
clinical_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/Clinical_vars.csv"

# Preprocessed spinal cord MRI measures (.csv, .xlsx, or .xls)
MR_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/MR_measures.csv"

# Directory containing the saved HC-based normative models
# Change to the folder where your spinal cord normative model .rds files are stored
model_dir <- "E:/Lifespan_SCT_results/Github/Test_results/spinalcord/normative_models"

# Output directory for individual deviation scores
savepath <- "E:/Lifespan_SCT_results/Github/Test_results/Individual_spinalcord"

# If MR_datapath is an Excel file and spinal cord measures are in a specific sheet,
# set the sheet name here; otherwise keep NULL for CSV files.
mr_sheet <- NULL

# If TRUE, export only non-HC subjects.
# If FALSE, export all eligible baseline subjects.
output_only_disease <- FALSE

# Prefix used in normative model filenames.
# The script will try both:
#   spinalcord_<feature>_normative_model.rds
#   <feature>_normative_model.rds
model_prefix <- "spinalcord"

# ============================================================
# Load source functions
# ============================================================
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")
source("ZZZ_function.R")

# ============================================================
# Load packages
# ============================================================
library(readxl)
library(dplyr)
library(openxlsx)
library(stringr)
library(gamlss)

# ============================================================
# Helper functions
# ============================================================
read_table_auto <- function(path, sheet = NULL) {
  ext <- tolower(tools::file_ext(path))

  if (ext %in% c("xlsx", "xls")) {
    if (is.null(sheet)) {
      return(as.data.frame(readxl::read_excel(path)))
    } else {
      return(as.data.frame(readxl::read_excel(path, sheet = sheet)))
    }
  }

  if (ext == "csv") {
    return(read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }

  stop("Unsupported file format. Please use a .csv, .xlsx, or .xls file.")
}

get_random_coef <- function(model, what = "mu") {
  smo_name <- paste0(what, ".coefSmo")
  smo <- model[[smo_name]]
  if (!is.null(smo) && length(smo) >= 1 && !is.null(smo[[1]]) && !is.null(smo[[1]]$coef)) {
    return(smo[[1]]$coef)
  }
  return(NULL)
}

align_site_for_random_effect <- function(new_data, model) {
  new_data$Site_ZZZ <- as.character(new_data$Site_ZZZ)

  coef_pool <- get_random_coef(model, "mu")
  if (is.null(coef_pool)) {
    coef_pool <- get_random_coef(model, "sigma")
  }

  if (!is.null(coef_pool)) {
    train_sites <- names(coef_pool)
    ref_site <- names(which.min(abs(coef_pool - mean(coef_pool, na.rm = TRUE))))
    unknown_site <- !(new_data$Site_ZZZ %in% train_sites)
    if (any(unknown_site)) {
      new_data$Site_ZZZ[unknown_site] <- ref_site
    }
  }

  new_data$Site_ZZZ <- as.factor(new_data$Site_ZZZ)
  new_data$Sex <- factor(new_data$Sex, levels = c("Female", "Male"))
  return(new_data)
}

load_normative_model <- function(feature_name, model_dir, model_prefix = "spinalcord") {
  candidates <- c(
    file.path(model_dir, paste0(model_prefix, "_", feature_name, "_normative_model.rds")),
    file.path(model_dir, paste0(feature_name, "_normative_model.rds"))
  )

  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop(
      paste0(
        "Normative model file not found for feature: ", feature_name,
        "\nChecked: ", paste(candidates, collapse = " ; ")
      )
    )
  }

  readRDS(existing[1])
}

calculate_deviation_for_feature <- function(all_data, model_object, feature_name) {
  if (is.list(model_object) && "m2" %in% names(model_object)) {
    model1 <- model_object$m2
  } else {
    model1 <- model_object
  }

  pred_data <- all_data[, c("Age", "Sex", "Site_ZZZ", "tem_feature")]
  pred_data$Sex <- factor(pred_data$Sex, levels = c("Female", "Male"))
  pred_data <- align_site_for_random_effect(pred_data, model1)

  mu <- predict(model1, newdata = pred_data, type = "response", what = "mu")
  sigma <- predict(model1, newdata = pred_data, type = "response", what = "sigma")
  nu <- predict(model1, newdata = pred_data, type = "response", what = "nu")

  if (length(mu) != nrow(pred_data) || length(sigma) != nrow(pred_data) || length(nu) != nrow(pred_data)) {
    stop(paste0("Prediction length mismatch for feature: ", feature_name))
  }

  z_score <- zzz_cent(
    obj = model1,
    type = c("z-scores"),
    mu = mu,
    sigma = sigma,
    nu = nu,
    xname = "Age",
    xvalues = pred_data$Age,
    yval = pred_data$tem_feature,
    calibration = FALSE,
    lpar = 3
  )

  quant_score <- zzz_cent(
    obj = model1,
    type = c("z-scores"),
    mu = mu,
    sigma = sigma,
    nu = nu,
    xname = "Age",
    xvalues = pred_data$Age,
    yval = pred_data$tem_feature,
    calibration = FALSE,
    lpar = 3,
    cdf = TRUE
  )

  output <- data.frame(
    ID = rownames(all_data),
    Freesufer_Path2 = all_data$Freesufer_Path2,
    Freesufer_Path3 = all_data$Freesufer_Path3,
    Diagnosis = all_data$Diagnosis,
    Age = all_data$Age,
    Sex = all_data$Sex,
    Site_ZZZ = all_data$Site_ZZZ,
    Feature = feature_name,
    Raw_value = as.numeric(all_data$tem_feature),
    Z_score = as.numeric(z_score),
    Quant_score = as.numeric(quant_score),
    stringsAsFactors = FALSE
  )

  return(output)
}

# ============================================================
# Load clinical data
# ============================================================
setwd(datapath)
clinical_data <- read_table_auto(clinical_datapath)

required_clinical_cols <- c(
  "Freesufer_Path2", "Freesufer_Path3", "Age", "Sex", "Site_ZZZ", "Diagnosis",
  "Database_included", "baseline", "Image_Quality_lab"
)

missing_clinical_cols <- setdiff(required_clinical_cols, colnames(clinical_data))
if (length(missing_clinical_cols) > 0) {
  stop(
    paste(
      "The clinical file is missing required columns:",
      paste(missing_clinical_cols, collapse = ", ")
    )
  )
}

clinical_data <- clinical_data %>%
  distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
rownames(clinical_data) <- paste0(clinical_data$Freesufer_Path2, clinical_data$Freesufer_Path3)

# ============================================================
# Load spinal cord MRI measures
# ============================================================
MRI <- read_table_auto(MR_datapath, sheet = mr_sheet)

required_mr_cols <- c("Freesufer_Path2", "Freesufer_Path3")
missing_mr_cols <- setdiff(required_mr_cols, colnames(MRI))
if (length(missing_mr_cols) > 0) {
  stop(
    paste(
      "The MR-measure file is missing required columns:",
      paste(missing_mr_cols, collapse = ", ")
    )
  )
}

MRI <- MRI %>%
  distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
MRI <- MRI[!is.na(MRI$Freesufer_Path3), , drop = FALSE]
MRI <- as.data.frame(MRI)
rownames(MRI) <- paste0(MRI$Freesufer_Path2, MRI$Freesufer_Path3)

# ============================================================
# Select spinal cord features
# ============================================================
# Default spinal cord feature columns used in your earlier script:
# c(4:6, 8:10, 12:14, 19:24)
# Please verify that these columns still correspond to:
# CSA, APD, and RLD at C1/C2/C3 plus averaged C1-C3 measures.

spinalcord_feature_cols <- c(4:6, 8:10, 12:14, 19:24)
if (max(spinalcord_feature_cols) > ncol(MRI)) {
  stop("The MR-measure file does not contain all expected spinal cord metric columns.")
}
spinalcord_features <- colnames(MRI)[spinalcord_feature_cols]

message("Selected spinal cord features: ", paste(spinalcord_features, collapse = ", "))

# ============================================================
# Match clinical and MRI tables
# ============================================================
common_ids <- intersect(rownames(clinical_data), rownames(MRI))
if (length(common_ids) == 0) {
  stop("No overlapping subjects were found between the clinical and MRI tables.")
}

clinical_data <- clinical_data[common_ids, , drop = FALSE]
MRI <- MRI[common_ids, , drop = FALSE]

# ============================================================
# Output directory
# ============================================================
out_dir <- file.path(savepath, "spinalcord_individual_deviation_scores")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Calculate individual deviation scores
# ============================================================
long_results <- list()
wide_results <- NULL

for (feature_name in spinalcord_features) {
  message("Calculating individual deviation scores for ", feature_name)

  model_object <- load_normative_model(feature_name, model_dir, model_prefix = model_prefix)

  feature_data <- cbind(clinical_data, tem_feature = MRI[rownames(clinical_data), feature_name])

  # Keep only included baseline scans that passed the same filtering logic
  # as in your previous spinal cord script.
  all_data <- feature_data[
    feature_data$Database_included == 1 &
      feature_data$baseline == "baseline" &
      (is.na(feature_data$Image_Quality_lab) | feature_data$Image_Quality_lab != 1),
    , drop = FALSE
  ]

  all_data <- all_data %>%
    distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
  rownames(all_data) <- paste0(all_data$Freesufer_Path2, all_data$Freesufer_Path3)

  feature_numeric <- suppressWarnings(as.numeric(all_data$tem_feature))
  age_numeric <- suppressWarnings(as.numeric(all_data$Age))

  valid_idx <-
    !is.na(feature_numeric) & is.finite(feature_numeric) &
    !is.na(age_numeric) & is.finite(age_numeric) &
    !is.na(all_data$Site_ZZZ) & all_data$Site_ZZZ != "" &
    !is.na(all_data$Sex) & all_data$Sex %in% c("Female", "Male")

  all_data <- all_data[valid_idx, , drop = FALSE]
  all_data$tem_feature <- feature_numeric[valid_idx]
  all_data$Age <- age_numeric[valid_idx]

  if (output_only_disease) {
    all_data <- all_data[all_data$Diagnosis != "HC", , drop = FALSE]
  }

  if (nrow(all_data) == 0) {
    warning(paste0("No eligible subjects found for feature: ", feature_name))
    next
  }

  all_data <- all_data[order(all_data$Age), , drop = FALSE]

  feature_output <- calculate_deviation_for_feature(all_data, model_object, feature_name)
  long_results[[feature_name]] <- feature_output

  # Save per-feature tables
  openxlsx::write.xlsx(
    feature_output,
    file = file.path(out_dir, paste0("spinalcord_", feature_name, "_individual_deviation_scores.xlsx")),
    rowNames = FALSE,
    overwrite = TRUE
  )

  # Save per-feature RDS in a structure similar to the individual brain script
  per_feature_result <- list(
    model = if (is.list(model_object) && "m2" %in% names(model_object)) model_object$m2 else model_object,
    feature_output = feature_output,
    all_data = all_data,
    feature_name = feature_name,
    dataset = "spinalcord"
  )

  saveRDS(
    per_feature_result,
    file = file.path(out_dir, paste0("spinalcord_", feature_name, "_individual_deviation_scores.rds"))
  )

  # Build wide-format summary
  feature_wide <- feature_output[, c(
    "ID", "Freesufer_Path2", "Freesufer_Path3", "Diagnosis", "Age", "Sex", "Site_ZZZ",
    "Raw_value", "Z_score", "Quant_score"
  )]

  colnames(feature_wide)[colnames(feature_wide) == "Raw_value"] <- paste0(feature_name, "_raw")
  colnames(feature_wide)[colnames(feature_wide) == "Z_score"] <- paste0(feature_name, "_Zscore")
  colnames(feature_wide)[colnames(feature_wide) == "Quant_score"] <- paste0(feature_name, "_Quant")

  if (is.null(wide_results)) {
    wide_results <- feature_wide
  } else {
    wide_results <- merge(
      wide_results,
      feature_wide,
      by = c("ID", "Freesufer_Path2", "Freesufer_Path3", "Diagnosis", "Age", "Sex", "Site_ZZZ"),
      all = TRUE
    )
  }
}

# ============================================================
# Save combined outputs
# ============================================================
if (length(long_results) == 0) {
  stop("No deviation-score outputs were generated. Please check your filtering rules, feature columns, and model files.")
}

long_results_df <- do.call(rbind, long_results)

openxlsx::write.xlsx(
  list(
    deviation_scores_long = long_results_df,
    deviation_scores_wide = wide_results
  ),
  file = file.path(out_dir, "spinalcord_individual_deviation_scores_all_features.xlsx"),
  rowNames = FALSE,
  overwrite = TRUE
)

saveRDS(
  list(
    deviation_scores_long = long_results_df,
    deviation_scores_wide = wide_results,
    feature_names = spinalcord_features,
    output_only_disease = output_only_disease,
    mr_sheet = mr_sheet
  ),
  file = file.path(out_dir, "spinalcord_individual_deviation_scores_all_features.rds")
)

message("Finished. Individual spinal cord deviation scores were saved to: ", out_dir)
