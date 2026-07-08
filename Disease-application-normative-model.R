rm(list = ls())

# ============================================================
# User-defined paths
# ============================================================
datapath <- "E:/Lifespan_SCT_results/Github/Source-codes/"              # Directory containing the source R scripts
clinical_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/Clinical_vars.csv"  # Preprocessed clinical variables
MR_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/MR_measures.csv"           # Preprocessed spinal cord MRI measures
savepath <- "E:/Lifespan_SCT_results/Github/Test_results/"              # Output directory

# For the public spinal cord example, only the spinalcord dataset is used.
var <- c("spinalcord")

# If TRUE, only non-HC subjects are exported.
# If FALSE, all eligible baseline subjects are exported.
output_only_disease <- TRUE

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

  if(ext %in% c("xlsx", "xls")) {
    if(is.null(sheet)) {
      return(as.data.frame(readxl::read_excel(path)))
    } else {
      return(as.data.frame(readxl::read_excel(path, sheet = sheet)))
    }
  }

  if(ext == "csv") {
    return(read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  }

  stop("Unsupported file format. Please use a .csv, .xlsx, or .xls file.")
}

get_random_coef <- function(model, what = "mu") {
  smo_name <- paste0(what, ".coefSmo")
  smo <- model[[smo_name]]
  if(!is.null(smo) && length(smo) >= 1 && !is.null(smo[[1]]$coef)) {
    return(smo[[1]]$coef)
  }
  return(NULL)
}

align_site_for_random_effect <- function(new_data, model) {
  new_data$Site_ZZZ <- as.character(new_data$Site_ZZZ)

  coef_pool <- get_random_coef(model, "mu")
  if(is.null(coef_pool)) {
    coef_pool <- get_random_coef(model, "sigma")
  }

  if(!is.null(coef_pool)) {
    train_sites <- names(coef_pool)
    ref_site <- names(which.min(abs(coef_pool - mean(coef_pool, na.rm = TRUE))))
    unknown_site <- !(new_data$Site_ZZZ %in% train_sites)
    if(any(unknown_site)) {
      new_data$Site_ZZZ[unknown_site] <- ref_site
    }
  }

  new_data$Site_ZZZ <- as.factor(new_data$Site_ZZZ)
  new_data$Sex <- factor(new_data$Sex, levels = c("Female", "Male"))
  return(new_data)
}

calculate_deviation_for_feature <- function(all_data, model_object, feature_name) {
  model1 <- model_object$m2

  pred_data <- all_data[, c("Age", "Sex", "Site_ZZZ", "tem_feature")]
  pred_data$Sex <- factor(pred_data$Sex, levels = c("Female", "Male"))
  pred_data <- align_site_for_random_effect(pred_data, model1)

  mu <- predict(model1, newdata = pred_data, type = "response", what = "mu")
  sigma <- predict(model1, newdata = pred_data, type = "response", what = "sigma")
  nu <- predict(model1, newdata = pred_data, type = "response", what = "nu")

  if(length(mu) != nrow(pred_data)) {
    warning(paste0("Prediction length mismatch for ", feature_name))
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
# Load preprocessed clinical data
# ============================================================
setwd(datapath)
clinical_data <- read_table_auto(clinical_datapath)

required_clinical_cols <- c(
  "Freesufer_Path2", "Freesufer_Path3", "Age", "Sex", "Site_ZZZ", "Diagnosis",
  "Database_included", "baseline", "Image_Quality_lab", "Data_baseline"
)

missing_clinical_cols <- setdiff(required_clinical_cols, colnames(clinical_data))
if(length(missing_clinical_cols) > 0) {
  stop(paste("The clinical file is missing required columns:", paste(missing_clinical_cols, collapse = ", ")))
}

# ============================================================
# Calculate disease deviation scores using saved HC-based normative models
# ============================================================
for(sheet in var) {
  setwd(datapath)
  MRI <- read_table_auto(MR_datapath, sheet = sheet)

  required_mr_cols <- c("Freesufer_Path2", "Freesufer_Path3")
  missing_mr_cols <- setdiff(required_mr_cols, colnames(MRI))
  if(length(missing_mr_cols) > 0) {
    stop(paste("The MR-measure file is missing required columns:", paste(missing_mr_cols, collapse = ", ")))
  }

  MRI <- MRI %>%
    distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
  MRI <- MRI[!is.na(MRI$Freesufer_Path3), ]
  MRI <- as.data.frame(MRI)
  rownames(MRI) <- paste0(MRI$Freesufer_Path2, MRI$Freesufer_Path3)

  if(str_detect(sheet, "spinalcord")) {
    # Select the columns corresponding to spinal cord metrics.
    # In this dataset, these columns include CSA, RLD, APD, and MUCCA measures.
    spinalcord_feature_cols <- c(4:6, 8:10, 12:14, 19:24)
    if(max(spinalcord_feature_cols) > ncol(MRI)) {
      stop("The MR-measure file does not contain all expected spinal cord metric columns.")
    }
    tem_feature <- colnames(MRI)[spinalcord_feature_cols]
  } else {
    stop("Unsupported dataset name. Please use a spinalcord dataset or modify the feature-selection block.")
  }

  str <- sheet
  model_dir <- paste0(savepath, "/", str, "/normative_models")
  deviation_dir <- paste0(savepath, "/", str, "/deviation_scores")
  dir.create(deviation_dir, recursive = TRUE, showWarnings = FALSE)

  long_results <- list()
  wide_results <- NULL

  for(i in tem_feature) {
    print(paste("Calculating deviation scores for", i))

    model_file <- paste0(model_dir, "/", str, "_", i, "_normative_model.rds")
    if(!file.exists(model_file)) {
      stop(paste0("Normative model file not found: ", model_file,
                  "\nPlease run 01_build_HC_normative_models.R first."))
    }

    model_object <- readRDS(model_file)

    data1 <- clinical_data %>%
      distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
    rownames(data1) <- paste0(data1$Freesufer_Path2, data1$Freesufer_Path3)

    inter_row <- intersect(rownames(data1), rownames(MRI))
    data1 <- cbind(data1[inter_row, ], MRI[inter_row, i])
    colnames(data1)[ncol(data1)] <- "tem_feature"

    # Keep only included baseline scans that passed image-quality control.
    all_data <- data1[data1$Database_included == 1 &
                        data1$baseline == "baseline" &
                        (is.na(data1$Image_Quality_lab) | data1$Image_Quality_lab != 1), ]
    all_data <- all_data %>%
      distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
    rownames(all_data) <- paste0(all_data$Freesufer_Path2, all_data$Freesufer_Path3)

    all_data <- all_data[all_data$tem_feature != "" &
                           !is.na(all_data$tem_feature) &
                           !is.infinite(all_data$tem_feature) &
                           all_data$Age != "" &
                           !is.na(all_data$Age) &
                           !is.infinite(all_data$Age) &
                           all_data$Site_ZZZ != "" &
                           !is.na(all_data$Site_ZZZ) &
                           all_data$Sex != "" &
                           !is.na(all_data$Sex), ]

    if(output_only_disease) {
      all_data <- all_data[all_data$Diagnosis != "HC", ]
    }

    if(nrow(all_data) == 0) {
      warning(paste0("No eligible subjects found for ", i))
      next
    }

    feature_output <- calculate_deviation_for_feature(all_data, model_object, i)
    long_results[[i]] <- feature_output

    openxlsx::write.xlsx(
      feature_output,
      file = paste0(deviation_dir, "/", str, "_", i, "_deviation_scores.xlsx"),
      rowNames = FALSE,
      overwrite = TRUE
    )

    feature_wide <- feature_output[, c("ID", "Freesufer_Path2", "Freesufer_Path3", "Diagnosis", "Age", "Sex", "Site_ZZZ",
                                       "Raw_value", "Z_score", "Quant_score")]
    colnames(feature_wide)[colnames(feature_wide) == "Raw_value"] <- paste0(i, "_raw")
    colnames(feature_wide)[colnames(feature_wide) == "Z_score"] <- paste0(i, "_Zscore")
    colnames(feature_wide)[colnames(feature_wide) == "Quant_score"] <- paste0(i, "_Quant")

    if(is.null(wide_results)) {
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

  if(length(long_results) > 0) {
    long_results_df <- do.call(rbind, long_results)

    openxlsx::write.xlsx(
      list(
        deviation_scores_long = long_results_df,
        deviation_scores_wide = wide_results
      ),
      file = paste0(deviation_dir, "/", str, "_disease_deviation_scores_all_features.xlsx"),
      rowNames = FALSE,
      overwrite = TRUE
    )

    saveRDS(
      list(
        deviation_scores_long = long_results_df,
        deviation_scores_wide = wide_results
      ),
      file = paste0(deviation_dir, "/", str, "_disease_deviation_scores_all_features.rds")
    )
  }
}
