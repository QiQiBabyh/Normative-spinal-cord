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
library(ggplot2)
library(doParallel)
library(foreach)
library(iterators)
library(openxlsx)
library(stringr)
library(gamlss)
library(reshape2)
library(patchwork)

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

average_by_age_grid <- function(values, num_length) {
  par <- values[1:num_length]
  seg <- length(values) / num_length
  if(seg > 1) {
    for(seg_i in 2:seg) {
      par <- par + values[((seg_i - 1) * num_length + 1):(seg_i * num_length)]
    }
    par <- par / seg
  }
  return(par)
}

predict_centile_grid <- function(model, reference_data, sex_levels, num_length = 5000) {
  age_grid <- seq(min(reference_data$Age), max(reference_data$Age), length.out = num_length)

  mu_coef <- get_random_coef(model, "mu")
  if(!is.null(mu_coef)) {
    new_mu <- expand.grid(Age = age_grid, Sex = sex_levels, Site_ZZZ = names(mu_coef))
  } else {
    new_mu <- expand.grid(Age = age_grid, Sex = sex_levels)
  }
  mu0 <- predict(model, newdata = new_mu, type = "response", what = "mu")

  sigma_coef <- get_random_coef(model, "sigma")
  if(!is.null(sigma_coef)) {
    new_sigma <- expand.grid(Age = age_grid, Sex = sex_levels, Site_ZZZ = names(sigma_coef))
  } else {
    new_sigma <- expand.grid(Age = age_grid, Sex = sex_levels)
  }
  sigma0 <- predict(model, newdata = new_sigma, type = "response", what = "sigma")

  nu_coef <- get_random_coef(model, "nu")
  if(!is.null(nu_coef)) {
    new_nu <- expand.grid(Age = age_grid, Sex = sex_levels, Site_ZZZ = names(nu_coef))
  } else {
    new_nu <- expand.grid(Age = age_grid, Sex = sex_levels)
  }
  nu0 <- predict(model, newdata = new_nu, type = "response", what = "nu")

  mu <- average_by_age_grid(mu0, num_length)
  sigma <- average_by_age_grid(sigma0, num_length)
  nu <- average_by_age_grid(nu0, num_length)

  centile_grid <- zzz_cent(
    obj = model,
    type = c("centiles"),
    mu = mu,
    sigma = sigma,
    nu = nu,
    cent = c(0.5, 2.5, 50, 97.5, 99.5),
    xname = "Age",
    xvalues = age_grid,
    calibration = FALSE,
    lpar = 3
  )
  centile_grid[, "sigma"] <- sigma
  return(centile_grid)
}

add_gradient <- function(centile_grid, reference_data, num_length = 5000) {
  colnames(centile_grid) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")
  step_age <- (max(reference_data$Age) - min(reference_data$Age)) / num_length
  gradient <- (centile_grid$median[2:nrow(centile_grid)] - centile_grid$median[1:(nrow(centile_grid) - 1)]) / step_age
  centile_grid$Gradient1 <- c(gradient, gradient[nrow(centile_grid) - 1])
  return(centile_grid)
}

fit_normative_model <- function(model_data, m0, i_rnd, j_rnd, con, include_sex = TRUE) {
  sex_term <- if(include_sex) "+ Sex" else ""

  if(i_rnd == 1) {
    mu_random <- "+ random(Site_ZZZ)"
  } else {
    mu_random <- ""
  }

  if(j_rnd == 1) {
    sigma_random <- "+ random(Site_ZZZ)"
  } else {
    sigma_random <- ""
  }

  mu_formula <- as.formula(paste0(
    "feature ~ bfpNA(Age, c(",
    paste(m0$mu.coefSmo[[1]]$power, collapse = ","),
    ")) ", sex_term, " ", mu_random
  ))

  sigma_formula <- as.formula(paste0(
    "feature ~ bfpNA(Age, c(",
    paste(m0$sigma.coefSmo[[1]]$power, collapse = ","),
    ")) ", sex_term, " ", sigma_random
  ))

  model <- gamlss(
    formula = mu_formula,
    sigma.formula = sigma_formula,
    control = con,
    family = GG(mu.link = "log", sigma.link = "log", nu.link = "identity"),
    data = model_data
  )

  return(model)
}

make_site_stratified_folds <- function(data, k = 10, strata_col = "Site_ZZZ", seed = 123) {
  set.seed(seed)
  folds <- vector("list", k)
  strata <- as.character(data[[strata_col]])
  strata[is.na(strata) | strata == ""] <- "UNKNOWN_SITE"

  for(s in unique(strata)) {
    idx <- which(strata == s)
    idx <- sample(idx)
    fold_id <- rep(seq_len(k), length.out = length(idx))
    for(ff in seq_len(k)) {
      folds[[ff]] <- c(folds[[ff]], idx[fold_id == ff])
    }
  }

  folds <- lapply(folds, sort)
  return(folds)
}

align_test_site_for_random_effect <- function(test_data, cv_model) {
  test_data$Site_ZZZ <- as.character(test_data$Site_ZZZ)

  coef_pool <- get_random_coef(cv_model, "mu")
  if(is.null(coef_pool)) {
    coef_pool <- get_random_coef(cv_model, "sigma")
  }

  if(!is.null(coef_pool)) {
    train_sites <- names(coef_pool)
    ref_site <- names(which.min(abs(coef_pool - mean(coef_pool, na.rm = TRUE))))
    unknown_site <- !(test_data$Site_ZZZ %in% train_sites)
    if(any(unknown_site)) {
      test_data$Site_ZZZ[unknown_site] <- ref_site
    }
  }

  test_data$Site_ZZZ <- as.factor(test_data$Site_ZZZ)
  test_data$Sex <- factor(test_data$Sex, levels = c("Female", "Male"))
  return(test_data)
}

run_hc_10fold_cv <- function(hc_data, feature_name, m0, i_rnd, j_rnd, con, k_cv = 10) {
  set.seed(123)
  folds_HC <- make_site_stratified_folds(hc_data, k = k_cv, strata_col = "Site_ZZZ", seed = 123)

  Z_score_folds_HC1 <- NULL
  Quant_score_folds_HC1 <- NULL
  HC_10fold_CV_failed <- data.frame()

  for(i_fold in seq_len(k_cv)) {
    cat("HC 10-fold CV:", feature_name, "fold", i_fold, "/", k_cv, "\n")

    test_idx <- folds_HC[[i_fold]]
    train_idx <- setdiff(seq_len(nrow(hc_data)), test_idx)

    train_data <- hc_data[train_idx, , drop = FALSE]
    test_data <- hc_data[test_idx, , drop = FALSE]

    train_data$Sex <- factor(train_data$Sex, levels = c("Female", "Male"))
    train_data$Site_ZZZ <- as.factor(train_data$Site_ZZZ)
    test_data$Sex <- factor(test_data$Sex, levels = c("Female", "Male"))
    test_data$Site_ZZZ <- as.factor(test_data$Site_ZZZ)

    cv_model_try <- try(
      fit_normative_model(
        model_data = train_data,
        m0 = m0,
        i_rnd = i_rnd,
        j_rnd = j_rnd,
        con = con,
        include_sex = TRUE
      ),
      silent = TRUE
    )

    if(inherits(cv_model_try, "try-error")) {
      warning(paste0("HC 10-fold CV failed for ", feature_name, ", fold ", i_fold))
      HC_10fold_CV_failed <- rbind(
        HC_10fold_CV_failed,
        data.frame(feature = feature_name, fold = i_fold, reason = as.character(cv_model_try))
      )
      next
    }

    cv_model <- cv_model_try
    test_data_pred <- align_test_site_for_random_effect(test_data, cv_model)

    mu_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "mu")
    sigma_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "sigma")
    nu_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "nu")

    Z_score_folds_HC <- zzz_cent(
      obj = cv_model, type = c("z-scores"),
      mu = mu_cv, sigma = sigma_cv, nu = nu_cv,
      xname = "Age", xvalues = test_data_pred$Age,
      yval = test_data_pred$feature,
      calibration = FALSE, lpar = 3
    )
    Z_score_folds_HC <- data.frame(Z_score = as.numeric(Z_score_folds_HC), Fold = i_fold)
    rownames(Z_score_folds_HC) <- rownames(test_data_pred)
    Z_score_folds_HC1 <- rbind(Z_score_folds_HC1, Z_score_folds_HC)

    Quant_score_folds_HC <- zzz_cent(
      obj = cv_model, type = c("z-scores"),
      mu = mu_cv, sigma = sigma_cv, nu = nu_cv,
      xname = "Age", xvalues = test_data_pred$Age,
      yval = test_data_pred$feature,
      calibration = FALSE, lpar = 3, cdf = TRUE
    )
    Quant_score_folds_HC <- data.frame(Quant_score = as.numeric(Quant_score_folds_HC), Fold = i_fold)
    rownames(Quant_score_folds_HC) <- rownames(test_data_pred)
    Quant_score_folds_HC1 <- rbind(Quant_score_folds_HC1, Quant_score_folds_HC)
  }

  if(!is.null(Z_score_folds_HC1) && nrow(Z_score_folds_HC1) > 0) {
    Z_score_folds_HC1 <- Z_score_folds_HC1[rownames(hc_data), , drop = FALSE]
  }
  if(!is.null(Quant_score_folds_HC1) && nrow(Quant_score_folds_HC1) > 0) {
    Quant_score_folds_HC1 <- Quant_score_folds_HC1[rownames(hc_data), , drop = FALSE]
  }

  summary <- data.frame(
    feature = feature_name,
    n_HC = nrow(hc_data),
    k = k_cv,
    n_success = ifelse(is.null(Z_score_folds_HC1), 0, sum(!is.na(Z_score_folds_HC1$Z_score))),
    n_failed_folds = nrow(HC_10fold_CV_failed),
    Z_mean = ifelse(is.null(Z_score_folds_HC1), NA, mean(Z_score_folds_HC1$Z_score, na.rm = TRUE)),
    Z_sd = ifelse(is.null(Z_score_folds_HC1), NA, sd(Z_score_folds_HC1$Z_score, na.rm = TRUE)),
    Z_median = ifelse(is.null(Z_score_folds_HC1), NA, median(Z_score_folds_HC1$Z_score, na.rm = TRUE)),
    Z_IQR = ifelse(is.null(Z_score_folds_HC1), NA, IQR(Z_score_folds_HC1$Z_score, na.rm = TRUE)),
    Quant_mean = ifelse(is.null(Quant_score_folds_HC1), NA, mean(Quant_score_folds_HC1$Quant_score, na.rm = TRUE)),
    Quant_sd = ifelse(is.null(Quant_score_folds_HC1), NA, sd(Quant_score_folds_HC1$Quant_score, na.rm = TRUE))
  )

  return(list(
    Z_score_folds_HC = Z_score_folds_HC1,
    Quant_score_folds_HC = Quant_score_folds_HC1,
    HC_10fold_CV_summary = summary,
    HC_10fold_CV_failed = HC_10fold_CV_failed,
    folds_HC = folds_HC
  ))
}

# ============================================================
# Load preprocessed clinical data
# ============================================================
setwd(datapath)
data1_original <- read_table_auto(clinical_datapath)

required_clinical_cols <- c(
  "Freesufer_Path2", "Freesufer_Path3", "Age", "Sex", "Site_ZZZ", "Diagnosis",
  "Database_included", "baseline", "Image_Quality_lab", "Data_baseline"
)

missing_clinical_cols <- setdiff(required_clinical_cols, colnames(data1_original))
if(length(missing_clinical_cols) > 0) {
  stop(paste("The clinical file is missing required columns:", paste(missing_clinical_cols, collapse = ", ")))
}

# ============================================================
# Build HC-based normative models
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
  figure_dir <- paste0(savepath, "/", str, "/figures")
  cv_dir <- paste0(savepath, "/", str, "/HC_10fold_CV")
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cv_dir, recursive = TRUE, showWarnings = FALSE)

  plot_list <- list()

  for(i in tem_feature) {
    print(i)

    model_file <- paste0(model_dir, "/", str, "_", i, "_normative_model.rds")
    if(file.exists(model_file)) {
      print("Normative model file already exists; skip this feature.")
      next
    }

    data1 <- data1_original %>%
      distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)

    site_count <- data1 %>%
      group_by(Site_ZZZ) %>%
      summarise(count = n(), .groups = "drop")
    site_count <- site_count[order(site_count$count), ]
    print(site_count)

    # Exclude sites with fewer than 10 participants from normative modeling.
    for(site in unique(site_count$Site_ZZZ)) {
      if(site_count[site_count$Site_ZZZ == site, "count"] < 10) {
        data1[data1$Site_ZZZ == site, "Database_included"] <- 0
      }
    }

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

    # Use eligible healthy controls to fit the normative model.
    hc_data <- all_data[!is.na(all_data$tem_feature) &
                          !is.na(all_data$Data_baseline) &
                          all_data$Diagnosis == "HC" &
                          all_data$Age >= 8 & all_data$Age <= 85, ]

    hc_data$Site_ZZZ <- as.factor(hc_data$Site_ZZZ)
    hc_data$Sex <- factor(hc_data$Sex, levels = c("Female", "Male"))
    hc_data <- hc_data[order(hc_data$Age), ]
    hc_data$feature <- hc_data$tem_feature

    hc_data <- hc_data[!is.na(hc_data$tem_feature), ]
    hc_data <- hc_data[hc_data$feature > (mean(hc_data$feature) - 3 * sd(hc_data$feature)) &
                         hc_data$feature < (mean(hc_data$feature) + 3 * sd(hc_data$feature)), ]
    hc_data <- hc_data[hc_data$feature > 0, ]
    hc_data <- hc_data[, c("Age", "Sex", "Site_ZZZ", "tem_feature", "feature")]

    # The source fit_model() function uses global objects named data1 and list_par.
    data1 <- hc_data

    list_par <- data.frame(matrix(0, 3 * 3 * 2 * 2, 4))
    colnames(list_par) <- c("mu_poly", "sigma_poly", "mu_random", "sigma_random")

    num <- 0
    for(i_poly in 1:3) {
      for(j_poly in 1:3) {
        for(i_rnd in 0:1) {
          for(j_rnd in 0:1) {
            num <- num + 1
            list_par[num, 1] <- i_poly
            list_par[num, 2] <- j_poly
            list_par[num, 3] <- i_rnd
            list_par[num, 4] <- j_rnd
          }
        }
      }
    }
    list_par <- na.omit(list_par)

    con <- gamlss.control()
    list_fit <- data.frame()
    n_cores <- max(1, parallel::detectCores() - 1)

    results_try <- try({
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      my_data <- foreach(
        num = 1:nrow(list_par),
        .combine = rbind,
        .packages = c("gamlss"),
        .errorhandling = "remove",
        .export = c("data1", "list_par", "fit_model")
      ) %dopar% {
        tryCatch({ fit_model(num) }, error = function(e) { return(NULL) })
      }
      stopCluster(cl)

      if(is.null(my_data) || nrow(my_data) == 0) {
        stop("All model-selection tasks failed. Please check the fit_model() function.")
      }

      list_fit <- my_data
      print(list_fit)
      model_ind <- which.min(list_fit$BIC)
      sel_mu_poly <- list_fit$mu_poly[model_ind]
      sel_sigma_poly <- list_fit$sigma_poly[model_ind]
      i_rnd <- list_fit$mu_random[model_ind]
      j_rnd <- list_fit$sigma_random[model_ind]
    }, silent = TRUE)

    if(inherits(results_try, "try-error")) {
      sel_mu_poly <- 2
      sel_sigma_poly <- 2
      i_rnd <- 1
      j_rnd <- 1
      con <- gamlss.control(c.crit = 0.01, n.cyc = 2, autostep = FALSE)
      list_fit <- data.frame(
        mu_poly = sel_mu_poly,
        sigma_poly = sel_sigma_poly,
        mu_random = i_rnd,
        sigma_random = j_rnd,
        BIC = NA_real_,
        note = "Fallback model used because automated model selection failed."
      )
    }

    m0 <- best_fit(sel_mu_poly, sel_sigma_poly, i_rnd, j_rnd)
    m2 <- fit_normative_model(data1, m0, i_rnd, j_rnd, con, include_sex = TRUE)
    m3 <- fit_normative_model(data1, m0, i_rnd, j_rnd, con, include_sex = FALSE)

    num_length <- 5000
    p2 <- predict_centile_grid(m3, data1, sex_levels = c("Female", "Male"), num_length = num_length)
    p2_all <- predict_centile_grid(m2, data1, sex_levels = c("Female", "Male"), num_length = num_length)
    male_p2 <- predict_centile_grid(m2, data1, sex_levels = c("Male"), num_length = num_length)
    female_p2 <- predict_centile_grid(m2, data1, sex_levels = c("Female"), num_length = num_length)

    p2 <- add_gradient(p2, data1, num_length)
    Male_p2 <- add_gradient(male_p2, data1, num_length)
    Female_p2 <- add_gradient(female_p2, data1, num_length)
    colnames(p2_all) <- c("Age", "lower99CI", "lower95CI", "median", "upper95CI", "upper99CI", "sigma")

    subtitle_map <- c(
      "CSA1" = "C1", "CSA2" = "C2", "CSA3" = "C3", "MUCCA" = "C1-3",
      "RLD1" = "C1", "RLD2" = "C2", "RLD3" = "C3", "RLD" = "C1-3",
      "APD1" = "C1", "APD2" = "C2", "APD3" = "C3", "APD" = "C1-3"
    )
    sub_title <- ifelse(i %in% names(subtitle_map), subtitle_map[i], i)

    if(str_detect(i, regex("CSA|MUCCA", ignore_case = TRUE))) {
      scale1 <- 1
      ylab1 <- "mm\u00B2"
      row_label <- "CSA"
    }
    if(str_detect(i, regex("^RLD", ignore_case = TRUE))) {
      scale1 <- 1
      ylab1 <- "mm"
      row_label <- "RLD"
    }
    if(str_detect(i, regex("^APD", ignore_case = TRUE))) {
      scale1 <- 1
      ylab1 <- "mm"
      row_label <- "APD"
    }

    left_col_features <- c("CSA1", "RLD1", "APD1")
    show_ylab <- i %in% left_col_features

    png(filename = paste0(figure_dir, "/", str, "_", i, "_all_with_sex_stratified.png"),
        width = 1480, height = 740, units = "px", bg = "white", res = 300)

    p_combined <- ggplot() +
      geom_point(data = data1[data1$Sex == "Female", ], aes(x = Age, y = tem_feature / scale1),
                 colour = "#E84935", shape = 16, size = 1.5, alpha = 0.08) +
      geom_point(data = data1[data1$Sex == "Male", ], aes(x = Age, y = tem_feature / scale1),
                 colour = "#4FBBD8", shape = 17, size = 1.5, alpha = 0.08) +
      geom_line(data = Female_p2, aes(x = Age, y = median / scale1),
                color = "#E84935", linewidth = 0.8, linetype = "solid") +
      geom_line(data = Female_p2, aes(x = Age, y = lower95CI / scale1),
                color = "#E84935", linewidth = 0.6, linetype = "dotted") +
      geom_line(data = Female_p2, aes(x = Age, y = upper95CI / scale1),
                color = "#E84935", linewidth = 0.6, linetype = "dotted") +
      geom_line(data = Male_p2, aes(x = Age, y = median / scale1),
                color = "#4FBBD8", linewidth = 0.8, linetype = "solid") +
      geom_line(data = Male_p2, aes(x = Age, y = lower95CI / scale1),
                color = "#4FBBD8", linewidth = 0.6, linetype = "dotted") +
      geom_line(data = Male_p2, aes(x = Age, y = upper95CI / scale1),
                color = "#4FBBD8", linewidth = 0.6, linetype = "dotted") +
      labs(title = sub_title, x = "", y = if(show_ylab) paste0(row_label, "\n(", ylab1, ")") else "") +
      theme_bw() +
      theme(
        plot.title = element_text(family = "Arial", size = 10, color = "black", hjust = 0.5, face = "bold"),
        axis.title.y = element_text(family = "Arial", size = 9, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(family = "Arial", size = 8, color = "black"),
        axis.text.y = element_text(family = "Arial", size = 8, color = "black"),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = c(8, 18, 60, 85), labels = c("8", "18", "60", "85"))

    print(p_combined)
    dev.off()
    plot_list[[i]] <- p_combined

    cv_results <- run_hc_10fold_cv(data1, i, m0, i_rnd, j_rnd, con, k_cv = 10)

    openxlsx::write.xlsx(
      list(
        Z_score_folds_HC = cv_results$Z_score_folds_HC,
        Quant_score_folds_HC = cv_results$Quant_score_folds_HC,
        HC_10fold_CV_summary = cv_results$HC_10fold_CV_summary,
        HC_10fold_CV_failed = cv_results$HC_10fold_CV_failed
      ),
      file = paste0(cv_dir, "/", str, "_", i, "_HC_10fold_CV.xlsx"),
      rowNames = TRUE,
      overwrite = TRUE
    )

    results <- list()
    results$feature_name <- i
    results$Female_p2 <- Female_p2
    results$Male_p2 <- Male_p2
    results$p2 <- p2
    results$peakage <- p2$Age[which.max(p2$median)]
    results$p2_all <- p2_all
    results$m2 <- m2
    results$m0 <- m0
    results$m3 <- m3
    results$list_fit <- list_fit
    results$hc_model_data <- data1
    results$str <- str
    results$Z_score_folds_HC <- cv_results$Z_score_folds_HC
    results$Quant_score_folds_HC <- cv_results$Quant_score_folds_HC
    results$HC_10fold_CV_summary <- cv_results$HC_10fold_CV_summary
    results$HC_10fold_CV_failed <- cv_results$HC_10fold_CV_failed
    results$folds_HC <- cv_results$folds_HC

    saveRDS(results, model_file)
  }

  row1_features <- c("CSA1", "CSA2", "CSA3", "MUCCA")
  row2_features <- c("RLD1", "RLD2", "RLD3", "RLD")
  row3_features <- c("APD1", "APD2", "APD3", "APD")
  all_features <- c(row1_features, row2_features, row3_features)

  missing_plots <- setdiff(all_features, names(plot_list))
  if(length(missing_plots) == 0) {
    row1_patch <- wrap_plots(plot_list[row1_features], nrow = 1)
    row2_patch <- wrap_plots(plot_list[row2_features], nrow = 1)
    row3_patch <- wrap_plots(plot_list[row3_features], nrow = 1)

    final_plot <- (row1_patch / row2_patch / row3_patch) +
      plot_annotation(
        caption = "Age (years)",
        theme = theme(plot.caption = element_text(family = "Arial", size = 10, hjust = 0.5, color = "black"))
      )

    png(filename = paste0(figure_dir, "/", str, "_combined_sex_stratified.png"),
        width = 3200, height = 2400, units = "px", bg = "white", res = 300)
    print(final_plot)
    dev.off()
  } else {
    warning(paste("Missing feature plots. The combined figure was not generated:", paste(missing_plots, collapse = ", ")))
  }
}
