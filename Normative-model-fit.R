rm(list = ls())

# ============================================================
# User-defined paths
# ============================================================
datapath <- "E:/Lifespan_SCT_results/Github/Source-codes/"              # Directory containing the source R scripts
clinical_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/Clinical_vars.csv"  # Preprocessed clinical variables
MR_datapath <- "E:/Lifespan_SCT_results/Github/Datasets/Dataset-norms/MR_measures.csv"           # Preprocessed spinal cord MRI measures
savepath <- "E:/Lifespan_SCT_results/Github/Test_results/"              # Output directory

# Select the MR-measure sheet(s) or dataset name(s) to be analyzed.
# For the public spinal cord example, only the spinalcord dataset is used.
var <- c("spinalcord")

# ============================================================
# Load required source functions
# ============================================================
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("300.variables.r")
source("301.functions.r")
source("ZZZ_function.R")

# ============================================================
# Load required packages
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
# Helper functions for loading preprocessed public data
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

# The clinical file is assumed to be preprocessed and to contain the required columns,
# including Site_ZZZ and quality-control fields. No additional numeric conversion,
# site-label construction, or local cohort exclusion is performed here.
setwd(datapath)
data1_original <- read_table_auto(clinical_datapath)

required_clinical_cols <- c(
  "Freesufer_Path2", "Freesufer_Path3", "Age", "Sex", "Site_ZZZ", "Diagnosis",
  "Database_included", "baseline", "Image_Quality_lab", "Data_baseline"
)

missing_clinical_cols <- setdiff(required_clinical_cols, colnames(data1_original))
if(length(missing_clinical_cols) > 0) {
  stop(paste(
    "The clinical file is missing required columns:",
    paste(missing_clinical_cols, collapse = ", ")
  ))
}

# ============================================================
# HC 10-fold cross-validation helper functions
# These functions are used only to generate out-of-fold Z-scores
# and quantiles in healthy controls. They do not change the final
# full-sample normative models or plotting workflow.
# ============================================================

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

fit_hc_cv_gamlss <- function(train_data, m0, i_rnd, j_rnd, con) {
  train_data$Sex <- factor(train_data$Sex, levels = c('Female', 'Male'))
  train_data$Site_ZZZ <- as.factor(train_data$Site_ZZZ)
  
  if(i_rnd == 1 & j_rnd == 1) {
    cv_model <- gamlss(
      formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
      sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
      control = con,
      family = GG(mu.link = 'log', sigma.link = 'log', nu.link = 'identity'),
      data = train_data
    )
  } else if(i_rnd == 1 & j_rnd == 0) {
    cv_model <- gamlss(
      formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
      sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex,
      control = con,
      family = GG(mu.link = 'log', sigma.link = 'log', nu.link = 'identity'),
      data = train_data
    )
  } else if(i_rnd == 0 & j_rnd == 1) {
    cv_model <- gamlss(
      formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex,
      sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex + random(Site_ZZZ),
      control = con,
      family = GG(mu.link = 'log', sigma.link = 'log', nu.link = 'identity'),
      data = train_data
    )
  } else if(i_rnd == 0 & j_rnd == 0) {
    cv_model <- gamlss(
      formula = feature ~ bfpNA(Age, c(m0$mu.coefSmo[[1]]$power)) + Sex,
      sigma.formula = feature ~ bfpNA(Age, c(m0$sigma.coefSmo[[1]]$power)) + Sex,
      control = con,
      family = GG(mu.link = 'log', sigma.link = 'log', nu.link = 'identity'),
      data = train_data
    )
  }
  
  return(cv_model)
}

align_test_site_for_random_effect <- function(test_data, cv_model) {
  test_data$Site_ZZZ <- as.character(test_data$Site_ZZZ)
  
  coef_pool <- NULL
  if(!is.null(cv_model$mu.coefSmo[[1]]) && !is.null(cv_model$mu.coefSmo[[1]]$coef)) {
    coef_pool <- cv_model$mu.coefSmo[[1]]$coef
  }
  if(is.null(coef_pool) && !is.null(cv_model$sigma.coefSmo[[1]]) && !is.null(cv_model$sigma.coefSmo[[1]]$coef)) {
    coef_pool <- cv_model$sigma.coefSmo[[1]]$coef
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
  test_data$Sex <- factor(test_data$Sex, levels = c('Female', 'Male'))
  return(test_data)
}
# ============================================================

for(sheet in var)
{ 
  setwd(datapath)
  MRI <- read_table_auto(MR_datapath, sheet = sheet)

  required_mr_cols <- c("Freesufer_Path2", "Freesufer_Path3")
  missing_mr_cols <- setdiff(required_mr_cols, colnames(MRI))
  if(length(missing_mr_cols) > 0) {
    stop(paste(
      "The MR-measure file is missing required columns:",
      paste(missing_mr_cols, collapse = ", ")
    ))
  }

  MRI_unique <- MRI %>%
    distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
  MRI<-MRI_unique
  
  MRI<-MRI[!is.na(MRI$Freesufer_Path3),]
  MRI<-as.data.frame(MRI)
  
  rownames(MRI)<-paste0(MRI$Freesufer_Path2,MRI$Freesufer_Path3)
  
  if(str_detect(sheet, 'spinalcord'))
  {
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
  
  if (!(dir.exists(paste0(savepath, '/', str)))) {
    dir.create(paste0(savepath, '/', str), recursive = TRUE)
  }
  
  setwd(paste0(savepath,'/',str))
  setwd(datapath)
  
  Z_data<-list()
  Quant_data<-list()
  
  # ============================================================
  # Initialize a list to store sex-stratified trajectory plots for each feature
  # ============================================================
  plot_list <- list()
  # ============================================================
  
  for(i in tem_feature[1:length(tem_feature)])
  {
    print(i)
    setwd(paste0(savepath,'/',str))
    rds_file <- paste0(str,'_',i,'_loop_our_model.rds')
    if(file.exists(rds_file)){
      old_results <- try(readRDS(rds_file), silent = TRUE)
      if(!inherits(old_results, 'try-error') &&
         !is.null(old_results$Z_score_folds_HC) &&
         !is.null(old_results$Quant_score_folds_HC)){
        print('file exist with HC 10-fold CV')
        next
      } else {
        print('file exist but without HC 10-fold CV; refit this feature')
      }
    }
    
    setwd(datapath)  
    data1 <- data1_original
    
    data1 <- data1 %>%
      distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
    
    library(dplyr)
    site_count <- data1 %>%
      group_by(Site_ZZZ) %>%  
      summarise(count = n())  
    site_count<-site_count[order(site_count$count),]
    print(site_count)
    
    # Exclude sites with fewer than 10 participants from normative modeling.
    for(site in unique(site_count$Site_ZZZ))
    {
      if(site_count[site_count$Site_ZZZ == site, 'count'] < 10)
      {
        data1[data1$Site_ZZZ == site, 'Database_included'] <- 0
      }
    }
    
    rownames(data1)<-paste0(data1$Freesufer_Path2,data1$Freesufer_Path3)
    
    setwd(paste0(savepath,'/',str))
    
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
    
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    
    
    # Keep only included baseline scans that passed image-quality control.
    all_data <- data1[data1$Database_included == 1 &
                        data1$baseline == "baseline" &
                        (is.na(data1$Image_Quality_lab) | data1$Image_Quality_lab != 1), ]
    all_data<- all_data %>%
      distinct(Freesufer_Path2, Freesufer_Path3, .keep_all = TRUE)
    
    rownames(all_data)<-paste0(all_data$Freesufer_Path2,all_data$Freesufer_Path3)
    
    all_data_original<-all_data
    
    data1<-all_data
    
    # Use eligible healthy controls to fit the normative model.
    data1 <- data1[!is.na(data1$tem_feature) &
                     !is.na(data1$Data_baseline) &
                     data1$Diagnosis == 'HC' &
                     data1$Age >= 8 & data1$Age <= 85, ]
    
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    data1$Sex<-as.factor(data1$Sex)
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))
    
    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature
    all_data[,'feature']<-all_data$tem_feature
    
    data1<-data1[!is.na(data1$tem_feature),]
    data1<-data1[data1$feature>(mean(data1$feature)-3*sd(data1$feature))&
                   data1$feature<(mean(data1$feature)+3*sd(data1$feature)),]
    
    data1<-data1[data1$feature>0,]
    data1<-data1[,c('Age','Sex','Site_ZZZ','tem_feature','feature')]
    
    
    list_par <- data.frame(matrix(0, 3*3*2*2, 4))
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
    
    print("Number of missing values in list_par:")
    print(colSums(is.na(list_par)))
    list_par <- na.omit(list_par)
    
    con <- gamlss.control()
    num <- 0
    list_fit <- data.frame()
    n_cores <- max(1, parallel::detectCores() - 1)
    
    results_try<-try({
      library(doParallel)
      library(foreach)
      
      cl<-makeCluster(n_cores)
      registerDoParallel(cl)
      
      my_data<-foreach(num=1:dim(list_par)[1],
                       .combine=rbind,
                       .packages = c('gamlss'),
                       .errorhandling = 'remove',
                       .export = c('data1', 'list_par', 'fit_model')) %dopar% {
                         tryCatch({
                           fit_model(num)
                         }, error = function(e) {
                           return(NULL)
                         })
                       }
      
      stopCluster(cl)
      
      if(is.null(my_data) || nrow(my_data) == 0) {
        stop("All model-selection tasks failed. Please check the fit_model() function.")
      }
      
      list_fit<-my_data
      print(list_fit)
      
      model_ind<-which.min(list_fit$BIC)
      sel_mu_poly=list_fit$mu_poly[model_ind]
      sel_sigma_poly=list_fit$sigma_poly[model_ind]
      i_rnd=list_fit$mu_random[model_ind]
      j_rnd=list_fit$sigma_random[model_ind]
    })
    
    if(inherits(results_try,'try-error')) 
    {
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
    
    
    m0<-best_fit(sel_mu_poly,sel_sigma_poly,i_rnd,j_rnd)    
    
    if(i_rnd==1&j_rnd==1){
      m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m2<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+Sex,
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+Sex,
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    
    if(i_rnd==1&j_rnd==1){
      m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                 sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data1)}else if(i_rnd==1&j_rnd==0){
                   m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power))+random(Site_ZZZ),
                              sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                              control=con,
                              family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                              data=data1)}else if(i_rnd==0&j_rnd==1){
                                m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                           sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power))+random(Site_ZZZ),
                                           control=con,
                                           family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                           data=data1)}else if(i_rnd==0&j_rnd==0){
                                             m3<-gamlss(formula=feature~bfpNA(Age,c(m0$mu.coefSmo[[1]]$power)),
                                                        sigma.formula = feature~bfpNA(Age,c(m0$sigma.coefSmo[[1]]$power)),
                                                        control=con,
                                                        family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                                                        data=data1)}
    
    # ============================================================
    # Prediction and plotting
    # ============================================================
    
    model1<-m3
    num_length=5000
    
    # --- Prediction using m3 without sex stratification ---
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male")) 
    }
    data4 <- do.call(what=expand.grid, args=data3)
    mu0 <- predict(model1, newdata=data4, type="response", what="mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    sigma0 <- predict(model1, newdata=data4, type="response", what="sigma")
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    nu0 <- predict(model1, newdata=data4, type="response", what="nu")
    
    # Average parameters across sex and site levels
    tem_par<-mu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; mu=par
    
    tem_par<-sigma0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; sigma=par
    
    tem_par<-nu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; nu=par
    
    p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                 cent=c(0.5,2.5,50,97.5,99.5),xname='Age',
                 xvalues=data4$Age[1:num_length],calibration=FALSE,lpar=3)
    p2[,'sigma']<-sigma
    
    # --- Prediction using m2 with sex included, averaged across levels ---
    model1<-m2
    
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    mu0 <- predict(model1, newdata=data4, type="response", what="mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    sigma0 <- predict(model1, newdata=data4, type="response", what="sigma")
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"),
                    Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female","Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    nu0 <- predict(model1, newdata=data4, type="response", what="nu")
    
    tem_par<-mu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; mu=par
    
    tem_par<-sigma0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; sigma=par
    
    tem_par<-nu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; nu=par
    
    p2_all<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                     cent=c(0.5,2.5,50,97.5,99.5),xname='Age',
                     xvalues=data4$Age[1:num_length],calibration=FALSE,lpar=3)
    p2_all[,'sigma']<-sigma
    
    # --- Prediction using m2 for males only ---
    model1<-m2
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"),
                    Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    mu0 <- predict(model1, newdata=data4, type="response", what="mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"),
                    Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    sigma0 <- predict(model1, newdata=data4, type="response", what="sigma")
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"),
                    Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Male"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    nu0 <- predict(model1, newdata=data4, type="response", what="nu")
    
    tem_par<-mu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; mu=par
    
    tem_par<-sigma0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; sigma=par
    
    tem_par<-nu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; nu=par
    
    male_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                      cent=c(0.5,2.5,50,97.5,99.5),xname='Age',
                      xvalues=data4$Age[1:num_length],calibration=FALSE,lpar=3)
    male_p2[,'sigma']<-sigma
    
    # --- Prediction using m2 for females only ---
    model1<-m2
    if(!is.null(model1$mu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"),
                    Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    mu0 <- predict(model1, newdata=data4, type="response", what="mu")
    
    if(!is.null(model1$sigma.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"),
                    Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    sigma0 <- predict(model1, newdata=data4, type="response", what="sigma")
    
    if(!is.null(model1$nu.coefSmo[[1]]$coef))
    {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"),
                    Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
    } else {
      data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),
                    Sex=c("Female"))
    }
    data4 <- do.call(what=expand.grid, args=data3)
    nu0 <- predict(model1, newdata=data4, type="response", what="nu")
    
    tem_par<-mu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; mu=par
    
    tem_par<-sigma0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; sigma=par
    
    tem_par<-nu0; par<-tem_par[1:num_length]; Seg=length(tem_par)/num_length
    if(Seg>1){ for(Seg1 in 2:Seg){ par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)] }; par=par/Seg }; nu=par
    
    female_p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                        cent=c(0.5,2.5,50,97.5,99.5),xname='Age',
                        xvalues=data4$Age[1:num_length],calibration=FALSE,lpar=3)
    female_p2[,'sigma']<-sigma
    
    # ============================================================
    # Format output column names
    # ============================================================
    library(reshape2)
    colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma')
    mydata<-melt(p2,id='Age'); colnames(mydata)<-c('Age','Percentile','Value')
    
    step_age<-(max(data1$Age)-min(data1$Age))/num_length
    
    Grad_p2<-(p2$median[2:dim(p2)[1]]-p2$median[1:(dim(p2)[1]-1)])/step_age
    Grad_p2<-data.frame(c(Grad_p2,Grad_p2[dim(p2)[1]-1]))
    p2<-cbind(p2,Grad_p2)
    colnames(p2)[dim(p2)[2]]<-c('Gradient1')
    
    colnames(female_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma')
    mydata<-melt(female_p2,id='Age'); colnames(mydata)<-c('Age','Percentile','Value')
    Female_p2<-female_p2
    
    colnames(male_p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma')
    mydata<-melt(male_p2,id='Age'); colnames(mydata)<-c('Age','Percentile','Value')
    Male_p2<-male_p2
    
    Grad_Female_p2<-(Female_p2$median[2:dim(Female_p2)[1]]-Female_p2$median[1:(dim(Female_p2)[1]-1)])/step_age
    Grad_Female_p2<-data.frame(c(Grad_Female_p2,Grad_Female_p2[dim(Female_p2)[1]-1]))
    Female_p2<-cbind(Female_p2,Grad_Female_p2)
    colnames(Female_p2)[dim(Female_p2)[2]]<-c('Gradient1')
    
    Grad_Male_p2<-(Male_p2$median[2:dim(Male_p2)[1]]-Male_p2$median[1:(dim(Male_p2)[1]-1)])/step_age
    Grad_Male_p2<-data.frame(c(Grad_Male_p2,Grad_Male_p2[dim(Male_p2)[1]-1]))
    Male_p2<-cbind(Male_p2,Grad_Male_p2)
    colnames(Male_p2)[dim(Male_p2)[2]]<-c('Gradient1')
    
    # ============================================================
    # Define units, y-axis labels, and panel titles
    # ============================================================
    
    # Mapping between feature names and panel titles
    subtitle_map <- c(
      'CSA1'  = 'C1',  'CSA2'  = 'C2',  'CSA3'  = 'C3',  'MUCCA' = 'C1-3',
      'RLD1'  = 'C1',  'RLD2'  = 'C2',  'RLD3'  = 'C3',  'RLD'   = 'C1-3',
      'APD1'  = 'C1',  'APD2'  = 'C2',  'APD3'  = 'C3',  'APD'   = 'C1-3'
    )
    sub_title <- ifelse(i %in% names(subtitle_map), subtitle_map[i], i)
    
    # Units and row labels
    if(str_detect(i, regex('CSA|MUCCA', ignore_case=TRUE)))
    {
      scale1    = 1
      ylab1     = 'mm\u00B2'
      row_label = 'CSA'
    }
    if(str_detect(i, regex('^RLD', ignore_case=TRUE)))
    {
      scale1    = 1
      ylab1     = 'mm'
      row_label = 'RLD'
    }
    if(str_detect(i, regex('^APD', ignore_case=TRUE)))
    {
      scale1    = 1
      ylab1     = 'mm'
      row_label = 'APD'
    }
    
    # Show the y-axis label only for the leftmost panel in each row
    left_col_features <- c('CSA1', 'RLD1', 'APD1')
    show_ylab <- i %in% left_col_features
    
    # ============================================================
    # Save the gradient plot without sex stratification
    # ============================================================
    png(filename=paste0(str,'_',i,'_all_without_sex_stratified_Gradient.png'),
        width=1480, height=740, units="px", bg="white", res=300)
    
    p3<-ggplot()+
      geom_line(data=p2, aes(x=Age, y=Gradient1/scale1),
                color='#262626', linewidth=1, linetype='solid')+
      labs(title=paste0(i,' ',ylab1), x='', y='')+
      theme_bw()+
      theme(
        plot.title  = element_text(family="Arial", size=12, color="black"),
        axis.title  = element_text(family="Arial", size=12, color="black"),
        axis.text.x = element_text(family="Arial", size=12, color="black"),
        axis.text.y = element_text(family="Arial", size=10, color="black")
      )+
      scale_x_continuous(breaks=c(8,18,60,85),
                         labels=c("8 yr","18 yr","60 yr","85 yr"))
    print(p3)
    dev.off()
    
    # ============================================================
    # Save the trajectory plot without sex stratification
    # ============================================================
    png(filename=paste0(str,'_',i,'_all_without_sex_stratified.png'),
        width=1480, height=740, units="px", bg="white", res=300)
    
    p3<-ggplot()+
      geom_point(data=data1[data1$Sex=='Female',], aes(x=Age, y=tem_feature/scale1),
                 colour='#E84935', shape=16, size=3, alpha=0.1)+
      geom_point(data=data1[data1$Sex=='Male',],   aes(x=Age, y=tem_feature/scale1),
                 colour='#4FBBD8', shape=17, size=3, alpha=0.1)+
      geom_line(data=p2, aes(x=Age, y=median/scale1),    color='#262626', linewidth=1, linetype='solid')+
      geom_line(data=p2, aes(x=Age, y=lower99CI/scale1), color='#262626', linewidth=1, linetype='dashed')+
      geom_line(data=p2, aes(x=Age, y=lower95CI/scale1), color='#262626', linewidth=1, linetype='dotted')+
      geom_line(data=p2, aes(x=Age, y=upper95CI/scale1), color='#262626', linewidth=1, linetype='dotted')+
      geom_line(data=p2, aes(x=Age, y=upper99CI/scale1), color='#262626', linewidth=1, linetype='dashed')+
      labs(title=paste0(i,' ',ylab1), x='', y='')+
      theme_bw()+
      theme(
        plot.title  = element_text(family="Arial", size=12, color="black"),
        axis.title  = element_text(family="Arial", size=12, color="black"),
        axis.text.x = element_text(family="Arial", size=12, color="black"),
        axis.text.y = element_text(family="Arial", size=10, color="black")
      )+
      scale_x_continuous(breaks=c(8,18,60,85),
                         labels=c("8 yr","18 yr","60 yr","85 yr"))
    print(p3)
    dev.off()
    
    # ============================================================
    # Save the sex-stratified gradient plot
    # ============================================================
    png(filename=paste0(str,'_',i,'_all_with_sex_stratified_Gradient.png'),
        width=1480, height=740, units="px", bg="white", res=300)
    
    p3<-ggplot()+
      geom_line(data=Female_p2, aes(x=Age, y=Gradient1/scale1),
                color='#E84935', linewidth=1, linetype='solid')+
      geom_line(data=Male_p2,   aes(x=Age, y=Gradient1/scale1),
                color='#4FBBD8', linewidth=1, linetype='solid')+
      labs(title=paste0(i,' ',ylab1), x='', y='')+
      theme_bw()+
      theme(
        plot.title  = element_text(family="Arial", size=12, color="black"),
        axis.title  = element_text(family="Arial", size=12, color="black"),
        axis.text.x = element_text(family="Arial", size=12, color="black"),
        axis.text.y = element_text(family="Arial", size=10, color="black")
      )+
      scale_x_continuous(breaks=c(8,18,60,85),
                         labels=c("8 yr","18 yr","60 yr","85 yr"))
    print(p3)
    dev.off()
    
    # ============================================================
    # Save the sex-stratified sigma plot
    # ============================================================
    png(filename=paste0(str,'_',i,'_all_with_sex_stratified_sigma.png'),
        width=1480, height=740, units="px", bg="white", res=300)
    
    p3<-ggplot()+
      geom_line(data=Female_p2, aes(x=Age, y=sigma),
                color='#E84935', linewidth=1, linetype='solid')+
      geom_line(data=Male_p2,   aes(x=Age, y=sigma),
                color='#4FBBD8', linewidth=1, linetype='solid')+
      labs(title=paste0(i,' ',ylab1), x='', y='')+
      theme_bw()+
      theme(
        plot.title  = element_text(family="Arial", size=12, color="black"),
        axis.title  = element_text(family="Arial", size=12, color="black"),
        axis.text.x = element_text(family="Arial", size=12, color="black"),
        axis.text.y = element_text(family="Arial", size=10, color="black")
      )+
      scale_x_continuous(breaks=c(8,18,60,85),
                         labels=c("8 yr","18 yr","60 yr","85 yr"))
    print(p3)
    dev.off()
    
    # ============================================================
    # Save the sex-stratified trajectory plot and store it for the combined figure
    # ============================================================
    png(filename=paste0(str,'_',i,'_all_with_sex_stratified.png'),
        width=1480, height=740, units="px", bg="white", res=300)
    
    p_combined <- ggplot()+
      geom_point(data=data1[data1$Sex=='Female',], aes(x=Age, y=tem_feature/scale1),
                 colour='#E84935', shape=16, size=1.5, alpha=0.08)+
      geom_point(data=data1[data1$Sex=='Male',],   aes(x=Age, y=tem_feature/scale1),
                 colour='#4FBBD8', shape=17, size=1.5, alpha=0.08)+
      geom_line(data=Female_p2, aes(x=Age, y=median/scale1),
                color='#E84935', linewidth=0.8, linetype='solid')+
      geom_line(data=Female_p2, aes(x=Age, y=lower95CI/scale1),
                color='#E84935', linewidth=0.6, linetype='dotted')+
      geom_line(data=Female_p2, aes(x=Age, y=upper95CI/scale1),
                color='#E84935', linewidth=0.6, linetype='dotted')+
      geom_line(data=Male_p2,   aes(x=Age, y=median/scale1),
                color='#4FBBD8', linewidth=0.8, linetype='solid')+
      geom_line(data=Male_p2,   aes(x=Age, y=lower95CI/scale1),
                color='#4FBBD8', linewidth=0.6, linetype='dotted')+
      geom_line(data=Male_p2,   aes(x=Age, y=upper95CI/scale1),
                color='#4FBBD8', linewidth=0.6, linetype='dotted')+
      # Center the panel title and show the y-axis label only in the leftmost panel
      labs(title = sub_title,
           x     = '',
           y     = if(show_ylab) paste0(row_label,'\n(',ylab1,')') else '')+
      theme_bw()+
      theme(
        plot.title       = element_text(family="Arial", size=10, color="black",
                                        hjust=0.5, face="bold"),
        axis.title.y     = element_text(family="Arial", size=9,  color="black"),
        axis.title.x     = element_blank(),
        axis.text.x      = element_text(family="Arial", size=8,  color="black"),
        axis.text.y      = element_text(family="Arial", size=8,  color="black"),
        panel.grid.minor = element_blank()
      )+
      scale_x_continuous(breaks=c(8,18,60,85),
                         labels=c("8","18","60","85"))
    
    print(p_combined)
    dev.off()
    
    # Store the plot for the combined figure
    plot_list[[i]] <- p_combined
    
    # ============================================================
    # Calculate Z-scores and quantiles
    # ============================================================
    Z_score_sum   <- NULL
    Quant_score_sum <- NULL
    
    all_data<-all_data[all_data$tem_feature!=''&
                         !is.null(all_data$tem_feature)&
                         !is.na(all_data$tem_feature)&
                         !is.infinite(all_data$tem_feature),]
    all_data<-all_data[all_data$Age!=''&
                         !is.null(all_data$Age)&
                         !is.na(all_data$Age)&
                         !is.infinite(all_data$Age),]
    all_data<-all_data[all_data$Site_ZZZ!=''&
                         !is.null(all_data$Site_ZZZ)&
                         !is.na(all_data$Site_ZZZ)&
                         !is.infinite(all_data$Site_ZZZ),]
    all_data<-all_data[all_data$Sex!=''&
                         !is.null(all_data$Sex)&
                         !is.na(all_data$Sex)&
                         !is.infinite(all_data$Sex),]
    
    all_data1<-all_data[,c('Age','Sex','Site_ZZZ','tem_feature')]
    all_data1$Sex <- factor(all_data1$Sex, levels = c('Female', 'Male'))
    all_data1$Site_ZZZ <- as.character(all_data1$Site_ZZZ)
    
    model1<-m2
    
    for(sub in 1:dim(all_data)[1])
    {
      if(!is.null(m2$mu.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$mu.coefSmo[[1]]$coef)))
        {
          all_data1$Site_ZZZ[sub]<-names(which.min(abs(m2$mu.coefSmo[[1]]$coef - mean(m2$mu.coefSmo[[1]]$coef, na.rm = TRUE))))
        }
      }
      if(!is.null(m2$sigma.coefSmo[[1]]))
      {
        if(!(all_data1$Site_ZZZ[sub] %in% names(m2$sigma.coefSmo[[1]]$coef)))
        {
          all_data1$Site_ZZZ[sub]<-names(which.min(abs(m2$sigma.coefSmo[[1]]$coef - mean(m2$sigma.coefSmo[[1]]$coef, na.rm = TRUE))))
        }
      }
    }
    
    all_data1$Site_ZZZ <- as.factor(all_data1$Site_ZZZ)
    
    mu    <- predict(model1, newdata=all_data1, type="response", what="mu")
    sigma <- predict(model1, newdata=all_data1, type="response", what="sigma")
    nu    <- predict(model1, newdata=all_data1, type="response", what="nu")
    
    if(length(mu)!=dim(all_data1)[1])
    {
      print("Error, Please Check Data!!!")
    }
    
    Z_score_sum<-zzz_cent(obj=model1, type=c("z-scores"),
                          mu=mu, sigma=sigma, nu=nu,
                          xname='Age', xvalues=all_data1$Age,
                          yval=all_data1$tem_feature,
                          calibration=FALSE, lpar=3)
    
    Quant_score_sum<-zzz_cent(obj=model1, type=c("z-scores"),
                              mu=mu, sigma=sigma, nu=nu,
                              xname='Age', xvalues=all_data1$Age,
                              yval=all_data1$tem_feature,
                              calibration=FALSE, lpar=3, cdf=TRUE)
    
    Z_score_sum<-data.frame(Z_score_sum)
    colnames(Z_score_sum)<-c('Z_score')
    rownames(Z_score_sum)<-rownames(all_data1)
    
    Quant_score_sum<-data.frame(Quant_score_sum)
    colnames(Quant_score_sum)<-c('Quant_score')
    rownames(Quant_score_sum)<-rownames(all_data1)
    
    Z_data[[i]]    <-Z_score_sum
    Quant_data[[i]]<-Quant_score_sum

    # ============================================================
    # HC 10-fold cross-validation
    # Generate out-of-fold Z-scores and quantiles in HCs for calibration and stability checks
    # Only the HC data used for normative modeling are used here; the final full-HC models and plots are unchanged
    # ============================================================
    set.seed(123)
    k_cv <- 10
    folds_HC <- make_site_stratified_folds(data1, k = k_cv, strata_col = "Site_ZZZ", seed = 123)
    
    Z_score_folds_HC1 <- NULL
    Quant_score_folds_HC1 <- NULL
    HC_10fold_CV_failed <- data.frame()
    
    for(i_fold in seq_len(k_cv)) {
      cat("HC 10-fold CV:", i, "fold", i_fold, "/", k_cv, "\n")
      
      test_idx <- folds_HC[[i_fold]]
      train_idx <- setdiff(seq_len(nrow(data1)), test_idx)
      
      train_data <- data1[train_idx, , drop = FALSE]
      test_data  <- data1[test_idx,  , drop = FALSE]
      
      train_data$Sex <- factor(train_data$Sex, levels = c('Female', 'Male'))
      train_data$Site_ZZZ <- as.factor(train_data$Site_ZZZ)
      test_data$Sex <- factor(test_data$Sex, levels = c('Female', 'Male'))
      test_data$Site_ZZZ <- as.factor(test_data$Site_ZZZ)
      
      cv_model_try <- try(
        fit_hc_cv_gamlss(
          train_data = train_data,
          m0 = m0,
          i_rnd = i_rnd,
          j_rnd = j_rnd,
          con = con
        ),
        silent = TRUE
      )
      
      if(inherits(cv_model_try, 'try-error')) {
        warning(paste0('HC 10-fold CV failed for ', i, ', fold ', i_fold))
        HC_10fold_CV_failed <- rbind(
          HC_10fold_CV_failed,
          data.frame(feature = i, fold = i_fold, reason = as.character(cv_model_try))
        )
        next
      }
      
      cv_model <- cv_model_try
      test_data_pred <- align_test_site_for_random_effect(test_data, cv_model)
      
      mu_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "mu")
      sigma_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "sigma")
      nu_cv <- predict(cv_model, newdata = test_data_pred, type = "response", what = "nu")
      
      if(length(mu_cv) != nrow(test_data_pred)) {
        warning(paste0('Prediction length mismatch for ', i, ', fold ', i_fold))
      }
      
      Z_score_folds_HC <- zzz_cent(
        obj = cv_model, type = c("z-scores"),
        mu = mu_cv, sigma = sigma_cv, nu = nu_cv,
        xname = 'Age', xvalues = test_data_pred$Age,
        yval = test_data_pred$feature,
        calibration = FALSE, lpar = 3
      )
      Z_score_folds_HC <- data.frame(Z_score = as.numeric(Z_score_folds_HC), Fold = i_fold)
      rownames(Z_score_folds_HC) <- rownames(test_data_pred)
      Z_score_folds_HC1 <- rbind(Z_score_folds_HC1, Z_score_folds_HC)
      
      Quant_score_folds_HC <- zzz_cent(
        obj = cv_model, type = c("z-scores"),
        mu = mu_cv, sigma = sigma_cv, nu = nu_cv,
        xname = 'Age', xvalues = test_data_pred$Age,
        yval = test_data_pred$feature,
        calibration = FALSE, lpar = 3, cdf = TRUE
      )
      Quant_score_folds_HC <- data.frame(Quant_score = as.numeric(Quant_score_folds_HC), Fold = i_fold)
      rownames(Quant_score_folds_HC) <- rownames(test_data_pred)
      Quant_score_folds_HC1 <- rbind(Quant_score_folds_HC1, Quant_score_folds_HC)
    }
    
    if(!is.null(Z_score_folds_HC1) && nrow(Z_score_folds_HC1) > 0) {
      Z_score_folds_HC1 <- Z_score_folds_HC1[rownames(data1), , drop = FALSE]
    }
    if(!is.null(Quant_score_folds_HC1) && nrow(Quant_score_folds_HC1) > 0) {
      Quant_score_folds_HC1 <- Quant_score_folds_HC1[rownames(data1), , drop = FALSE]
    }
    
    HC_10fold_CV_summary <- data.frame(
      feature = i,
      n_HC = nrow(data1),
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
    
    # Save the HC 10-fold cross-validation results for each feature
    try({
      openxlsx::write.xlsx(
        list(
          Z_score_folds_HC = Z_score_folds_HC1,
          Quant_score_folds_HC = Quant_score_folds_HC1,
          HC_10fold_CV_summary = HC_10fold_CV_summary,
          HC_10fold_CV_failed = HC_10fold_CV_failed
        ),
        file = paste0(str, '_', i, '_HC_10fold_CV.xlsx'),
        rowNames = TRUE,
        overwrite = TRUE
      )
    }, silent = TRUE)
    # ============================================================

    results<-list()
    results$Female_p2       <- Female_p2
    results$Male_p2         <- Male_p2
    results$p2              <- p2
    results$peakage         <- p2$Age[which.max(p2$median)]
    results$p2_all          <- p2_all
    results$m2              <- m2
    results$m0              <- m0
    results$m3              <- m3
    results$list_fit        <- list_fit
    results$Zscore          <- Z_data
    results$Quant_data      <- Quant_data
    results$data1           <- data1
    results$all_data        <- all_data1
    results$str             <- str
    results$i               <- i
    results$all_data_original <- all_data_original
    results$Z_score_folds_HC <- Z_score_folds_HC1
    results$Quant_score_folds_HC <- Quant_score_folds_HC1
    results$HC_10fold_CV_summary <- HC_10fold_CV_summary
    results$HC_10fold_CV_failed <- HC_10fold_CV_failed
    results$folds_HC <- folds_HC
    
    saveRDS(results, paste0(str,'_',i,'_loop_our_model.rds'))
    
  } # End of the feature loop
  
  # ============================================================
  # Combine all panels into a 3 x 4 figure after the feature loop
  # ============================================================
  library(patchwork)
  
  row1_features <- c('CSA1',  'CSA2',  'CSA3',  'MUCCA')
  row2_features <- c('RLD1',  'RLD2',  'RLD3',  'RLD')
  row3_features <- c('APD1',  'APD2',  'APD3',  'APD')
  all_features  <- c(row1_features, row2_features, row3_features)
  
  # Check whether all expected feature plots were generated
  missing_plots <- setdiff(all_features, names(plot_list))
  if(length(missing_plots) > 0){
    warning(paste("Missing feature plots. Please check:", paste(missing_plots, collapse = ", ")))
    message("The combined 3 x 4 figure was not generated because one or more panels are missing.")
    next
  }
  
  # Combine panels by row
  row1_patch <- wrap_plots(plot_list[row1_features], nrow=1)
  row2_patch <- wrap_plots(plot_list[row2_features], nrow=1)
  row3_patch <- wrap_plots(plot_list[row3_features], nrow=1)
  
  # Combine the three rows and add a shared x-axis label
  final_plot <- (row1_patch / row2_patch / row3_patch) +
    plot_annotation(
      caption = 'Age (years)',
      theme   = theme(
        plot.caption = element_text(family="Arial", size=10,
                                    hjust=0.5, color="black")
      )
    )
  
  # Save the combined figure
  setwd(paste0(savepath,'/',str))
  png(filename = paste0(str,'_combined_sex_stratified.png'),
      width  = 3200,
      height = 2400,
      units  = "px",
      bg     = "white",
      res    = 300)
  print(final_plot)
  dev.off()
  
  message("Combined figure saved: ", paste0(str, "_combined_sex_stratified.png"))
  # ============================================================
  
} # End of the sheet loop 
 
