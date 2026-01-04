library(tidyverse)
library(glmnet)
library(sf)
library(sp)
library(raster)
library(blockCV)

set.seed(123)


nov_df <- readr::read_csv("../../data/data/nov_df_reduced.csv") %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::filter(!is.na(avg_Log10_NoV))

#coverage probability function-------------------------------------------------
COV <- function(z, lower=NULL, upper=NULL, coverage=NULL) {
  if(!is.null(lower) && !is.null(upper)){
    z <- as.matrix(z)
    lower <- as.matrix(lower)
    upper <- as.matrix(upper)
    x <- z>=lower & z<=upper
    u <- x[!is.na(x)]
    round(sum(u)/length(u) * 100, 4)
  }else if(!is.null(coverage)){
    round(mean(coverage, na.rm = T),4)
  }
}

#scale function-------------------------------------------------
covariates_to_scale <- c("school_den", "carehome_den", "imd", "bame", "mob_7day", "prop_agri", "prop_urb")

scale_covariates <- function(df, covariates_to_scale) {
  # Define a custom scaling function
  scale <- function(x) {
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }
  
  # Apply the custom scaling function across specified covariates
  df <- df %>% mutate(across(all_of(covariates_to_scale), scale, .names = "scale_{.col}"))
  
  return(df)
}

nov_df<- scale_covariates(nov_df, covariates_to_scale)
nov_df<- as.data.frame(nov_df)

# ensure no missing coordinates
coords_master <- nov_df %>%
  filter(!is.na(Easting) & !is.na(Northing)) %>%
  distinct(site_code, Easting, Northing)

nov_df <- nov_df %>%
  dplyr::select(-Easting, -Northing) %>%
  left_join(coords_master, by = "site_code")

nov <- nov_df
nov <- st_as_sf(nov_df, coords= c("Easting", "Northing"), crs= 27700)
nov <- st_transform(nov,  crs = "+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +units=km +no_defs")


nov <- nov %>%
  arrange(site_code) %>%
  mutate(site_code= as.factor(site_code)) %>%
  mutate(s_index=as.numeric(site_code)) %>%
  mutate(site_code= as.character(site_code))

nov$Week_date<- as.Date(nov$Week_date)

nov<- nov %>%
  arrange(Week_date) %>%
  mutate(f_index= as.numeric(Week_date))

nov <- nov %>%
  mutate(avg_Log10_NoV = ifelse(GI_ND | GII_ND, NA, avg_Log10_NoV))



# create fold blocks------------------------------------------------------------
# unique sites

sites<- nov %>%
  dplyr::select(site_code) %>%
  unique()

# shapefile
england<- st_read("../../data/data/england_crop.shp")
england<- st_transform(england, crs= st_crs(nov))

raster<- raster(england)

folds<- cv_spatial(
  x= sites,
  r= raster,
  k= 10,
  seed= 12,
  plot = FALSE
)

fold_blocks<- folds$blocks

fold_blocks<- st_transform(fold_blocks, crs= st_crs(nov))

# assign folds to data  
nov <- nov %>%
  st_intersection(fold_blocks)

#cross-validation---------------------------------------------------
fit <- list()
metrics<- data.frame()
predictions <- data.frame()

for (k in 1:10) {
  
    # subset data
    train <- subset(nov, folds != k)  # All folds except k
    val <- subset(nov, folds == k)   # Only fold 
    
    vars_used <- c(
      "avg_Log10_NoV","lockdown_step3","lockdown_step4","lockdown_planB","scale_school_den", "scale_carehome_den", "scale_mob_7day",
      "scale_bame", "scale_imd", "scale_prop_agri", "scale_prop_urb","temp_2", "temp_3", "temp_4", "temp_5",
      "rain_2", "rain_3", "rain_4", "rain_5","week", "site_code")

    train <- train %>%
      st_drop_geometry() %>%
      dplyr::filter(complete.cases(dplyr::select(., all_of(vars_used))))

    val <- val %>%
      st_drop_geometry() %>%
      dplyr::filter(complete.cases(dplyr::select(., all_of(vars_used))))


    y_train <- train$avg_Log10_NoV
    X_train <- model.matrix(~ -1 + lockdown_step3 + lockdown_step4 + lockdown_planB + scale_school_den + scale_carehome_den + scale_mob_7day + scale_bame + scale_imd + scale_prop_agri + scale_prop_urb + 
      temp_2 + temp_3 + temp_4 + temp_5 + rain_2 + rain_3 + rain_4 + rain_5 + week * site_code, data = train)
                              

    X_test <- model.matrix( ~ -1 +lockdown_step3 + lockdown_step4 + lockdown_planB + scale_school_den + scale_carehome_den + scale_mob_7day + scale_bame + scale_imd + scale_prop_agri + scale_prop_urb + 
      temp_2 + temp_3 + temp_4 + temp_5 + rain_2 + rain_3 + rain_4 + rain_5 + week * site_code, data = val)                       

    y_test <- val$avg_Log10_NoV   

    # Align test matrix columns to training matrix
    missing_cols <- setdiff(colnames(X_train), colnames(X_test))
    extra_cols <- setdiff(colnames(X_test), colnames(X_train))

    # Remove extra columns from test set if any
    if(length(extra_cols) > 0) {
      X_test <- X_test[, !(colnames(X_test) %in% extra_cols), drop = FALSE]
    }

    # Add missing columns as zeros
    if(length(missing_cols) > 0) {
      zero_mat <- matrix(0, nrow = nrow(X_test), ncol = length(missing_cols))
      colnames(zero_mat) <- missing_cols
      X_test <- cbind(X_test, zero_mat)
    }

    # Reorder columns to match training set
    X_test <- X_test[, colnames(X_train), drop = FALSE]

    fit[[k]] <- cv.glmnet(X_train, y_train, alpha = 1, family = "gaussian")                    

    fit.fold<- fit[[k]]
  
    print(summary(fit.fold))

    best_lambda <- fit[[k]]$lambda.min
    
    # Predictions
    val$predicted <- predict(fit.fold, s= best_lambda, newx = X_test)

    # bootstrap predictions for uncertainty estimation
    n_boot <- 500
    boot_preds <- matrix(NA, nrow = n_boot, ncol = nrow(val))

    for (i in 1:n_boot) {
      set.seed(123 + i)
      boot_idx <- sample(seq_len(nrow(train)), replace = TRUE)
      X_train_boot <- X_train[boot_idx, ]
      y_train_boot <- y_train[boot_idx]
      boot_model <- glmnet(X_train_boot, y_train_boot, alpha = 1, lambda = best_lambda, family = "gaussian")
      boot_preds[i, ] <- as.numeric(predict(boot_model, s = best_lambda, newx = X_test))
    }

    val$q0.025 <- apply(boot_preds, 2, quantile, probs = 0.025)
    val$q0.975 <- apply(boot_preds, 2, quantile, probs = 0.975)
        
    val <- val %>%
      dplyr::select(site_code, s_index, Week_date, week, avg_Log10_NoV, predicted, q0.025, q0.975) %>%
      st_drop_geometry()
      
    predictions <- rbind(predictions, val)

    # Calculate metrics
    mse <- mean((val$avg_Log10_NoV - val$predicted)^2, na.rm = TRUE)
    rmse <- sqrt(mse)
    mae <- mean(abs(val$avg_Log10_NoV - val$predicted), na.rm = TRUE)
    mape <- mean(abs((val$avg_Log10_NoV - val$predicted) / val$avg_Log10_NoV), na.rm = TRUE) * 100
    bias <- mean(val$predicted - val$avg_Log10_NoV, na.rm = TRUE)
    pbias <- (bias / mean(val$avg_Log10_NoV, na.rm = TRUE)) * 100
    corr <- as.numeric(cor(val$avg_Log10_NoV, val$predicted, use="complete.obs", method="spearman"))  
    cov <- COV(val$avg_Log10_NoV, lower = val$q0.025, upper = val$q0.975)
    
    metrics_fold <- data.frame(
    Fold = k,
    MSE = mse,
    RMSE = rmse,
    MAE = mae,
    MAPE = mape,
    BIAS = bias,
    pBIAS = pbias,
    CORR = corr,
    COV = cov
  )
    metrics <- rbind(metrics, metrics_fold)

    print(paste(k, "folds done"))
   
}

# 95 % coverage probability
# Coverage
coverage_vector <- numeric(nrow(predictions))

for (i in 1:nrow(predictions)) {
  # Get the confidence interval bounds for the i-th observation
  lower_bound <- predictions$q0.025[i]
  upper_bound <- predictions$q0.975[i]
  
  # Get the observed value for the i-th observation
  observed <- predictions$avg_Log10_NoV[i]
  
  # Check if the observed value is within the confidence interval bounds
  coverage <- ifelse(observed >= lower_bound & observed <= upper_bound, 1, 0)
  
  # Store the coverage result (1 if within interval, 0 if outside) for the i-th observation
  coverage_vector[i] <- coverage
}

# Calculate the overall coverage probability
coverage_probability <- mean(coverage_vector, na.rm = TRUE) * 100
print(paste("Coverage Probability (%):", coverage_probability))



summary_metrics <- metrics %>%
  summarise(
    Mean_MSE = mean(MSE, na.rm = TRUE),
    SE = sd(MSE, na.rm = TRUE) / sqrt(length(MSE)),
    CI_Lower = Mean_MSE - qt(0.975, df = length(MSE) - 1) * SE,
    CI_Upper = Mean_MSE + qt(0.975, df = length(MSE) - 1) * SE,
    Mean_rmse = mean(RMSE, na.rm = TRUE),
    SE_rmse = sd(RMSE, na.rm = TRUE) / sqrt(length(RMSE)),
    CI_Lower_rmse = Mean_rmse - qt(0.975, df = length(RMSE) - 1) * SE_rmse,
    CI_Upper_rmse = Mean_rmse + qt(0.975, df = length(RMSE) - 1) * SE_rmse,
    Mean_MAE = mean(MAE, na.rm = TRUE),
    Mean_MAPE = mean(MAPE, na.rm = TRUE),
    Mean_BIAS = mean(BIAS, na.rm = TRUE),
    Mean_pBIAS = mean(pBIAS, na.rm = TRUE),
    Mean_CORR = mean(CORR, na.rm = TRUE),
    COV = coverage_probability
  )
  

# Print the summary of metrics
print(summary_metrics)


write.csv(metrics, "../../output/ML_lasso_model_metrics_full_final.csv")
write.csv(predictions, "../../output/ML_lasso_model_predictions_full_final.csv")
