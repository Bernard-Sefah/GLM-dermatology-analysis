library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(mi)
library(easyalluvial)
library(viridis)
library(janitor)


# na.strings converts the "?" in the age column to NA at read time
raw_data <- read.csv("dermatology.csv", na.strings = c("?", "NA", ""))

score_labels <- c("Low severity", "Higher severity")

collapse_score <- function(x) {
  factor(ifelse(x <= 1, score_labels[1], score_labels[2]),
         levels = score_labels)
}

clinical_vars  <- c("erythema", "scaling", "definite_borders", "itching")
disease_labels <- c("Psoriasis",
                    "Seborrheic dermatitis",
                    "Lichen planus",
                    "Pityriasis rosea",
                    "Chronic dermatitis",
                    "Pityriasis rubra pilaris")


# Check raw clinical scores are valid BEFORE grouping them
for (v in clinical_vars) {
  bad <- !is.na(raw_data[[v]]) & !(raw_data[[v]] %in% 0:3)
  if (any(bad)) {
    stop("Unexpected value(s) in ", v, "; expected scores 0, 1, 2, or 3.")
  }
}

# Check disease codes BEFORE converting to a factor
bad_class <- is.na(raw_data$class) | !(raw_data$class %in% 1:6)
if (any(bad_class)) {
  stop("Unexpected or missing disease class code(s) were found. Check the class column.")
}

skin <- raw_data %>%
  dplyr::mutate(
    age = as.numeric(age),
    disease = factor(class, levels = 1:6, labels = disease_labels),
    dplyr::across(tidyselect::all_of(clinical_vars), collapse_score)) %>%
  dplyr::select(tidyselect::all_of(clinical_vars), age, disease)

str(skin)
summary(skin)

## checking for missing data 
colSums(is.na(skin))
round(colMeans(is.na(skin)) * 100, 2)

missing_summary <- skin %>%
  summarise(
    Total_Observations = n(),
    Missing_Age        = sum(is.na(age)),
    Percentage_Missing = round(mean(is.na(age)) * 100, 2)
  )
print(missing_summary)

skin %>%
  group_by(disease) %>%
  summarise(
    Total              = n(),
    Missing_Age        = sum(is.na(age)),
    Percentage_Missing = round(mean(is.na(age)) * 100, 2),
    .groups = "drop"
  )

skin_mf <- missing_data.frame(as.data.frame(skin))
image(skin_mf)

## Final dataset 

analysis_data <- raw_data %>%
  dplyr::mutate(
    age = as.numeric(age),
    disease = factor(class, levels = 1:6, labels = disease_labels),
    dplyr::across(tidyselect::all_of(clinical_vars), collapse_score)) %>%
  dplyr::select(disease, age, tidyselect::all_of(clinical_vars))

model_vars <- c("disease", "age", clinical_vars)

analysis_complete <- analysis_data %>%
  filter(if_all(all_of(model_vars), ~ !is.na(.x))) %>%
  droplevels() %>%
  mutate(age10 = age / 10)   # age in decades: coefficients are per 10 years

# Make PRP the reference category
analysis_complete$disease <- relevel(analysis_complete$disease,
                                     ref = "Psoriasis")

cat("Original observations:", nrow(analysis_data), "\n")
cat("Complete observations used in primary models:", nrow(analysis_complete), "\n")
cat("Records excluded due to missing selected variables:",
    nrow(analysis_data) - nrow(analysis_complete), "\n")

# Which diseases lost patients to missing age?
analysis_data %>%
  count(disease, name = "n_total") %>%
  left_join(count(analysis_complete, disease, name = "n_used"), by = "disease") %>%
  mutate(n_used     = coalesce(n_used, 0L),
         n_excluded = n_total - n_used)

# Sparse-cell check: small cells can cause separation in a multinomial model
analysis_complete %>%
  pivot_longer(all_of(clinical_vars), names_to = "Predictor", values_to = "Level") %>%
  count(Predictor, Level, disease) %>%
  arrange(n) %>%
  head(15)



### EDA
p1 <- skin %>%
  dplyr::select(all_of(clinical_vars), disease) %>%
  alluvial_wide(fill_by = "last_variable")+ 
  labs(title = "Alluvial diagram, \n colored by sin disease bins")
p1

# Disease frequency
disease_frequency <- skin %>%
  count(disease) %>%
  mutate(Percentage = round(n / sum(n) * 100, 2))
print(disease_frequency)

ggplot(disease_frequency,
       aes(x = reorder(disease, -n), y = n, fill = disease)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  labs(title = "Frequency Distribution of Dermatological Diseases",
       x = "Disease Category", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Long format of the clinical predictors (full 0-3 scores)
skin_long <- skin %>%
  tidyr::pivot_longer(cols = tidyselect::all_of(clinical_vars),
    names_to = "Predictor",values_to = "Level")

# (categorical_frequency was used but never defined in the original script)
categorical_frequency <- skin_long %>%
  count(Predictor, Level, name = "n") %>%
  rename(Category = Level)


# Age boxplot
ggplot(skin, aes(y = age)) +
  geom_boxplot(fill = "red", na.rm = TRUE) +
  labs(title = "Boxplot of Patient Age", y = "Age (years)") +
  theme_minimal()

# Conditional density plots of diagnosis given age (bandwidth = 5 years)
p2 <- ggplot(skin, aes(x = age, fill = disease)) +
  geom_density(position = "fill", bw = 5) +
  geom_rug(aes(col = disease), alpha = 0.3) +
  labs(x = "Age (years)", y = "Proportion of patients",
       fill = "Disease", col = "Disease")
p2


# Proportional bar charts: clinical feature levels by disease
ggplot(skin_long, aes(x = disease, fill = Level)) +
  geom_bar(position = "fill", color = "white", linewidth = 0.3) +
  facet_wrap(~ Predictor) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = c("Low severity"      = "#A8DADC",
                               "Higher severity"  = "#D95D39")) +
  labs(title = "Proportional distribution of clinical features by disease",
       x = "Disease Category", y = "Proportion", fill = "Score group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

#### Statistical Modelling 

additive_formula <- disease ~age10 +erythema +scaling +definite_borders + itching

model_additive <- nnet::multinom(additive_formula,data = analysis_complete,
  trace = FALSE,maxit = 5000,Hess = TRUE,summ = 0)
summary(model_additive)


# Model coef and relative risk 

coef_mat <- coef(model_additive)

vcov_mat <- tryCatch(vcov(model_additive),
  error = function(e) NULL
)

if (is.null(vcov_mat)) {
  
  cat(
    "The covariance matrix could not be obtained. ",
    "Wald intervals are unavailable.\n"
  )
  
} else {
  
  se_matrix <- matrix(
    sqrt(diag(vcov_mat)),
    nrow = nrow(coef_mat),
    ncol = ncol(coef_mat),
    byrow = TRUE,
    dimnames = dimnames(coef_mat)
  )
  
  z_matrix <- coef_mat / se_matrix
  
  p_matrix <- 2 * pnorm(
    abs(z_matrix),
    lower.tail = FALSE
  )
  
  coefficient_table <- as.data.frame(
    as.table(coef_mat)
  )
  
  names(coefficient_table) <- c(
    "Disease_vs_reference",
    "Term",
    "Estimate"
  )
  
  coefficient_table$SE <- as.vector(se_matrix)
  
  coefficient_table$RRR <- exp(
    coefficient_table$Estimate
  )
  
  coefficient_table$CI_lower <- exp(
    coefficient_table$Estimate -
      1.96 * coefficient_table$SE
  )
  
  coefficient_table$CI_upper <- exp(
    coefficient_table$Estimate +
      1.96 * coefficient_table$SE
  )
  
  coefficient_table$z <- with(
    coefficient_table,
    Estimate / SE
  )
  
  coefficient_table$p_value <- 2 * pnorm(
    abs(coefficient_table$z),
    lower.tail = FALSE
  )
  
  coefficient_table <- coefficient_table %>%
    mutate(
      across(
        c(
          Estimate, SE, RRR,
          CI_lower, CI_upper,
          z, p_value
        ),
        ~ signif(.x, 4)))
  
  # Keep an unrounded numeric copy for any later calculations
  coefficient_table_raw <- coefficient_table
  
  # Display version: 4 decimals everywhere, formatted p-values
  coefficient_table <- coefficient_table %>%
    mutate(
      across(c(Estimate, SE, RRR, CI_lower, CI_upper, z),
             ~ sprintf("%.4f", .x)),
      p_value = ifelse(p_value < 0.0001, "<0.0001", sprintf("%.4f", p_value))
    )
  
  print(coefficient_table, row.names = FALSE)
}



# Null model and the additive model comparison

null_model <- nnet::multinom(disease ~ 1,data = analysis_complete,trace = FALSE,
  maxit = 5000,Hess = TRUE)

model_fit_table <- data.frame(
  
  Model = c(
    "Null Model",
    "Additive Multinomial Model"
  ),
  
  Log_likelihood = c(as.numeric(logLik(null_model)),
    as.numeric(logLik(model_additive))),
  
  Parameters = c(attr(logLik(null_model), "df"),
    attr(logLik(model_additive), "df")),
  
  Residual_deviance = c(deviance(null_model),
    deviance(model_additive)),
  
  AIC = c(AIC(null_model),AIC(model_additive)),
  
  BIC = c(BIC(null_model),BIC(model_additive)))

model_fit_table$McFadden_R2 <- c(NA_real_,
  1 - as.numeric(logLik(model_additive)) /
    as.numeric(logLik(null_model))
)

print(model_fit_table %>%mutate(across(
        where(is.numeric),
        ~ round(.x, 3) )))


### LRT

lrt_compare <- function(reduced, full, label) {
  
  ll0 <- as.numeric(logLik(reduced))
  ll1 <- as.numeric(logLik(full))
  
  df0 <- attr(logLik(reduced), "df")
  df1 <- attr(logLik(full), "df")
  
  statistic <- 2 * (ll1 - ll0)
  
  df_diff <- df1 - df0
  
  if (df_diff <= 0) {
    stop("The full model must have more parameters than the reduced model.")
  }
  
  data.frame(
    Predictor = label,
    LR_chisq = statistic,
    df = df_diff,
    p_value = pchisq(
      statistic,
      df = df_diff,
      lower.tail = FALSE
    )
  )
}


## LRT for the preditors 

predictor_terms <- c(
  "age10",
  clinical_vars
)

predictor_lrt_raw <- bind_rows(lapply(predictor_terms, function(term) {
  reduced_formula <- reformulate(
    setdiff(predictor_terms, term),
    response = "disease")
  
  reduced_fit <- nnet::multinom(reduced_formula,
                                data = analysis_complete, trace = FALSE, maxit = 5000)
  
  lrt_compare(reduced_fit, model_additive, term)
}))

# Display version: 4 decimals, formatted p-values
predictor_lrt <- predictor_lrt_raw %>%
  mutate(
    LR_chisq = sprintf("%.4f", LR_chisq),
    p_value  = ifelse(p_value < 0.0001, "<0.0001", sprintf("%.4f", p_value))
  )

print(predictor_lrt, row.names = FALSE)

### Is age Linear, cheking to answer  the 4the research question

# Linear model
fit_lin <- nnet::multinom(disease ~ age + erythema + scaling + definite_borders + itching,
                          data = analysis_complete, maxit = 5000, trace = FALSE)

# quadratic
fit_quad <- nnet::multinom(disease ~ age + I(age^2) + erythema + scaling +
                             definite_borders + itching,
                           data = analysis_complete, maxit = 5000, trace = FALSE)
 
# Likelihood-ratio test
anova(fit_lin, fit_quad)

# AIC and BIC
AIC(fit_lin, fit_quad)
BIC(fit_lin, fit_quad)


##### Assessment of  linearity in age assumption

# Compare linear age with a natural spline.

model_age_spline <- nnet::multinom(disease ~splines::ns(age10, df = 3) +
    erythema +scaling +definite_borders +itching,data = analysis_complete,
  trace = FALSE,maxit = 5000)

age_linearity_test <- lrt_compare(model_additive,model_age_spline,
  "Natural spline age vs linear age")

print(age_linearity_test %>%
    mutate(across( where(is.numeric),~ signif(.x, 4))))


# Age by itching interaction model

model_interaction <- nnet::multinom(disease ~age10 * itching +erythema +
    scaling + definite_borders,
  data = analysis_complete,trace = FALSE, maxit = 5000,Hess = TRUE)

# Compare additive and interaction models

interaction_lrt <- lrt_compare(model_additive,model_interaction,"Age-by-itching interaction")

print(interaction_lrt %>%
    mutate(p_value = signif(p_value, 4)))

### model comparison

comparison_table <- data.frame(Model = c(
    "Additive Model",
    "Age-by-itching Interaction Model"),
  Log_likelihood = c(
    as.numeric(logLik(model_additive)),
    as.numeric(logLik(model_interaction)) ),
  
  Parameters = c(attr(logLik(model_additive), "df"),
    attr(logLik(model_interaction), "df")),
  
  AIC = c(AIC(model_additive),AIC(model_interaction)),
  
  BIC = c(BIC(model_additive),BIC(model_interaction)))

print(comparison_table %>%mutate(
      across(where(is.numeric),~ round(.x, 3))))


##### Assessing if the model has converged

convergence_table <- data.frame(
  
  Model = c("Additive Model","Interaction Model","Age Spline Model" ),
  
  Convergence_code = c(model_additive$convergence,model_interaction$convergence,
  model_age_spline$convergence))

print(convergence_table)



## observed vs fitted probabilities 

predicted_all <- predict(
  model_additive,
  newdata = analysis_complete,
  type = "probs"
)

predicted_all <- as.data.frame(
  predicted_all
)

calibration_data <- bind_rows(
  lapply(
    disease_labels,
    function(category) {
      
      tibble(
        disease = factor(
          category,
          levels = disease_labels
        ),
        
        observed = as.integer(
          analysis_complete$disease == category
        ),
        
        predicted = predicted_all[[category]]
      )
      
    }
  )
) %>%
  group_by(disease) %>%
  mutate(
    probability_group = ntile(
      predicted,
      10
    )
  ) %>%
  group_by(
    disease,
    probability_group
  ) %>%
  summarise(
    mean_predicted = mean(predicted),
    observed_rate = mean(observed),
    n = n(),
    .groups = "drop"
  )

ggplot(
  calibration_data,
  aes(
    x = mean_predicted,
    y = observed_rate
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = 2,
    color = "gray40"
  ) +
  geom_point(
    aes(size = n),
    alpha = 0.8
  ) +
  geom_line() +
  facet_wrap(
    ~ disease,
    scales = "free"
  ) +
  labs(
    title = "Observed and Predicted Probabilities by Disease",
    x = "Mean Predicted Probability",
    y = "Observed Proportion",
    size = "Group n"
  ) +
  theme_minimal()


#### Classification perfomance 

# Stratified 75% training and 25% testing split.

test_prop <- 0.25

set.seed(123)

train_indices <- unlist(
  lapply(
    split(
      seq_len(nrow(analysis_complete)),
      analysis_complete$disease
    ),
    function(idx) {
      
      n_train <- max(
        1,
        floor((1 - test_prop) * length(idx))
      )
      
      sample(
        idx,
        size = n_train
      )
      
    }
  ),
  use.names = FALSE
)

train_indices <- sort(train_indices)

train_data <- analysis_complete[
  train_indices,
]

test_data <- analysis_complete[
  -train_indices,
]


### fitting model using training data 

train_model <- nnet::multinom(
  additive_formula,
  data = train_data,
  trace = FALSE,
  maxit = 5000
)

# Predict disease categories for test observations.

predicted_class <- predict(
  train_model,
  newdata = test_data,
  type = "class"
)

predicted_class <- factor(
  predicted_class,
  levels = levels(analysis_complete$disease)
)

observed_class <- factor(
  test_data$disease,
  levels = levels(analysis_complete$disease)
)


### Confusion Matrix

confusion <- table(
  Observed = observed_class,
  Predicted = predicted_class
)

print(confusion)


### Classification Metrics 

metric_by_class <- bind_rows(
  lapply(
    disease_labels,
    function(category) {
      
      tp <- sum(
        observed_class == category &
          predicted_class == category
      )
      
      fp <- sum(
        observed_class != category &
          predicted_class == category
      )
      
      fn <- sum(
        observed_class == category &
          predicted_class != category
      )
      
      support <- sum(
        observed_class == category
      )
      
      precision <- ifelse(
        tp + fp == 0,
        NA_real_,
        tp / (tp + fp)
      )
      
      recall <- ifelse(
        tp + fn == 0,
        NA_real_,
        tp / (tp + fn)
      )
      
      f1 <- ifelse(
        is.na(precision) ||
          is.na(recall) ||
          precision + recall == 0,
        NA_real_,
        2 * precision * recall /
          (precision + recall)
      )
      
      tibble(
        Disease = category,
        Support = support,
        Precision = precision,
        Recall = recall,
        F1 = f1
      )
      
    }
  )
)

# Overall accuracy.

accuracy <- mean(
  predicted_class == observed_class
)

# Balanced accuracy.

balanced_accuracy <- mean(
  metric_by_class$Recall,
  na.rm = TRUE
)

# Majority-class baseline.

majority_baseline <- max(
  prop.table(
    table(train_data$disease)
  )
)

performance_summary <- data.frame(
  
  Metric = c(
    "Test Accuracy",
    "Balanced Accuracy",
    "Training Majority-Class Baseline"
  ),
  
  Value = c(
    accuracy,
    balanced_accuracy,
    as.numeric(majority_baseline)
  )
)

print(
  performance_summary %>%
    mutate(
      Value = round(Value, 3)
    )
)

print(
  metric_by_class %>%
    mutate(
      across(
        c(Precision, Recall, F1),
        ~ round(.x, 3)
      )
    )
)


#### Multiclass ROC analysis

# One-versus-rest ROC analysis.

test_prob <- predict(
  train_model,
  newdata = test_data,
  type = "probs"
)

test_prob <- as.data.frame(
  test_prob
)

roc_results <- lapply(
  disease_labels,
  function(category) {
    
    if (!category %in% names(test_prob)) {
      stop(
        "Predicted-probability column not found for: ",
        category
      )
    }
    
    # Check that both outcome classes are represented.
    
    binary_outcome <- factor(
      ifelse(
        observed_class == category,
        "Case",
        "Control"
      ),
      levels = c(
        "Control",
        "Case"
      )
    )
    
    if (length(unique(binary_outcome)) < 2) {
      stop(
        "Both cases and controls are required for ROC analysis: ",
        category
      )
    }
    
    roc_fit <- pROC::roc(
      binary_outcome,
      test_prob[[category]],
      levels = c("Control", "Case"),
      direction = "<",
      quiet = TRUE
    )
    
    data.frame(
      Disease = category,
      Specificity = roc_fit$specificities,
      Sensitivity = roc_fit$sensitivities,
      AUC = as.numeric(
        pROC::auc(roc_fit)
      )
    )
    
  }
)

roc_data <- bind_rows(
  roc_results
)

# AUC for each disease.

auc_table <- roc_data %>%
  distinct(
    Disease,
    AUC
  ) %>%
  arrange(
    desc(AUC)
  )

print(
  auc_table %>%
    mutate(
      AUC = round(AUC, 3)
    )
)

# Macro-average AUC.

macro_auc <- mean(
  auc_table$AUC,
  na.rm = TRUE
)

print(
  data.frame(
    Macro_average_AUC = round(macro_auc, 3)
  )
)


#### ROC curves 

ggplot(roc_data,
  aes(x = 1 - Specificity,y = Sensitivity,
    colour = Disease)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1,intercept = 0,linetype = 2,
    colour = "gray50"
  ) +
  coord_equal() +
  labs(
    title = "One-versus-Rest ROC Curves",
    x = "False-Positive Rate (1 - Specificity)",
    y = "True-Positive Rate (Sensitivity)",
    colour = "Disease Category"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom"
  )









