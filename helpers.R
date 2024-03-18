library(policytree)
library(labelled)
library(fuzzyjoin)
# library(MASS) # For linear regression lm()
# library(glmnet) # For LASSO regression
library(rpart)
library(rpart.plot)
library(rsample)
library(boot)


impute_covars <- function(X) {
  X_imp <- X %>%
    missMethods::impute_median()
  
  X_imp[, colSums(is.na(X_imp)) < nrow(X_imp), drop = FALSE]
}


trim_distribution <- function(x, pct) {
  quants <- ntile(x, 100)
  map_lgl(quants, ~.x >= pct/2 & .x <= (100 - (pct/2)))
}


create_cowplot <- function(plot_data, variable_list, tau1, pctile_trim = 4, var_labels = NULL) {
  if(!is.data.frame(tau1)) tau1 <- as.data.frame(tau1)
  names(tau1) <- c("Group 2", "Group 3", "Group 4")
  
  if(is.null(var_labels)) var_labels <- variable_list
  var_labels <- str_wrap(var_labels, 25)
  
  num_treats <- ncol(tau1)
  plots <- vector(mode = "list", length = num_treats * length(variable_list))
  
  # Use an inner and outer map to map over variables then treatments within variables
  new_plots <- map2(
    variable_list, var_labels,
    function(variable, label) {
      # Set number of options in question for later use
      num_options <- length(unique(pull(plot_data[variable])))
      
      # Trim top an bottom X values
      if (num_options > 10) {
        trim_selector <- trim_distribution(plot_data[variable] %>% pull(),pctile_trim)
        data1 <- plot_data[trim_selector,]
        tau2 <- tau1[trim_selector,]
      } else {
        data1 <- plot_data
        tau2 <- tau1
      }
      
      map2(tau2, names(tau2), function(tau, name) {
        if (num_options > 10) {
          ggplot(data = data1, aes(x=!!rlang::sym(variable), y = tau)) +
            geom_jitter(alpha = 0.5) +
            geom_smooth(color = "red") +
            geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
            labs(x = label, y = "Treatment effect", subtitle = name) +
            theme(text = element_text(size = 9))
          
        } else {
          ggplot(data = data1, aes(x= !!rlang::sym(variable), y = tau)) +
            geom_jitter(alpha = 0.5) +
            geom_violin(aes(group = !!rlang::sym(variable)), alpha = 0.25) +
            stat_summary(fun = "mean",
                         geom = "point",
                         color = "red",
                         size = 3) +
            geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
            labs(x = label, y = "Treatment effect", subtitle = name) +
            theme(text = element_text(size = 9))
        }
      })
    }
  )
  

  plot_grid(plotlist = new_plots %>% unlist(recursive = F), 
            ncol = num_treats)
}


qinish <- function(forest, treat_costs, learning_vars = NULL) {
  # Get net benefit (benefit minus cost)
  treat_tau <- forest$predictions %>% as.data.frame()
  
  dfs <- map(treat_tau, function(tau) {
    if (!is.vector(treat_tau)) treat_tau <- pull(treat_tau)
    if (length(treat_costs) == 1) treat_costs <- rep.int(treat_costs, length(treat_tau))
    cb_ratio = treat_costs / treat_tau
    cb_ratio[cb_ratio < 0] <- Inf
    df <- data.frame(cb_point = cb_ratio, port_treat = 0)
    for (i in 1:length(cb_ratio)) {
      selector <- cb_ratio[i] <= cb_ratio
      df[i, "port_treat"] <- sum(selector) / length(cb_ratio)
    }
    df
  })
  
  names(dfs) <- paste0("treat", ncol(treat_tau))
  
  first_degree <- bind_rows(dfs, .id = "treatment")
  
  # Look at policy learner
  
  # Get covariates
  if (is.null(learning_vars)) {
    learning_vars <- forest$X.orig
  }
  
  # Get DR scores
  dr_scores <- double_robust_scores(forest)
  
  # Rerun analysis for learnt policy
  first_degree %>%
    mutate(
      pct_treated = hybrid_policy_tree(
        learning_vars %>% impute_covars(),
        dr_scores[,c(1, as.numeric(treatment) + 1)],
        search.depth = 1
      ) %>%
        predict(learning_vars %>% impute_covars()) %>%
        as.tibble() %>%
        mutate(dr = dr_scores[,c(1, as.numeric(treatment) + 1)]) %>%
        filter(value == 2) %>%
        nrow()
        ,
      pct_treated1 = pct_treated / nrow(treat_tau)
    )
}



# Join maths to baseline based on names
join_on_names <- function(baseline, endline, maths) {

  maths$hhid_endline <- maths$hhid_endline %>%
    as_factor() %>%
    as.character()
  maths %>%
    filter(!is.na(hhid_endline) & hhid_endline != "")
  
  endline$hhid_endline <- endline$hhid_endline %>%
    as_factor() %>%
    as.character()
  
  baseline$endline
  
  out_df <- endline %>%
    filter(!is.na(hhid_endline) & hhid_endline != "") %>%
    labelled::look_for_and_select("Child's first name") %>%
    bind_cols(endline["hhid_endline"]) %>%
    gather(key = "Variable", value = "Name", -hhid_endline) %>%
    left_join(maths, by = "hhid_endline") %>%
    group_by(hhid_endline) %>%
    mutate(str_dist = stringdist::stringdist(prenom_enf_test, Name, method = "lv") / nchar(prenom_enf_test)) %>%
    arrange(str_dist) %>%
    slice(1) %>%
    ungroup()
    
  out_df_filtered <- out_df %>%  
    filter(str_dist <= 0.25)
  
  if (nrow(out_df) != nrow(out_df_filtered)) {
    warning(paste0("Some rows did not have a valid match at threshold of 0.25. ",
                   nrow(out_df) - nrow(out_df_filtered),
                   " rows did not have a match."))
  }
  
  out_df_filtered
  
}


get_var_labels <- function(orig_data, variables) {
  data_names <- names(orig_data)
  print("Getting labels, this could take a little while")
  cond_1 <- variables %in% data_names
  cond_2 <- variables %>% str_replace("_self", "_1") %in% data_names
  cond_3 <- (variables %>% str_detect("^c_")) &
    (variables %>% 
      str_replace("_self", "") %>% 
      str_replace("^c_", "c1_")
    ) %in% data_names
  
  input_list <- list(variable = variables, cond1 = cond_1, cond2 = cond_2, cond3 = cond_3)
  
  pmap_chr(input_list, 
    function(variable, cond1, cond2, cond3) {
      # If variable is already labelled, just get existing label
      if (cond1) {
        orig_data %>% select(variable) %>% var_label() %>% unlist()
      } else if (cond2) {
        orig_data %>% 
          select(variable %>% str_replace("_self", "_1")) %>% 
          var_label() %>% 
          unlist()
      } else if (cond3) {
        orig_data %>% 
          select(variable %>% 
                    str_replace("_self", "") %>% 
                    str_replace("^c_", "c1_")) %>%
          var_label()%>% 
          unlist()
      } else {
        NA
      }
  }, .progress = TRUE)
}


get_vars_by_regex <- function(df, selector, num_selector, person_num, suffix = "_self") {
  df[paste0(person_num, "join")] <- df[person_num] %>% 
    pull() %>%
    as.numeric()
  
  df %>%
    select(hhid, matches(selector)) %>%
    gather("name", "value", -hhid)%>%
    transmute(value = value, 
              person = str_extract(name, num_selector) %>% str_extract("[:digit:]+") %>% as.numeric(),
              variable = str_replace(name, num_selector, ""),
              variable = paste0(variable, suffix),
              hhid = hhid
              ) %>%
    pivot_wider(names_from = variable, values_from = value) %>%
    right_join(df, by = c("hhid", "person" = paste0(person_num, "join")))
}

basu_forest <- function(X, Y, W, clusters, num.trees, iters=Inf, ret_list = FALSE, num_vars = 0) {
  for_list <- list()
  What <- regression_forest(X, W, num.trees) %>% .$predictions %>% as.numeric()
  Yhat <- regression_forest(X, Y, num.trees) %>% .$predictions %>% as.numeric()
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = iters, style = 3)
  i <- 1
  
  while(ncol(X) > num_vars & i < iters + 1) {
    # Update progress bar
    setTxtProgressBar(pb, i)
    
    if (length(unique(W)) > 2) {
      W <- as.factor(W)
      
      forest <- multi_arm_causal_forest (
        X, 
        Y, 
        W,
        num.trees = num.trees, # Ensure this uses the function argument 'num.trees' correctly
        clusters = clusters,
        W.hat = What,
        Y.hat = Yhat
      )
    }
    else {
      forest <- causal_forest (
        X, 
        Y, 
        W,
        num.trees = num.trees, # Ensure this uses the function argument 'num.trees' correctly
        clusters = clusters,
        W.hat = What,
        Y.hat = Yhat
      )
    }
    
    # Assuming you are using dplyr for 'select' and 'which', make sure the library is loaded
    X <- X %>%
      dplyr::select(
        which(
          variable_importance(forest) > mean(variable_importance(forest))
        )
      )
    for_list[[i]] <- forest
    
    i <- i + 1
  }
  
  # Close progress bar
  close(pb)
  
  if (ret_list) return(for_list)
  for_list[[iters]]
  
  print(paste0("Finished in ", i, " iterations."))
  
  for_list
}



display_variable_imp <- function(forest, labelled_data) {
  all_labels <- labelled_data %>% 
    labelled::get_variable_labels() %>%
    unlist()
  
  test <- get_var_imp(forest) %>%
    arrange(desc(Importance))
  
  test$Label <- get_var_labels(labelled_data, test$Variable)
  
  test %>% select(Label, Variable, Importance)
}

get_var_imp <- function(forest) {
  varimp <- variable_importance(forest) %>% as.data.frame.matrix()
  var_names <- names(forest$X.orig)
  
  varimp <- bind_cols(var_names, format(varimp %>% pull() %>% round(3), nsmall = 3))
  names(varimp) <- c("Variable", "Importance")
  
  varimp
}


median_split <- function(x) {
  if (x %>% na.omit() %>% unique() %>% length() == 2) return(x)
  
  median_value <- median(x, na.rm = TRUE)
  ifelse(x >= median_value, 1, 0)
}


graph_single_hte <-function(variable1, name, forest, trim = NULL) {
  num_options <- variable1 %>% unique() %>% length()
  tau1 <- forest$predictions %>% as.numeric()
  dr_scores <- forest %>% get_scores()
  data1 <- data.frame(variable = variable1, tau = tau1, dr_scores = forest %>% get_scores())
  
  if (!is.null(trim)) {
    lower_crit <- quantile(morocco$h4_1, trim / 2, na.rm = T)
    upper_crit <- quantile(morocco$h4_1, 1 - (trim / 2), na.rm = T)
    data1 <- data1 %>% filter(variable <= upper_crit, variable >= lower_crit)
  }
  
  if (num_options > 10) {
    ggplot(data = data1, aes(x = variable, y = tau)) +
      geom_jitter(alpha = 0.5) +
      geom_smooth(aes(y = tau), color = "red", se = FALSE) +
      geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
      labs(x = "", y = "Treatment effect", subtitle = name) +
      theme(text = element_text(size = 9))
    
  } else {
    ggplot(data = data1, aes(x = variable, y = tau)) +
      geom_jitter(alpha = 0.5) +
      geom_violin(aes(group = variable), alpha = 0.25) +
      stat_summary(fun = "mean",
                   geom = "point",
                   color = "red",
                   size = 3) +
      geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
      labs(x = "", y = "Treatment effect", subtitle = name) +
      theme(text = element_text(size = 9))
  }
}

get_covariate_variance <- function(forest) {
  forest$X.orig %>% 
    map_dbl(~.x %>% var(na.rm = T))
}



# Define the function for estimating CATE
dr_adaptive_lasso <- function(D, covariates, gamma = 1, lambda = 0.01) {
  # Step 1 and 2 are assumed to be done outside this function and D is the doubly robust score
  
  # Step 4(a): Run a linear regression of D on V (covariates)
  lm_model <- lm(D ~ ., data = covariates)
  betas <- coef(lm_model)[-1] # Exclude the intercept
  
  # Step 4(b): Define the weights
  weights <- 1 / (abs(betas) ^ gamma)
  
  # Step 4(c): Run a LASSO regression of D on V with weights as the penalty factor
  lasso_model <- glmnet(as.matrix(covariates), D, weights = weights, lambda = lambda, alpha = 1, standardize = FALSE)
  
  # Determine the lambda value to use (for simplicity we're using the one supplied)
  # In practice, cross-validation should be used to find an optimal value for lambda
  # This would be done with cv.glmnet()
  
  # Step 4(d): Extract non-zero coefficients (the selected effect modifiers)
  lasso_coefficients <- predict(lasso_model, type = "coefficients", s = lambda)[1:nrow(covariates), , drop = FALSE]
  nonzero_indices <- which(lasso_coefficients != 0)
  
  # Step 5: Calculate the final estimate of the CATE for each subject
  intercept <- lasso_coefficients[1, 1] # Extract the intercept from the LASSO model
  cate_estimates <- intercept + as.matrix(covariates) %*% lasso_coefficients[nonzero_indices, , drop = FALSE]
  
  return(list(cate_estimates = cate_estimates, selected_modifiers = nonzero_indices))
}


predict_nodes <-
  function (object, newdata, na.action = na.pass) {
    where <-
      if (missing(newdata)) 
        object$where
    else {
      if (is.null(attr(newdata, "terms"))) {
        Terms <- delete.response(object$terms)
        newdata <- model.frame(Terms, newdata, na.action = na.action, 
                               xlev = attr(object, "xlevels"))
        if (!is.null(cl <- attr(Terms, "dataClasses"))) 
          .checkMFClasses(cl, newdata, TRUE)
      }
      rpart:::pred.rpart(object, rpart:::rpart.matrix(newdata))
    }
    as.integer(row.names(object$frame))[where]
  }


distilled_causal_tree <- function(forest, maxdepth = 3, seed = 123, replicates = 500, num_candidate_trees = 100) {
  if (multi_arm <- "multi_arm_causal_forest" %in% class(forest)) {
    distilled_causal_tree_multi(forest, maxdepth, seed, replicates, num_candidate_trees)
  } else {
    distilled_causal_tree_binary(forest, maxdepth, seed, replicates, num_candidate_trees)
  }
}

distilled_causal_tree_multi <- function(forest, maxdepth, seed, replicates, num_candidate_trees) {
  unique_treat <- forest$W.orig %>% as.numeric() %>% unique()
  unique_treat <- unique_treat[-1]
  dr_scores <- get_scores(forest) %>% as.data.frame()
  preds <- forest$predictions %>% as.data.frame()
  map2(preds, dr_scores, ~distilled_causal_tree_binary(forest, maxdepth, seed, replicates, preds = .x, dr_scores = .y, num_candidate_trees = num_candidate_trees))
}



distilled_causal_tree_binary <- function(forest, maxdepth, seed, replicates, num_candidate_trees, preds = forest$predictions, dr_scores = get_scores(forest)) {
  # Pull from forest
  X <- forest$X.orig
  
  
  # Set up DF
  data <- data.frame(y = preds) %>%
    bind_cols(X)
  
  # Partition data
  # Assuming 'data' is your dataset and you want a 70-30 split
  set.seed(seed) # For reproducibility
  total_rows <- nrow(data)
  shuffled_indices <- sample(1:total_rows)
  
  # Calculate the index to split on for a 70-30 split
  split_index <- round(total_rows * 0.5)
  
  # Create the training and testing datasets
  train_data <- data[shuffled_indices[1:split_index], ]
  test_data <- data[shuffled_indices[(split_index + 1):total_rows], ]
  
  # Fit tree
  fit <- fit_best_tree(train_data, formula = y ~ ., num_candidate_trees, seed = seed, maxdepth = maxdepth)
  fit$fit_data <- shuffled_indices[1:split_index]
  fit$est_data <- shuffled_indices[(split_index + 1):total_rows]
  
  out <- get_honest_estimates(fit, data, dr_scores, replicates)
  
  list(estimates = out, model = fit)
}



get_honest_estimates <- function(tree, data = NULL, dr_scores = NULL, replicates, forest = NULL) {
  if(!is.null(forest)) {
    dr_scores <- get_scores(forest)
    data <- forest$X.orig
  }
  
  dr_scores <- dr_scores[tree$est_data]
  data <- data[tree$est_data,]
  
  
  # Get estimates for each leaf
  held_out_nodes <- findDataPointsInNodes(tree, data = data)
  
  dr_point_est <- data.frame(node = names(held_out_nodes))
  dr_point_est$est <- map_dbl(dr_point_est$node, ~mean(dr_scores[held_out_nodes[[.x]]], na.rm = T))
  dr_point_est$n <- map_int(dr_point_est$node, ~length(held_out_nodes[[.x]]))
  
  # Bootstrap confidence intervals
  base_bs <- data.frame(
    dr_scores = dr_scores
    ) %>%
    bind_cols(expand_node_encoding(held_out_nodes))
  
  bs <- rsample::bootstraps(base_bs, times = replicates)
  
  out <- bs$splits %>%
    map(function (x) {
      dr_scores_bs <- analysis(x)
      
      bs_df <- data.frame(node = names(held_out_nodes))
      bs_df$est <- map_dbl(dr_point_est$node, ~mean(dr_scores_bs[dr_scores_bs[,.x],"dr_scores"], na.rm = T))
      bs_df
    }) 
  
  out <- out %>%
    bind_rows() %>%
    group_by(node) %>%
    summarise(se = sd(est)) %>%
    left_join(dr_point_est, by = "node")
  
  paths <- path.rpart(tree, as.numeric(names(held_out_nodes))) %>% map(~paste(.x[-1], collapse = "    "))
  
  out <- out %>%
    mutate(
      sig = ifelse(abs(est) > 1.96 * se, "*", "")
      )
  
  paths <- data.frame(
    node = names(paths),
    rules = paths %>% unlist()
  ) %>%
    mutate(
      depth = str_count(rules, "<") + str_count(rules, ">")
    )
    
  out %>%
    left_join(paths, by = "node") %>%
    dplyr::select(node, depth, est, se, sig, n, rules) %>%
    arrange(as.numeric(node))
}




# Function to find best tree across bootstraps
fit_best_tree <- function(data, formula, B, maxdepth, seed = 123) {
  
  best_mse <- Inf
  
  # Loop over B subsamples
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (i in 1:B) {
    # Generating a subsample
    bootstrap_indices <- sample(round(nrow(data) / 2), replace = FALSE)
    bootstrap_data <- data[bootstrap_indices, ]
    
    # Fitting a tree to the subsample
    model <- rpart(formula = formula, data = bootstrap_data, method = "anova", maxdepth = maxdepth, usesurrogate = 0)
    
    # Assess performance on the left out samples
    # Calculating average MSE for the current model
    val_data <- data[-bootstrap_indices, ]
    predictions <- predict(model, val_data, type = "vector")
    avg_mse <- mean((val_data$y - predictions)^2)
    
    if (avg_mse < best_mse) best_model <- model
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  # Close progress bar
  close(pb)
  

  best_model
}




# Function to find data points in each node of an rpart tree
findDataPointsInNodes <- function(model, data) {
  paths <- path.rpart(model, row.names(model$frame))
  map(paths, ~get_df_indices_based_on_rules(data, .x[-1]))
}



get_df_indices_based_on_rules <- function(df, rules) {
  if (length(rules) == 0) {
    # Return all indices if there are no rules
    return(1:nrow(df))
  } else {
    indices <- 1:nrow(df)
    
    for (rule in rules) {
      # Parse the rule into column, operator, and value
      match <- regmatches(rule, gregexpr("<|<=|>|>=|==|!=", rule))[[1]]
      parts <- strsplit(rule, "<|<=|>|>=|==|!=")[[1]]
      column_name <- trimws(parts[1])
      operator <- match[1]
      value <- as.numeric(trimws(parts[2]))
      
      # Construct the filtering expression and get indices
      current_indices <- which(eval(parse(text = sprintf("df$`%s` %s %f", column_name, operator, value))))
      indices <- intersect(indices, current_indices)
    }
    
    return(indices)
  }
}


expand_node_encoding <- function(held_out_nodes) {
  list_df <- enframe(held_out_nodes, name = "item", value = "index") %>%
    unnest(index) %>%
    mutate(value = TRUE) # This column indicates presence
  
  # Create a full expansion of item by index
  full_df <- crossing(item = names(held_out_nodes), index = 1:length(unique(unlist(held_out_nodes))))
  
  # Left join and replace NA with FALSE
  final_df <- full_df %>%
    left_join(list_df, by = c("item", "index")) %>%
    replace_na(list(value = FALSE)) %>%
    pivot_wider(names_from = item, values_from = value, values_fill = list(value = FALSE))
}


node.fun.se <- function(x, labs, digits, varlen)
{
  split_labs <- labs %>% str_split("\n")
  s <- round(as.numeric(x$frame$se), 3) # round sd to 2 digits
  outside_info <- list(est = x$frame$yval, se = s, n = x$frame$n) %>% transpose()
  map2_chr(split_labs, outside_info, function(vec, info) {
    sig_star <- ifelse(abs(info$est) > 1.96 * info$se, "*", "")
    c(paste0(format(round(info$est, 3), nsmall = 3L), sig_star), paste0("(", format(round(info$se, 3), nsmall = 3L), ")"), info$n) %>%
      paste(collapse = "\n")
  })
}


plot_ddrct <- function(ddrct) {
  # Edit model yval w DR estimates
  results <- ddrct$estimates %>%
    # mutate(est = round(est, 3), se = round(se, 3)) %>%
    arrange(as.numeric(node))
  
  sorted_results <- results[as.numeric(row.names(ddrct$model$frame)),]
  ddrct$model$frame$yval <- sorted_results$est
  ddrct$model$frame$se <- sorted_results$se
  ddrct$model$frame$n <- sorted_results$n
  
  # Get model scaffold
  rpart.plot(ddrct$model, node.fun = node.fun.se, nn = TRUE, box.col  = "lightblue")
}
