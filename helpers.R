library(policytree)
library(labelled)
library(fuzzyjoin)


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

basu_forest <- function(X, Y, W, clusters, num.trees, iters=1, ret_list = FALSE) {
  W <- as.factor(W)
  for_list <- list()
  
  for(i in 1:iters) {
    forest <- multi_arm_causal_forest (
      X, 
      Y, 
      W,
      num.trees = NUM_TREES,
      clusters = clusters
    )
    
    X <- X %>%
      select(
        which(
          variable_importance(forest) > mean(variable_importance(forest))
          )
        )
    for_list[[i]] <- forest
  }
  
  if (ret_list) return(for_list)
  for_list[[iters]]
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

