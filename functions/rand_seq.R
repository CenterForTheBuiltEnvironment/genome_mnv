require(base)

# Function to get next decision
get_next_decision <- function(last_decision, consecutive_count, interval) {
  if (consecutive_count >= interval * 3) {
    # Switch decision after 9 consecutive occurrences
    return(ifelse(last_decision == 1, 2, 1))
  } else {
    # Randomly choose between A and B, but ensure it's different after three timesteps
    if (consecutive_count %% interval == 0) {
      return(sample(c(1, 2), 1))
    } else {
      return(last_decision)
    }
  }
}


# Function to generate the schedule
generate <- function(start_date, end_date, interval){
  
  # Get date
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  # Generate the sequence of dates
  dates <- seq(from = start_date, to = end_date, by = "day")
  
  # Initialize the decisions vector
  decisions <- rep(NA, length(dates))
  
  # Initialize counters
  last_decision <- sample(c(1, 2), 1)  # Start with a random decision
  consecutive_count <- 0
  
  # Fill the decisions
  for (i in 1:length(dates)) {
    decisions[i] <- get_next_decision(last_decision, consecutive_count, interval)
    
    # Update last decision and count
    if (decisions[i] == last_decision) {
      consecutive_count <- consecutive_count + 1
    } else {
      last_decision <- decisions[i]
      consecutive_count <- 1
    }
  }
  
  # Combine dates and decisions into a data frame
  decision_data <- data.frame(date = dates, strategy = decisions)
  
  # Print the first few rows to check
  return(decision_data)
  
}

# Functions to evaluate the schedule
evaluate <- function(decision_data, exclude){
  
  wday_df <- decision_data %>% 
    mutate(flag = strategy,
           switch = ifelse(flag == lag(flag, 1), 0, flag),
           consec = ifelse(flag == lag(flag, 1), flag, 0)) %>%
    mutate(type = ifelse(switch != 0, "non-consec", "consec")) %>%
    mutate(type = replace(type, 1, "start")) %>% 
    mutate(type = as.factor(type)) %>% 
    mutate(strategy = as.factor(strategy)) %>% 
    filter(type != "non-consec") %>% 
    filter(!date %in% exclude$date) %>% 
    group_by(weekday = wday(date, week_start = 1),
             strategy) %>%
    count() %>% 
    ungroup()
  
  margin_wday <- max(wday_df$n) - min(wday_df$n)
  
  strategy_df <- decision_data %>% 
    mutate(flag = strategy,
           switch = ifelse(flag == lag(flag, 1), 0, flag),
           consec = ifelse(flag == lag(flag, 1), flag, 0)) %>%
    mutate(type = ifelse(switch != 0, "non-consec", "consec")) %>%
    mutate(type = replace(type, 1, "start")) %>% 
    mutate(type = as.factor(type)) %>% 
    mutate(strategy = as.factor(strategy)) %>% 
    filter(type != "non-consec") %>% 
    filter(!date %in% exclude$date) %>% 
    group_by(strategy) %>%
    count() %>% 
    ungroup()
  
  margin_strategy <- abs(diff(strategy_df %>% .$n))
  
  return(list(margin_wday = margin_wday, 
              margin_strategy = margin_strategy,
              wday_df = wday_df))
    
}

# Function to formulate the process
rand_seq <- function(start_date, end_date, interval, threshold, exclude){
  regenerate <- T
  i <- 0
  success <- F
  while (regenerate & i <= 5000){
    decision_data <- generate(start_date, end_date, interval)
    eva <- evaluate(decision_data, exclude)
    if (eva$margin_wday > threshold | eva$margin_strategy > threshold){
      i <- i + 1
    } else {
      regenerate <- F
      success <- T
      print("success!")
    }
  }
  
  if (!success){
    print("try again!")
    threshold_loose <- threshold + 1
    regenerate <- T
    i <- 0
    while (regenerate & i <= 5000){
      decision_data <- generate(start_date, end_date, interval)
      eva <- evaluate(decision_data, exclude)
      if (eva$margin_wday > threshold_loose | eva$margin_strategy > threshold_loose){
        i <- i + 1
      } else {
        regenerate <- F
        success <- T
        print("success in loose!")
      }
    }
  }
  
  if (!success){
    print("fail!")
  }
  
  return(list(decision_data = decision_data, eva = eva))
}
