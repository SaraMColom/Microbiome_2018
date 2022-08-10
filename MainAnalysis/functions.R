# Additional functions

# Code to further tidy results table

tidy_more <- function(x){
  x %>% 
  rename(Term = term,
         SS = sumsq,
         DF = df,
         `F-value` = statistic,
  ) %>% 
    mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
    mutate(P = `p.value`) %>% 
    mutate(P = case_when(`p.value` < 0.001 ~ "<0.001***",
                         `p.value` < 0.01 ~ paste(`p.value`, "**", sep = ""),
                         `p.value` < 0.05 ~ paste(`p.value`, "*", sep = ""),
                         TRUE ~ P %>% as.character())) 
}

tidy_more2 <- function(x){
  x %>% 
    rename(Term = term,
          SE = `std.error`,
          `F-value` = statistic,
    ) %>% 
    mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
    mutate(P = `p.value`) %>% 
    mutate(P = case_when(`p.value` < 0.001 ~ paste("<0.001", "***"),
                         `p.value` < 0.01 ~ paste(`p.value`, "**"),
                         `p.value` < 0.05 ~ paste(`p.value`, "*"),
                         TRUE ~ P %>% as.character())) 
}


tidy_p <- function(x){
  x %>%
    mutate_if(is.numeric, function(x) round(x, 3)) %>% 
    mutate(`p-value` = case_when(p_value < 0.001 ~ paste("<0.001***"),
                                 p_value < 0.01 ~ paste(p_value, "**", sep = ""),
                                 p_value < 0.05 ~ paste(p_value, "*", sep = ""),
                                 p_value > 0.99 ~ ">0.99",
                                 TRUE ~ paste(p_value)))
}