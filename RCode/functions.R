# Additional functions

# Code to further tidy results table

tidy_more <- function(x){
  x %>% 
  rename(Term = term,
         SS = sumsq,
         DF = df,
         `t-value` = statistic,
  ) %>% 
    mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
    mutate(P = `p.value`) %>% 
    mutate(P = case_when(P < 0.001 ~ paste(P, "***"),
                         P < 0.01 ~ paste(P, "**"),
                         P < 0.05 ~ paste(P, "*"),
                         TRUE ~ P %>% as.character())) 
}

tidy_more2 <- function(x){
  x %>% 
    rename(Term = term,
          SE = `std.error`,
          `t-value` = statistic,
    ) %>% 
    mutate_if(is.numeric, function(x){round(x, 3)}) %>% 
    mutate(P = `p.value`) %>% 
    mutate(P = case_when(P < 0.001 ~ paste(P, "***"),
                         P < 0.01 ~ paste(P, "**"),
                         P < 0.05 ~ paste(P, "*"),
                         TRUE ~ P %>% as.character())) 
}