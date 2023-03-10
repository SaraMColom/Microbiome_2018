# Additional functions


# Aesthetics
Tx<-theme(axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 20)) +
  theme(axis.text.x = element_text(vjust = 1, hjust=1, angle=0, size = 20),
        axis.title.x = element_text(angle=0, size = 12),
        plot.title=element_text(size = 25,hjust=0))

# Aesthetics
Tx2<-theme(axis.text.y = element_text(size = 12),
           axis.title.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(vjust = 1, hjust=1, size = 25),
        axis.title.x = element_text(size = 25),
        plot.title=element_text(size = 25,hjust=0))

GoldGrey <- c("#F1CE63", "#79706E")
GreenBlue <- c("#59A14F", "#4E79A7")

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