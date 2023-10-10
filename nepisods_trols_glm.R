load(file="Diary/trols_diary_data.rdata")
library(broom)
library(gtools)

myepi = function(X, var, var2, categorical=FALSE, interaction=FALSE) {
  if (categorical) {
    X_trols <- left_join(diary_trols_all, X %>% mutate(abcno = as.numeric(abcno)))
    X_trols[[var]] <- quantcut(X_trols[[var]], q=2)
    X_trols[[var]]  <- as.factor(as.numeric(X_trols[[var]])-1)
    X_trols$var <- X_trols[[var]] 
    X_trols$var2 <- X_trols[[var2]]
    X_trols$asthma_mother
  } else {
    X_trols <- left_join(diary_trols_all, X %>% mutate(abcno = as.numeric(abcno)))
    X_trols$var <- X_trols[[var]] 
    X_trols$var2 <- X_trols[[var2]]
  }
  test <- X_trols %>% dplyr::select(c(abcno, var, var2, asthma_mother,
                                      trols_episodes,
                                      trols_episodes_year_1,
                                      trols_episodes_year_2,
                                      trols_episodes_year_3,
                                      total_days,
                                      days_year_1,
                                      days_year_2,
                                      days_year_3)) %>%
    gather(variable, value, -abcno, -var, -var2, -asthma_mother, -total_days, -days_year_1, -days_year_2, -days_year_3)
  test[grep("_1", test$variable),c("year","days")] <- c("1",test[grep("_1", test$variable),"days_year_1"])
  test[grep("_2", test$variable),c("year","days")] <- c("2",test[grep("_2", test$variable),"days_year_2"])
  test[grep("_3", test$variable),c("year","days")] <- c("3",test[grep("_3", test$variable),"days_year_3"])
  test[is.na(test$year),c("year","days")] <- c("all",test[is.na(test$year),"total_days"])
  
  if (categorical & !interaction) {
    ged1 <- test %>% group_by(variable) %>%
      do(glm(as.numeric(value) ~ var,
             family = quasipoisson,
             na.action = na.exclude,
             data = .) %>%
           tidy(exponentiate = T, conf.int = T) %>%
           filter(term != "(Intercept)"))  
  } else if (!categorical & !interaction){
    ged1 <- test %>% group_by(variable) %>%
      do(glm(as.numeric(value) ~ scale(var) + offset(log(days+1)),
             family = quasipoisson,
             na.action = na.exclude,
             data = .) %>%
           tidy(exponentiate = T, conf.int = T) %>%
           filter(term == "scale(var)"))
  } else if (!categorical & interaction){
    ged1 <- test %>% group_by(variable) %>%
      do(glm(as.numeric(value) ~ scale(var)*scale(var2) + offset(log(days+1)),
             family = quasipoisson,
             na.action = na.exclude,
             data = .) %>%
           tidy(exponentiate = T, conf.int = T))
    print(ged1 |> 
            filter(term == "scale(var):scale(var2)"))
    ged1 <- ged1 |> 
           filter(term == "scale(var)")
  }
  
  ged1$name <- c("Total", "Year 1", "Year 2", "Year 3")
  ged1$type <- rep("Episodes", 4)
  ged1$color <- c("All", rep("Yearly", 3))
  
  qpois_fig <- ggplot(ged1, aes(x=name, y=estimate, ymin=conf.low, ymax=conf.high,
                                color=color)) +
    geom_hline(aes(yintercept=1), lty=2, color="black") +
    geom_errorbar(width=0.2) +
    geom_point(shape=18, size=3) +
    geom_text(aes(y=max(conf.high)+0.025,
                  label=ifelse(p.value > 0.001,
                               paste0("P = ", roundex(p.value, 3)),
                               paste0("P < 0.001"))),
              color="black",
              size=3,
              position=position_dodge(width=0.80)) +
    scale_color_tableau() +
    ylab("Incidence Risk Ratio (95% CI)") +
    xlab(NULL) +
    theme_bw(base_size=12) +
    scale_y_log10() +
    theme(panel.grid.major = element_line(colour = "grey", size=0.1),
          legend.position = "none",
          axis.title.y = element_text(size=12))
  return(qpois_fig)
}
