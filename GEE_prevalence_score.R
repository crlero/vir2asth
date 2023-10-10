library(dplyr)
library(ggplot2)
library(gdata)
library(dplyr)
library(gtools)
library(rio)

theme_pufa <-   theme_bw() +
  theme(axis.title = element_text(margin = ggplot2::margin(20,20,20,20))) +
  theme(panel.grid.minor = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(legend.title =element_blank()) +
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  theme(legend.background = element_rect(colour = NA, fill = NA)) +
  theme(legend.key = element_rect(colour = NA, fill = NA)) +
  theme(plot.margin = unit(c(2,2,2,2),"points")) +
  theme(strip.background = element_blank(), strip.text=element_text(size=14)) +
  theme(legend.key.width = unit(0.5, "inches"))


# Last visit

Last_planned <- read.xls('/Volumes/UserFolders/Jakob/Projects/J45_AD_cox_cross/raw/Jakob20220511.xlsx',sheet = 'PlanLangt')
Last_planned <- Last_planned[order(rev(as.Date(Last_planned$DATO))), ] #Fjerner laveste dublet;
Last_planned <- Last_planned[!duplicated(Last_planned$ABCNO),]
colnames(Last_planned)[4:5] <- c("DATO_PLANNED","VISIT_PLANNED")
Last_planned <- Last_planned[c("ABCNO","BIRTHDATE","DATO_PLANNED","VISIT_PLANNED")] 

Last_adhoc <- read.xls('/Volumes/UserFolders/Jakob/Projects/J45_AD_cox_cross/raw/Jakob20220511.xlsx',sheet = 'Ad_hoc')
Last_adhoc <- Last_adhoc[order(rev(as.Date(Last_adhoc$DATO))), ] #Fjerner laveste dublet;
Last_adhoc <- Last_adhoc[!duplicated(Last_adhoc$ABCNO),]
colnames(Last_adhoc)[4:5] <- c("DATO_ADHOC","VISIT_ADHOC")
Last_adhoc <- Last_adhoc[c("ABCNO","DATO_ADHOC","VISIT_ADHOC")] 
Last_adhoc <- Last_adhoc[Last_adhoc$DATO_ADHOC!="00:00:00.00",]
Last_visit <- merge(Last_planned,Last_adhoc,by = 'ABCNO', all.x = T, all.y = T)

Last_visit$last_visit <- ifelse(((as.Date(Last_visit$DATO_PLANNED) > as.Date(Last_visit$DATO_ADHOC)) | is.na(Last_visit$DATO_ADHOC)),as.character(Last_visit$VISIT_PLANNED), as.character(Last_visit$VISIT_ADHOC))
Last_visit$last_date <- ifelse(((as.Date(Last_visit$DATO_PLANNED) > as.Date(Last_visit$DATO_ADHOC)) | is.na(Last_visit$DATO_ADHOC)),as.character(Last_visit$DATO_PLANNED), as.character(Last_visit$DATO_ADHOC))
Last_visit$visit_age <- as.Date(Last_visit$last_date) - as.Date(Last_visit$BIRTHDATE)
Last_visit_all <- Last_visit[c("ABCNO","BIRTHDATE","last_visit","last_date","visit_age")] 
Last_visit <- Last_visit_all[Last_visit_all$visit_age > 0,]
Last_visit <- Last_visit[!is.na(Last_visit$ABCNO),]

# Follow up

Last_visit$Follow_0yr<-as.numeric(1)
Last_visit$Follow_1yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= 365 -30| Last_visit$last_visit == "1 yr" | Last_visit$last_visit == "18 mth") & (Last_visit$visit_age >= 365 -30| Last_visit$last_visit == "1 yr" | Last_visit$last_visit == "18 mth"), "Follow_1yr"] <- as.numeric(1)
Last_visit$Follow_2yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (2*365)-30 | Last_visit$last_visit == "2 yrs" ) & (Last_visit$visit_age >= (2*365)-30 | Last_visit$last_visit == "2 yrs" ), "Follow_2yr"] <- as.numeric(1)
Last_visit$Follow_3yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (3*365)-30 | Last_visit$last_visit == "3 yrs") & (Last_visit$visit_age >= (3*365)-30 | Last_visit$last_visit == "3 yrs") , "Follow_3yr"] <- as.numeric(1)
Last_visit$Follow_4yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (4*365)-30 | Last_visit$last_visit == "4 yrs") & (Last_visit$visit_age >= (4*365)-30 | Last_visit$last_visit == "4 yrs"), "Follow_4yr"] <- as.numeric(1)
Last_visit$Follow_5yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (5*365)-30 | Last_visit$last_visit == "5 yrs") & (Last_visit$visit_age >= (5*365)-30 | Last_visit$last_visit == "5 yrs"), "Follow_5yr"] <- as.numeric(1)
Last_visit$Follow_6yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (6*365)-30 | Last_visit$last_visit == "6 yrs") & (Last_visit$visit_age >= (6*365)-30 | Last_visit$last_visit == "6 yrs"), "Follow_6yr"] <- as.numeric(1)
Last_visit$Follow_7yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (7*365)-30 | Last_visit$last_visit == "7 yrs") & (Last_visit$visit_age >= (7*365)-30 | Last_visit$last_visit == "7 yrs"), "Follow_7yr"] <- as.numeric(1)
Last_visit$Follow_8yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (8*365)-30 | Last_visit$last_visit == "8 yrs") & (Last_visit$visit_age >= (8*365)-30 | Last_visit$last_visit == "8 yrs"), "Follow_8yr"] <- as.numeric(1)
Last_visit$Follow_9yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (9*365)-30 | Last_visit$last_visit == "9 yrs") & (Last_visit$visit_age >= (9*365)-30 | Last_visit$last_visit == "9 yrs"), "Follow_9yr"] <- as.numeric(1)
Last_visit$Follow_10yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (10*365)-30 | Last_visit$last_visit == "10 yrs") & (Last_visit$visit_age >= (10*365)-30 | Last_visit$last_visit == "10 yrs"), "Follow_10yr"] <- as.numeric(1)
Last_visit$Follow_11yr<-as.numeric(0)
Last_visit[!is.na(Last_visit$visit_age >= (11*365)-30 | Last_visit$last_visit == "11 yrs") & (Last_visit$visit_age >= (11*365)-30 | Last_visit$last_visit == "11 yrs"), "Follow_11yr"] <- as.numeric(1)

Last_visit$Follow<-"1 year"
Last_visit[!is.na(Last_visit$visit_age >= 365 -30| Last_visit$last_visit == "1 yr" | Last_visit$last_visit == "18 mth") & (Last_visit$visit_age >= 365 -30| Last_visit$last_visit == "1 yr" | Last_visit$last_visit == "18 mth"), "Follow"] <- "2 years"
Last_visit[!is.na(Last_visit$visit_age >= (2*365)-30 | Last_visit$last_visit == "2 yrs" | Last_visit$last_visit == "30 mth") & (Last_visit$visit_age >= (2*365)-30 | Last_visit$last_visit == "2 yrs" | Last_visit$last_visit == "30 mth"), "Follow"] <- "3 years"
Last_visit[!is.na(Last_visit$visit_age >= (3*365)-30 | Last_visit$last_visit == "3 yrs") & (Last_visit$visit_age >= (3*365)-30 | Last_visit$last_visit == "3 yrs"), "Follow"] <- "4 years"
Last_visit[!is.na(Last_visit$visit_age >= (4*365)-30 | Last_visit$last_visit == "4 yrs") & (Last_visit$visit_age >= (4*365)-30 | Last_visit$last_visit == "4 yrs"), "Follow"] <- "5 years"
Last_visit[!is.na(Last_visit$visit_age >= (5*365)-30 | Last_visit$last_visit == "5 yrs") & (Last_visit$visit_age >= (5*365)-30 | Last_visit$last_visit == "5 yrs"), "Follow"] <- "6 years"
Last_visit[!is.na(Last_visit$visit_age >= (6*365)-30 | Last_visit$last_visit == "6 yrs") & (Last_visit$visit_age >= (6*365)-30 | Last_visit$last_visit == "6 yrs"), "Follow"] <- "7 years"
Last_visit[!is.na(Last_visit$visit_age >= (7*365)-30 | Last_visit$last_visit == "7 yrs") & (Last_visit$visit_age >= (7*365)-30 | Last_visit$last_visit == "7 yrs"), "Follow"] <- "8 years"
Last_visit[!is.na(Last_visit$visit_age >= (8*365)-30 | Last_visit$last_visit == "8 yrs") & (Last_visit$visit_age >= (8*365)-30 | Last_visit$last_visit == "8 yrs"), "Follow"] <- "9 years"
Last_visit[!is.na(Last_visit$visit_age >= (9*365)-30 | Last_visit$last_visit == "9 yrs") & (Last_visit$visit_age >= (9*365)-30 | Last_visit$last_visit == "9 yrs"), "Follow"] <- "10 years"
Last_visit[!is.na(Last_visit$visit_age >= (10*365)-30 | Last_visit$last_visit == "10 yrs") & (Last_visit$visit_age >= (10*365)-30 | Last_visit$last_visit == "10 yrs"), "Follow"] <- "11 years"
Last_visit[!is.na(Last_visit$visit_age >= (11*365)-30 | Last_visit$last_visit == "11 yrs") & (Last_visit$visit_age >= (11*365)-30 | Last_visit$last_visit == "11 yrs"), "Follow"] <- "12 years"

Last_visit$Follow <- factor(Last_visit$Follow, levels=c( "1 year","2 years", "3 years", "4 years","5 years","6 years","7 years","8 years","9 years","10 years","11 years","12 years"))

# Age at break-the-code day
Last_visit$breakthecode_age <- as.numeric(as.Date("2014/03/28")-as.Date(Last_visit$BIRTHDATE))
Last_visit$Follow_breakthecode <- 0
Last_visit[!is.na(Last_visit$visit_age >= Last_visit$breakthecode_age-30) & (Last_visit$visit_age >= Last_visit$breakthecode_age-30), "Follow_breakthecode"] <- as.numeric(1)

Follow <-  Last_visit[c("ABCNO","Follow_0yr","Follow_1yr","Follow_2yr","Follow_3yr","Follow_4yr","Follow_5yr","Follow_6yr","Follow_7yr","Follow_8yr","Follow_9yr","Follow_10yr","Follow_11yr","Follow", "Follow_breakthecode")]


Diagnoses <- rio::import('/Volumes/UserFolders/Jakob/Projects/J45_AD_cox_cross/raw/Jakob20220511.xlsx',which = 'J45JJJ')
unreliable_6yr <- rio::import('/Volumes/UserFolders/Jakob/Projects/J45_AD_cox_cross/raw/6yrs_missingdata_final.xlsx')
valid_6to10yr <- read.xls('/Volumes/UserFolders/Jakob/Projects/J45_AD_cox_cross/raw/Valid_10yr.xlsx', stringsAsFactors = F)
J45 <- Diagnoses[Diagnoses$DIAGNOSCODE == 'J45' | Diagnoses$DIAGNOSCODE == 'J45.990',] #Keep asthma + excercise asthma
#J45 <- Diagnoses[Diagnoses$DIAGNOSCODE == 'J45',] #Only keep J45 - asthma
J45$DIAGNOSCODE <- NULL
J45 <- J45[!is.na(as.Date(J45$J45S)), ]

J45_start <- J45[order(J45$ABCNO , as.Date(J45$J45S)), ] #Fjerner højeste dublet;
J45_start <- J45_start[!duplicated(J45_start$ABCNO),]
J45_start$J45E <-  NULL

J45_end <- J45[!is.na(as.Date(J45$J45E)), ]
J45_end <- J45_end[order(rev(as.Date(J45_end$J45E))), ] #Fjerner laveste dublet;
J45_end <- J45_end[!duplicated(J45_end$ABCNO),]
J45_end$J45S <-  NULL

#Nyt dataset med første start og sidste slut)
J45_cox <- merge(J45_start,J45_end,by = 'ABCNO',all =T)
J45_cox <- merge(Last_visit,J45_cox, by = 'ABCNO',all =T)
#J45_cox <- merge(J45_cox, Follow, by = 'ABCNO',all =T)


#Correct last age for Xyr visit before age Xyr.
J45_cox[J45_cox$visit_age < (5*365) & J45_cox$last_visit == "5 yrs" , "visit_age"] <- (5*365)
J45_cox[J45_cox$visit_age > (5*365) - 30 & J45_cox$visit_age < (5*365), "visit_age"] <- (5*365)
J45_cox[J45_cox$visit_age < (6*365) & J45_cox$last_visit == "6 yrs" , "visit_age"] <- (6*365)
J45_cox[J45_cox$visit_age > (6*365) - 30 & J45_cox$visit_age < (6*365), "visit_age"] <- (6*365)
J45_cox[J45_cox$visit_age < (8*365) & J45_cox$last_visit == "8 yrs" , "visit_age"] <- (8*365)
J45_cox[J45_cox$visit_age < (10*365) & J45_cox$last_visit == "10 yrs" , "visit_age"] <- (10*365)
J45_cox[J45_cox$visit_age > (10*365) - 30 & J45_cox$visit_age < (10*365), "visit_age"] <- (10*365)


# Adjust visit_age for children who recently came back.
J45_cox <- merge(J45_cox, unreliable_6yr[, c("abcno", "asthma_valid_days", "asthma_last_valid_visit")], by.x = "ABCNO", by.y = "abcno", all.x = T, all.y = F, sort = F)
J45_cox$last_visit_reentry <- J45_cox$last_visit
J45_cox$visit_age_reentry <- J45_cox$visit_age
J45_cox <- J45_cox %>%
  mutate(	visit_age = ifelse(is.na(asthma_valid_days) | asthma_valid_days > visit_age, visit_age,
                             asthma_valid_days),
          last_visit = ifelse(is.na(asthma_valid_days) | asthma_valid_days > visit_age, last_visit, asthma_last_valid_visit))

# Adjust visit_age for children who came back at 10yrs.
J45_cox <- merge(J45_cox, valid_6to10yr[valid_6to10yr$valid==0, c("ABCNO", "valid", "valid_days", "last_valid_visit")], by.x = "ABCNO", by.y = "ABCNO", all.x = T, all.y = F, sort = F)
J45_cox$last_visit_reentry <- J45_cox$last_visit
J45_cox$visit_age_reentry <- J45_cox$visit_age
J45_cox <- J45_cox %>%
  mutate(	visit_age = ifelse(is.na(valid_days) | valid_days > visit_age, visit_age,valid_days),
          last_visit = ifelse(is.na(valid_days) | valid_days > visit_age, last_visit, last_valid_visit))


#Remove diagnosis end_dates if after last_visit
J45_cox[!is.na(J45_cox$J45E) & as.Date(J45_cox$J45E) > as.Date(J45_cox$BIRTHDATE)+J45_cox$visit_age,"J45E"]<-NA


# Set up predictor

mygee = function(X, var, cutoff) {
  X$predictorX <- quantcut(X[[var]], q=2)
  prev <- merge(X[!is.na(X$predictorX),c("abcno","predictorX")], J45_cox, by.x="abcno", by.y="ABCNO", all.x=F, all.y=F)
  
  
  prev$predictor <- as.numeric(factor(prev$predictorX))-1
  prev$visit_age <- as.numeric(prev$visit_age)
  prev$J45S <- as.Date(prev$J45S)
  prev$J45E <- as.Date(prev$J45E)
  prev$BIRTHDATE <- as.Date(prev$BIRTHDATE)
  prev$startage <- as.numeric(with(prev, J45S - BIRTHDATE))
  prev$endage <- as.numeric(with(prev, J45E - BIRTHDATE))
  prev[is.na(prev$endage) & !is.na(prev$startage),"endage"] <- prev[is.na(prev$endage) & !is.na(prev$startage),"visit_age"]  # Close diagnosis when completing followup
  
  
  # Yearly point prevalences
  startyear <- 0
  endyear <- cutoff
  intlen <- 365/12*cutoff  # Interval length
  stoptime <- cutoff*365
  pointprev <- expand.grid(years = startyear:endyear, predictor = 0:1)
  pointprev$prev <- NA
  prev[prev$Follow_10yr == TRUE & prev$visit_age < stoptime-intlen,]
  for(i in startyear:endyear){
    ncases <- sum(prev[prev$predictor == 1,]$startage <= i*365 & prev[prev$predictor == 1,]$endage >= i*365, na.rm=T)
    ntot <- sum(prev[prev$predictor == 1,]$visit_age >= i*365, na.rm=T)
    pointprev[with(pointprev,years == i & predictor == 1),"ncases"] <- ncases
    pointprev[with(pointprev,years == i & predictor == 1),"ntot"] <- ntot
    pointprev[with(pointprev,years == i & predictor == 1),"prev"] <- ncases/ntot
    
    ncases <- sum(prev[prev$predictor == 0,]$startage <= i*365 & prev[prev$predictor == 0,]$endage >= i*365, na.rm=T)
    ntot <- sum(prev[prev$predictor == 0,]$visit_age >= i*365, na.rm=T)
    pointprev[with(pointprev,years == i & predictor == 0),"ncases"] <- ncases
    pointprev[with(pointprev,years == i & predictor == 0),"ntot"] <- ntot
    pointprev[with(pointprev,years == i & predictor == 0),"prev"] <- ncases/ntot
  }
  
  # Table
  library(geepack)

  intlen <- 365 # Interval length
  stoptime <- cutoff*365
  cuts <- seq(365/2,stoptime,intlen)
  cutlist <- list()
  for(i in 2:length(cuts)){
    listobj <- prev[,c("abcno", "predictor", "visit_age", "last_visit", "Follow_10yr", "startage", "endage")]
    listobj$period <- i-1
    listobj$time <- paste0(cuts[i-1]," - ",cuts[i])
    listobj$curr_asthma <- with(listobj, ifelse(!is.na(startage) & (startage < cuts[i] & endage > cuts[i-1]),1,0))
    
    listobj[with(listobj,!visit_age > cuts[i-1] & !is.na(curr_asthma)), "curr_asthma"] <- NA
    cutlist[[i-1]] <- listobj
  }
  geedat <- do.call("rbind",cutlist)
  modeldat <- geedat[with(geedat, order(abcno,period)),]
  fit <- geeglm(curr_asthma ~ factor(period) + factor(predictor), id=abcno, data=modeldat, family = binomial(link="logit"), corstr="independence")
  fits <- summary(fit)
  hr <- exp(fits$coefficients["factor(predictor)1",1])
  hr_low <- exp(fits$coefficients["factor(predictor)1",1] - 1.96*fits$coefficients["factor(predictor)1",2])
  hr_high <- exp(fits$coefficients["factor(predictor)1",1] + 1.96*fits$coefficients["factor(predictor)1",2])
  pval <- fits$coefficients["factor(predictor)1",4]
  tfig3 <- data.frame(text = "Interval prevalence GEE (HR is OR here)",fig = "fig3", var = "PUFA", HR = hr, HRlow = hr_low, HR_high = hr_high, p = pval)
  # Plot
  HRtext <- paste0("GEE OR ", sprintf("%.2f",round(hr,2)), " [", sprintf("%.2f",round(hr_low,2)), "-",sprintf("%.2f",round(hr_high,2)),"], P ",ifelse(pval<0.001,"< 0.001",paste0("= ",format.pval(pval,1,0.001,nsmall=2))))
  
  
  gee_fig <- pointprev %>% mutate(predictor = ifelse(predictor == 0, "Low score", "High score")) %>%
    ggplot(.,aes(x = years, y = prev, linetype = factor(predictor),
                 color=factor(predictor))) + 
    geom_point() +
    geom_line(size=0.7) +
    scale_color_manual(values = c("Low score" = "#fdbb63", "High score" = "#f6511d")) +
    scale_y_continuous("Prevalence of persistent\nwheeze / asthma",labels=scales::percent, limits = c(0,0.40))+
    scale_x_continuous(breaks=seq(0,10,1)) +
    scale_linetype_discrete(name=NULL) +
    xlab("Age (years)") +
    guides(linetype="none") +
    annotate("text", x = 2.5 ,y=0.4, label=HRtext, size = 4) +
    theme_pufa
  
  return(gee_fig)
}