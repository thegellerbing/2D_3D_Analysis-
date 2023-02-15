#Library and Directory ####
library(pacman)
p_load(tidyverse, data.table, survival, survminer, magrittr, readxl, survivalROC, pROC, OptimalCutpoints)


setwd("C:/Users/PC/OneDrive/Documents/R/TCGA ABC Transporters/Microarray")

#Remove Rows with NA

#Splitting Survival Status ####
df %<>%
  mutate(status = as.numeric(
    ifelse(
      CDEvitalstatus == 'DECEASED',
      '1',
      '0'
    ))) %>%
  print()

#Median ####
mdn = median(df[,FN1]) %>% print()

#Splitting Into Groups Based on Z-Score/Median ####
df %<>%
  mutate(expression = as.factor(
    ifelse(
      FN1 >= mdn,
      'High',
      'Low'
    ))) %>%
  print()

#Kaplan Meier Survival Curve #### ALWAYS run the custom_theme first before plotting the survival plot if you're using my code
survdf <- Surv(time = df$CDEsurvivaltime, event = df$status) 
fit1 = surv_fit(survdf ~ expression, data = df) %>% print()

custom_theme <- function() {
  theme_survminer() %+replace%
    theme(
      plot.title=element_text(hjust=0.5, size = 18),
      legend.box = "vertical"
    )}

ggsurvplot(fit1, data=df, 
           pval = TRUE, pval.coord = c(1250,0.3), pval.size = 5.5,
           break.time.by = 500, 
           xlab = "Time (Days)", ylab = "Survival Probability",
           palette = c( "firebrick1","steelblue2"),  
           title = "CD44 Survival Curve", ggtheme = custom_theme(),
           legend.labs = c("High Expression (n = 77) : 323 days",
                           "Low Expression (n = 77) : 459 days"), 
           legend.title = "Median Survival", legend = c(0.58,0.75),
           font.main = c(18,"bold"), font.legend = c(16),
           font.x = c(16, "bold"), font.y = c(16, "bold"), font.tickslab = c(15)
           )  
ggsave("CD44.png")

for (i in unique(1:length(list)))
{
  x = list[i]
  df %>% dplyr::select(x) %>% unlist() -> y
  mdn = median(y) 
  df %<>%
    mutate(expression = as.factor(
      ifelse(
        y >= mdn,
        'High',
        'Low'
      ))) 
  fit1 = surv_fit(survdf ~ expression, data = df)
  shigh = round(unname(summary(fit1)$table[,'median'][1]),0)
  slow = round(unname(summary(fit1)$table[,'median'][2]),0)
  nhigh = unname(summary(fit1)$table[,'records'][1])
  nlow = unname(summary(fit1)$table[,'records'][2])
  graph = ggsurvplot(fit1, data=df, 
             pval = TRUE, pval.coord = c(1250,0.3), pval.size = 5.5,
             break.time.by = 500, 
             xlab = "Time (Days)", ylab = "Survival Probability",
             palette = c( "firebrick1","steelblue2"),  
             title = (paste(x, " Survival Curve", sep = "")), ggtheme = custom_theme(),
             legend.labs = c(paste("High Expression (n= ", nhigh, ")", " : ", shigh, " days", sep = ""),
                             paste("Low Expression (n= ", nlow, ")", " : ", slow, " days", sep = "")), 
             legend.title = "Median Survival", legend = c(0.58,0.75),
             font.main = c(18,"bold"), font.legend = c(16),
             font.x = c(16, "bold"), font.y = c(16, "bold"), font.tickslab = c(15))
  print(graph)
  ggsave(file = paste(x, ".png", sep = ""))
}

#Backup with Risk Table
plot = ggsurvplot(fit1, data=df, 
           pval = TRUE, pval.coord = c(2500,0.3), pval.size = 4.5,
           break.time.by = 500, 
           risk.table = TRUE, risk.table.y.text = FALSE, 
           xlab = "Time (Days)", ylab = "Survival Probability",
           palette = c( "firebrick1","steelblue2"),
           surv.median.line = "hv",
           title = "CD44 Survival Curve", ggtheme = custom_theme(),
           legend.labs = c("High Expression","Low Expression"), 
           legend.title = "Median Survival", legend = c(0.65,0.75),
           font.main = c(13,"bold"), font.legend = c(12),
           font.x = c(12, "bold"), font.y = c(12, "bold")
            )  
print(plot)
png("CD44 RT.png")
print(plot, newpage = FALSE)
dev.off()


#Cox Proportional Hazards Model ####
cox = coxph(survdf ~ Expression + `IDH1 Status` + `MGMT Status` + `Therapy Type`, data = df)
ggforest(cox, data = df)

for (i in unique(1:length(list2)))
{
  x = list2[i]
  name = paste(x, 'Expression', sep = "")
  df %>% dplyr::select(x) %>% unlist() -> y
  mdn = median(y)
  df %<>% mutate("{name}" := as.factor(
    ifelse(
      y >= mdn,
      'High',
      'Low'
    )))}

list3 = colnames(df)[37:40]
for (i in list3){
  df[[i]] <- relevel(factor(df[[i]]), ref = 'Low')
}
f1 <- as.formula(paste('survdf ~', paste(list3, collapse = "+")))
cox = coxph(f1, data = df) %>% print()

#Large Scale Univariate Analysis in a Dataframe #####
unicolumn = c('Gene', 'Hazard Ratio', '95% CI', 'P-value')
unidf = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
colnames(unidf) = unicolumn

for (i in unique(1:length(list2)))
{
  x = list2[i]
  y = list3[i]
  f2 <- as.formula(paste('survdf ~', paste(y)))
  unicox = coxph(f2, data = df)
  ucoef = summary(unicox)$coef[1,2]
  upvalue = summary(unicox)$coef[1,5]
  ci = confint(unicox, level = 0.95)
  ci <- exp(ci[1,])
  ci <- round(ci, digits = 2)
  uci <- paste(ci[1], '~', ci[2], sep = " ")
  unidf2 = data.frame(x, ucoef, uci, upvalue)
  colnames(unidf2) <- unicolumn
  unidf <- rbind(unidf, unidf2)
}
rm(unicolumn, ucoef,upvalue, unidf2, ci, uci)

#Survival ROC ####
roc = survivalROC(Stime = df$CDEsurvivaltime, status = df$status, marker = df$risk_factor,
                  predict.time = 365, method = 'KM')
plot(roc$FP, roc$TP, type = "l", xlim = c(0,1), ylim = c(0,1), 
     xlab = paste( "Specificity \n AUC = ", round(roc$AUC, 3)),
     ylab = "Sensitivity", main = "24 Months",
     log = "",
     abline(0,1))

roc = timeROC(df$CDEsurvivaltime, df$status, df$risk_score3, other_markers = NULL, )

#Factor Count ####
df %>%
  pull(CD44_expression) %>%
  fct_count() %>% as.data.table() 



for (i in unique(1:length(list))){
  x = list[i]
  df %>% dplyr::select(x) %>% unlist() -> y
  mdn = median(y)
  if(mdn == 0){df %<>% dplyr::select(-x)}
}
