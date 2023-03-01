tidydata = function(expr, md, hgnccolnum, genestostudy){
  cat("Combining expression data with clinical data ....\n")
  if(exists("genestostudy")) {expr %<>% filter(expr[,hgnccolnum] %in% genestostudy)}
  expr <- expr %>% tibble::column_to_rownames(colnames(df)[hgnccolnum])
  if(all(sapply(expr[,1], is.character))){
    expr <- expr[, -1]
  }
  expr <- as.data.frame(t(expr))
  colnames(md)[1] <- 'id'
  if(!identical (rownames(expr), md$id)) {
    expr %<>% filter(rownames(expr) %in% md$id)
  }
  tryCatch(
    stopifnot(all(rownames(expr) %in% md$id)),
    error = function(e){
      message("Samples of the expression and metadata aren't equal")
    }
  )
  md <- cbind(md, expr)
  
  cat("Checking for NA values ....\n")
  md %<>% filter(md$os != 'NA')
  if(all(md$os != 'NA')) {cat("Data is clean. Congrats :p \n")}
  
  md %<>%
    mutate(status = as.numeric(
      ifelse(
        vitalstatus == 'DECEASED'| os == '1:DECEASED',
        '1',
        '0'
      )))
  
  colnames(md) <- gsub("-", "", colnames(md))
  return(md)
}

KMgraph = function(df, genestostudy){
  cat("Executing Kaplan Meier Analysis .... \n")
  df %<>% as.data.table()
  survdf <- Surv(time = df$os, event = df$status)
  custom_theme <- function() {
    theme_survminer() %+replace%
      theme(
        plot.title=element_text(hjust=0.5, size = 18),
        legend.box = "vertical"
      )}
  for (i in unique(1:length(genestostudy)))
  {
    x = genestostudy[i]
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
    ggsave(file = paste(x, ".png", sep = ""))
  }
  df %<>% as.data.frame()
  cat("Analysis completed and saved... \n")
}

KMunitable = function(kmdf, genestostudy){
  at("Executing Kaplan Meier analysis in table form .... \n")
  unicolumn = c('Gene', 'Hazard Ratio', '95% CI', 'P-value')
  unidf = data.frame(matrix(nrow=0, ncol = length(unicolumn)))
  colnames(unidf) = unicolumn
  
  for (i in unique(1:length(genestostudy)))
  {
    z = genestostudy[i]
    name = paste(z, 'Expression', sep = "")
    kmdf %>% dplyr::select(z) %>% unlist() -> c
    mdn = median(c)
    kmdf %<>% mutate("{name}" := as.factor(
      ifelse(
        c >= mdn,
        'High',
        'Low'
      )))}
  
  gtslist = grep("Expression", colnames(kmdf), value = TRUE)
  for (o in gtslist){
    kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
  }
  
  survdf <- Surv(time = kmdf$os, event = kmdf$status)
  
  for (p in unique(1:length(genestostudy)))
  {
    x = genestostudy[p]
    y = gtslist[p]
    f2 <- as.formula(paste('survdf ~', paste(y)))
    unicox = coxph(f2, data = kmdf)
    ucoef = summary(unicox)$coef[1,2]
    upvalue = summary(unicox)$coef[1,5]
    ci = confint(unicox, level = 0.95)
    ci <- exp(ci[1,])
    ci <- round(ci, digits = 2)
    uci <- paste(ci[1], '~', ci[2], sep = " ")
    unidf2 = data.frame(x, ucoef, uci, upvalue)
    colnames(unidf2) <- unicolumn
    unidf <- rbind(unidf, unidf2)}
  cat("Analysis Complete :> \n")
  return(unidf)
}

multicox = function(kmdf, genestostudy){
  cat("Executing Multivariate Cox Analysis ... \n")
  for(col_name in colnames(kmdf)) {
    if(tolower(col_name) %in% tolower(genestostudy)) {
      colnames(kmdf)[colnames(kmdf) == col_name] <- tolower(col_name)
        }
  }
  
  for (i in unique(1:length(genestostudy)))
  {
    z = genestostudy[i]
    l = tolower(genestostudy[i])
    name = paste(z)
    kmdf %>% dplyr::select(l) %>% unlist() -> c
    mdn = median(c)
    kmdf %<>% mutate("{name}" := as.factor(
      ifelse(
        c >= mdn,
        'High',
        'Low'
      )))}
  
  for (o in genestostudy){
    kmdf[[o]] <- relevel(factor(kmdf[[o]]), ref = 'Low')
  }
  
  survdf <- Surv(time = kmdf$os, event = kmdf$status)
  f1 <- as.formula(paste('survdf ~', paste(genestostudy, collapse = "+")))
  cox = coxph(f1, data = kmdf) %>% print()
  ggforest(cox, data = kmdf, refLabel = "Reference") %>% print()
  cat("Analysis Complete :> \n")
  return(cox)
}