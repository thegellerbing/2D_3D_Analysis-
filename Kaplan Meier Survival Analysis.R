library(pacman)
p_load(tidyverse, data.table, survival, survminer, magrittr, readxl)

#Make sure the columns containing os and deceased status are named 'os' and 'vitalstatus' respectively
source("C:/Users/PC/OneDrive/Desktop/PhD Folder/Projects/R/Codes/Functions/KM_func.R") 

exprdata #Normalised Expression Data
md #Survival Data
hgnccolnum <- 2 #Column number with hgnc symbol

df <- tidydata(exprdata, md, hgnccolnum, genestostudy)

#If graphs are needed ####
KMGraph(df, genestostudy) #Graphs will be saved in the current directory

#Univariate KM Analysis ####
unidf <- KMunitable(df, genestostudy)

#Multivariate Cox Analysis ####
cox <- multicox(df, genestostudy)
