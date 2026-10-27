#**********************************************************# 
#
# 'assoc_test_severe_gba1eff_controls_adjBMI_2026-10-26.R'
#* created: 2026-10-26
#* association tests between gba1 status and lipids
#* outfile: 
#* output file contains: 
# 
#**********************************************************# 
#
#0. initialize ---- 
library(dplyr)
library(tidyverse)
library(stringr)
library(data.table)
library(readr)
library(ggplot2)
library(ggforce)
library(patchwork)
## for table
library(table1)
require(survival)

describe = psych::describe
rename = dplyr::rename
select = dplyr::select

inormal <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

#1. read in data ----
wd="~/Library/CloudStorage/OneDrive-Personal/00_Work/GermanLipidomics/"
setwd(wd)

load("data/cleaned_lipid_data_list_and_covariate_data2026-10-08.Rdata",verbose = T)

severe.muts = c('IVS2+1g>a','H255Q','G202R','R463C','L444P','R359X','R398X','c.203delC')
nonsevere.muts = c("N370S", "N370S + S271G", "E326K", "A446P", "T369M", "T297S", 
                   "S271G", "N370S + T369M")
for(i in 4){
  chr_colnames_i = intersect(names(d_list[[i]]),chr_colnames)
  mnames = d_rm_high.missing_sd0[[i]] %>% select(-all_of(chr_colnames_i)) %>% names
  fit_res=c();for(j in 1:length(mnames)){
    mj = mnames[j]
    
    analdati = covariate_data %>% 
      left_join(
        d_rm_high.missing_sd0[[i]] %>% select(1,all_of(mj))
      ) %>% mutate(bmi.c = scale(bmi)/2) %>%mutate(bmi.c2 = bmi.c^2)
    #lipid variable - untransformed (y) and INT-transformed (m)
    analdati[['y']] = (analdati[[mj]])
    analdati[['m']] = inormal(analdati[[mj]])
    
    analdati = analdati %>% data.frame
    rownames(analdati) = analdati$shortname
    
    #*************** SELECT ONLY healthy controls ***************# 
    #** n=11 GBA1pos vs 136 GBA1neg
    sel_ind = !is.na(analdati$GBA_Mutation) & (analdati$GBA_Mutation %in% severe.muts);
    sel_ind = sel_ind|is.na(analdati$GBA_Mutation);print(sum(sel_ind))
    
    analdati = analdati %>% filter(sel_ind) %>%
      # filter(on_lipid_fin=="no") %>% 
      filter(PD_status==1) 
    table(is.na(analdati$GBA_Mutation),useNA='a')#11(positive)+72(negative)
    #*************************************************************#
    # all
    fit.lipid = lm(inormal(y)~(age.c + age.c2 + bmi.c + bmi.c2)*sex, 
                   data=analdati,
                   na.action = na.exclude)
    summary(fit.lipid)
    
    p = analdati %>% 
      ggplot(
        aes(x=Age_at_Visit,
            y=resid(fit.lipid),color=sex)) + 
      geom_point() + geom_smooth(method='gam') +
      ggpubr::stat_cor(show.legend = F)
    
    p2 = analdati %>% 
      ggplot(
        aes(x=bmi,
            y=resid(fit.lipid),color=sex)) + 
      geom_point() + geom_smooth(method='gam') +
      ggpubr::stat_cor(show.legend = F)
    p+p2
    fit1 = lm(inormal(resid(fit.lipid))~cohort2+(age.c + age.c2 + bmi.c + bmi.c2)*sex,
              data=analdati)
    summary(fit1)
    # 
    fit2 = lm(inormal(y)~cohort2+(age.c + age.c2 + bmi.c + bmi.c2)*sex, 
              data=analdati,na.action = na.exclude)
    summary(fit2)
    
    fit_res = rbind(fit_res,c(
      lipidID = mj, 
      summary(fit2)$coef['cohort2MIGAP',],
      n=nobs(fit2)))
  }
  colnames(fit_res)[2:5] = c('beta','se','z','pval')
  resi = fit_res %>% 
    data.table %>% 
    mutate(across(-lipidID,as.numeric)) %>%
    mutate(lipidID = as.character(lipidID)) %>%
    arrange(lipidID) %>%
    mutate(lipid.order = 1:length(mnames)) %>%
    mutate(fdr.p = p.adjust(pval,method = 'fdr')) %>%
    # mutate(fdr.p_interact = p.adjust(pval_interact,method = 'fdr')) %>%
    select(-contains("z"))
  
  resi_OR_file = paste0('results/assoc_test_gba1eff_severe_',names(d_rm_high.missing_sd0)[i],"_PDcases_adjBMI_",Sys.Date(),".txt")
  write_tsv(resi,resi_OR_file)
}
resi %>% arrange(pval) %>% head %>% print
resi %>% arrange(-beta) %>% filter(fdr.p<0.05) %>% print
# rm(resi,resi_OR_file)
