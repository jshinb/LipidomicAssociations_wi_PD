#-----------------------------------------@
#* File name: [generate]_forestplot_252lipids_GBA1eff_severe_nonsevere_cases_or_controls_unadjBMI_2026-10-22.R
#* Created: 2026-10-22
#* Output files:
#*   forestplot_252lipids_XX_XXXX-XX-XX.pdf
#*   forestplot_252lipids_XX_XXXX-XX-XX.pdf
#-----------------------------------------@

# ============================================================@
# 1. Working directory ----
# ============================================================@
wd="~/Library/CloudStorage/OneDrive-Personal/00_Work/GermanLipidomics/"
setwd(wd)

# ============================================================@
# 2. Libraries and helper functions ----
# ============================================================@
library(dplyr)
library(stringr)
library(data.table)
library(ggplot2)
library(patchwork)
library(readr)
library(tidyverse)
library(ggforce)
library(ggstats)
library(ggforestplot)
library(ggpattern)#geom_stripped_rows
select = dplyr::select
describe = psych::describe

clean_lipid_id <- function(x) {
  x |>
    as.character() |>
    str_remove("\\*") |>
    str_remove_all(";\\d+") |>
    str_replace("_O-", "-O") |>
    str_replace("_", "") |>
    str_replace_all("_", "/")
}
generate_forest_col = function(data,is.OR=FALSE,my.pointsize=1.25,
                               my.width=0.5,my.stroke=0.5,dodge.width=0){
  
  column_forestplot <- ggplot(
    data,
    aes(x = .estimate, y = lipidID)
  ) +
    # 1. Background, Guidelines, and Base Theme from ggforestplot
    theme_forest() +
    geom_stripped_rows(
      odd = "grey85", 
      even = "white", 
      alpha = 0.2
    ) +
    geom_vline(
      xintercept = if(is.OR) 1 else 0, 
      linetype = "solid", 
      linewidth = 0.4, 
      colour = "black"
    ) +
    
    # 2. FIXED Connecting Lines Layer (No loops, smooth class transitions)
    geom_line(
      aes(
        x = beta,
        y = lipidID,
        colour = cohort,
        group = cohort         # Groups strictly by cohort to allow TAG -> DAG connection
      ),
      orientation = "y",       # Forces R to sort vertically along the Y-axis and stop at the ends
      position = position_dodge(width = dodge.width),
      linewidth = 0.5,
      alpha = 0.35
    ) +
    
    # 3. Shaded Polygon Layer (Thickened transparent range lines)
    geom_linerange(
      aes(
        xmin = .xmin, 
        xmax = .xmax, 
        colour = cohort,
        alpha = p.signif,
        y = lipidID      
      ),
      position = position_dodge(width = dodge.width), 
      # alpha = 0.5,              
      linewidth = 0.5           
    ) + scale_alpha_manual(values = c('p<0.05'=0.9,'p>=0.05'=0.2)) + 
    
    # 4. Central Points Layer
    geom_point(
      aes(
        x = beta,           
        y = lipidID,        
        colour = cohort,
        shape = fdr.signif,
        fill = fdr.signif
      ),
      position = position_dodge(width = dodge.width),
      size = my.pointsize,     
      stroke = my.stroke) +
    
    # 5. Professional Color Customization
    scale_color_manual(values = plot.cols) +  
    
    # 6. Fine-Tuning Text Layout Styles
    theme(
      axis.text.y = element_text(size = 5), 
      legend.position = "none"
    )
  return(column_forestplot)
}
generate_forest_col_by_adj = function(data,is.OR=FALSE,my.pointsize=1.25,
                                      my.width=0.5,my.stroke=0.5,dodge.width=0.5){
  
  column_forestplot <- ggplot(
    data,
    aes(x = .estimate, y = lipidID)
  ) +
    # 1. Background, Guidelines, and Base Theme from ggforestplot
    theme_forest() +
    geom_stripped_rows(
      odd = "grey85", 
      even = "white", 
      alpha = 0.2
    ) +
    geom_vline(
      xintercept = if(is.OR) 1 else 0, 
      linetype = "solid", 
      linewidth = 0.4, 
      colour = "black"
    ) +
    
    # 2. FIXED Connecting Lines Layer (No loops, smooth class transitions)
    geom_line(
      aes(
        x = beta,
        y = lipidID,
        colour = adj,
        group = adj         # Groups strictly by cohort to allow TAG -> DAG connection
      ),
      orientation = "y",       # Forces R to sort vertically along the Y-axis and stop at the ends
      position = position_dodge(width = dodge.width),
      linewidth = 0.5,
      alpha = 0.35
    ) +
    
    # 3. Shaded Polygon Layer (Thickened transparent range lines)
    geom_linerange(
      aes(
        xmin = .xmin, 
        xmax = .xmax, 
        colour = adj,
        alpha = p.signif,
        y = lipidID      
      ),
      position = position_dodge(width = dodge.width), 
      # alpha = 0.5,              
      linewidth = 0.5           
    ) + scale_alpha_manual(values = c('p<0.05'=0.9,'p>=0.05'=0.2)) + 
    
    # 4. Central Points Layer
    geom_point(
      aes(
        x = beta,           
        y = lipidID,        
        colour = adj,
        shape = fdr.signif,
        fill = fdr.signif
      ),
      position = position_dodge(width = dodge.width),
      size = my.pointsize,     
      stroke = my.stroke) +
    
    # 5. Professional Color Customization
    scale_color_manual(values = plot.cols) +  
    
    # 6. Fine-Tuning Text Layout Styles
    theme(
      axis.text.y = element_text(size = 5), 
      legend.position = "none"
    )
  return(column_forestplot)
}

# 3. lipid class definition ----
lipid_class = fread('Class252lipids.txt')
class.levels = c(
  "TAG",#1
  "DAG",#2
  "CE",#3
  "Chol",#4
  "PC",#5
  "PC-O",#6
  "LPC",#10
  "PE",#7
  "PE-O",#8
  "LPE",#11
  "LPE-O",#12
  "PI",#9
  "Cer",#13
  "HexCer",#14
  "SM"#15
)
lipid_class$lipidClass = str_replace(lipid_class$lipidClass,"_O[-]","-O")

lipid_class$lipidClass %>% table(useNA = 'a')
lipid_class = lipid_class %>% mutate(lipidClass = factor(lipidClass,levels=class.levels))
(table(lipid_class$lipidClass))
head(lipid_class)
lipid_class = lipid_class %>% mutate(lipidID2 = clean_lipid_id(lipidID))
lipid_class %>% head
lipid_class %>% tail
lipid_class %>% filter(lipidClass=="PC") %>% select(lipidID2)

# 4. load BMI-adjusted and -unadjusted results ----
# n=73
d_PD_GBA1_severe = fread("results/assoc_test_gba1eff_severe_Species-Class_PDcases_unadjBMI_2026-10-26.txt")
range(d_PD_GBA1_severe %>% pull(n))#54-77
# n=111
d_PD_GBA1_nonsevere = fread("results/assoc_test_gba1eff_nonsevere_Species-Class_PDcases_unadjBMI_2026-10-26.txt")
range(d_PD_GBA1_nonsevere %>% pull(n))#73-101
# n
d_HC_GBA1_severe = fread("results/assoc_test_gba1eff_severe_Species-Class_Controls_unadjBMI_2026-10-26.txt")
range(d_HC_GBA1_severe %>% pull(n))#50-71 (60-neg vs. 18 pos)
#
d_HC_GBA1_nonsevere = fread("results/assoc_test_gba1eff_nonsevere_Species-Class_Controls_unadjBMI_2026-10-26.txt")
range(d_HC_GBA1_nonsevere %>% pull(n))#58-78 (60-neg vs. 18 pos)

mycohorts = c(
  'GBA1eff: PD cases severe vs. no mutations',
  'GBA1eff: PD cases non-severe vs. no mutations',
  'GBA1eff: Controls severe vs. no mutations',
  'GBA1eff: Controls non-severe vs. no mutations'
)

df1 = ## PD effects on lipids #NOTE fdr.p = fdr.cohort
  d_PD_GBA1_severe %>% select(lipidID,beta,se,n,pval,fdr.p) %>% 
  mutate(cohort="GBA1eff: PD cases severe vs. no mutations") %>%
  mutate(fdr.cohort = p.adjust(pval,method="fdr")) %>%
  
  bind_rows(
    d_PD_GBA1_nonsevere %>% select(lipidID,beta,se,n,pval,fdr.p) %>% 
      mutate(cohort="GBA1eff: PD cases non-severe vs. no mutations") %>%
      mutate(fdr.cohort = p.adjust(pval,method="fdr"))
  ) %>%
  
  bind_rows(
    d_HC_GBA1_severe %>% select(lipidID,beta,se,n,pval,fdr.p) %>% 
      mutate(cohort='GBA1eff: Controls severe vs. no mutations') %>%
      mutate(fdr.cohort = p.adjust(pval,method="fdr"))
  ) %>%
  
  bind_rows(
    d_HC_GBA1_nonsevere %>% select(lipidID,beta,se,n,pval,fdr.p) %>% 
      mutate(cohort='GBA1eff: Controls non-severe vs. no mutations') %>%
      mutate(fdr.cohort = p.adjust(pval,method="fdr"))
  ) %>%
  
  right_join(lipid_class, join_by(lipidID)) %>% 
  
  relocate(lipidClass) %>% 
  
  mutate(cohort = factor(cohort,levels=rev(mycohorts))) %>%
  
  mutate(lipidClass = factor(lipidClass,levels=class.levels)) %>%
  
  select(-fdr.cohort) %>% #NOTE fdr.p = fdr.cohort
  
  mutate(adj = 'unadjBMI')

df1 %>% head
df1 %>% tail
df1 %>% pull(cohort) %>% table
df1 %>% pull(cohort) %>% head

# below not wroking any longer
# graphics.off();df1 %>% select(fdr.p,fdr.cohort) %>% plot(pch=20);abline(0,1,col="red")

## 4b. un-adjusted for BMI ----
# load(file="results/df1_unadjBMI.Rdata",verbose = T); df1_unadjBMI = 
#   df1_unadjBMI %>%
#   mutate(adj="unadjBMI",
#          lipidID2 = clean_lipid_id(lipidID)) 
# 
# ## 4c. combine across adjustments ----
# setdiff(names(df1),names(df1_unadjBMI))
# df_wi_wo_adjBMI = df1 %>% 
#   bind_rows(
#     df1_unadjBMI %>% select(all_of(names(df1)))
#     )
# 
# df_wi_wo_adjBMI %>% select(cohort) %>% head
# df_wi_wo_adjBMI %>% select(cohort,adj) %>% table (useNA = 'a')
# rm(df1,df1_unadjBMI)

# 5. Re-organized forest plots: fdr-p calc within groups (2026-10-13) ----
plot.mycohorts = c('GBA1eff: PD cases severe vs. no mutations',
                   'GBA1eff: Controls severe vs. no mutations')
myadjs = c("unadjBMI")
plot.cols = c("#00BA38","#619CFF");names(plot.cols) = plot.mycohorts
lipid_class = lipid_class %>% arrange(lipidClass)
lipidID.levels = lipid_class %>% pull(lipidID2)

## format plot data 
is.OR=FALSE;df_intervals <- df1 %>% 
  filter(adj=='unadjBMI') %>%
  filter(cohort %in% plot.mycohorts) %>% 
  #mutate(lipidID = clean_lipid_id(lipidID)) %>%
  mutate(
    lipidID = factor(lipidID2,levels=rev(lipidID.levels)),
    lipidClass = factor(lipidClass, levels=(class.levels)),
    cohort = factor(cohort,levels = rev(mycohorts)),
    adj = factor(adj, levels=rev(myadjs))
  ) %>%
  dplyr::mutate(
    # Use standard 1.96 * se for 95% CI
    .xmin = if(is.OR) exp(beta - 1.96 * se) else (beta - 1.96 * se),
    .xmax = if(is.OR) exp(beta + 1.96 * se) else (beta + 1.96 * se),
    .estimate = if(is.OR) exp(beta) else beta,
    # Match the ggforestplot hollow circle logic for non-significant points
    .filled = pval < 0.05
  ) %>% 
  mutate(cohort=factor(cohort,levels=rev(mycohorts)),
         fdr.signif = ifelse(fdr.p <0.05,"p.fdr<0.05","p.fdr>=0.05"),
         p.signif = ifelse(pval<0.05,"p<0.05","p>=0.05")) %>%
  arrange(lipidClass,lipidID,adj,cohort)
df_intervals %>% head

## format plot attributes
my.pointsize=2
my.width = 0.5
my.alpha = 0.6
my.stroke = 0.6

## 5a. column 1 ----
plot.lipid.classes=c("TAG","DAG"); 
p1 = generate_forest_col(data = df_intervals %>% 
                           filter(cohort %in% plot.mycohorts) %>%
                           filter(lipidClass %in% plot.lipid.classes) %>% 
                           arrange(lipidClass, lipidID) %>% 
                           slice(rev(row_number())));p1

## 5b. column 2: (50 lipids from CE and PC) ----
plot.lipid.classes=c("CE","PC"); 
plot.lipidIDs = df_intervals %>% filter(lipidClass %in%plot.lipid.classes)
plot.lipidIDs  = plot.lipidIDs  %>% arrange(lipidID) %>% pull(lipidID) %>% unique 
plot.lipidIDs_sel = rev(plot.lipidIDs)[1:50]
plot.lipidIDs_unsel = rev(plot.lipidIDs)[-c(1:50)]

p2 <- generate_forest_col(data = df_intervals %>% 
                            filter(cohort %in% plot.mycohorts) %>%
                            filter(lipidID %in% plot.lipidIDs_sel)); p2 

## 5c. column 3: (50 lipids from PC and PC-O) ----
plot.lipid.classes=c("PC-O"); 
plot.lipidIDs = df_intervals %>% filter(lipidClass %in%plot.lipid.classes)
plot.lipidIDs  = c(plot.lipidIDs  %>% pull(lipidID) %>% unique)
plot.lipidIDs = c(plot.lipidIDs_unsel, plot.lipidIDs) %>% sort
plot.lipidIDs_sel = rev(plot.lipidIDs)[1:50]
plot.lipidIDs_unsel = rev(plot.lipidIDs)[-c(1:50)]

p3 <- generate_forest_col(
  data = df_intervals %>% 
    filter(cohort %in% plot.mycohorts) %>%
    filter(lipidID %in% plot.lipidIDs_sel));p3
## 5d. column 4: (50 lipids from PC-O, LPC, PE) ----
plot.lipid.classes=c("LPC","PE"); 
plot.lipidIDs = df_intervals %>% filter(lipidClass %in%plot.lipid.classes)
plot.lipidIDs  = c(plot.lipidIDs  %>% pull(lipidID) %>% unique)
plot.lipidIDs = c(plot.lipidIDs_unsel, plot.lipidIDs) %>% sort
plot.lipidIDs_sel = rev(plot.lipidIDs)[1:50]
plot.lipidIDs_unsel = rev(plot.lipidIDs)[-c(1:50)]#0

p4 <- generate_forest_col(
  data = df_intervals %>% 
    filter(cohort %in% plot.mycohorts) %>%
    filter(lipidID %in% plot.lipidIDs_sel));p4

## 5e. column 5: ( lipids from PE-O, LPE, LPE-O, PI, Cer,SM) ----
plot.lipid.classes=c("PE-O","LPE","PI","Cer","SM"); 
plot.lipidIDs = df_intervals %>% filter(lipidClass %in%plot.lipid.classes)
plot.lipidIDs  = c(plot.lipidIDs  %>% pull(lipidID) %>% unique)
plot.lipidIDs = c(plot.lipidIDs_unsel, plot.lipidIDs) %>% sort
plot.lipidIDs_sel = plot.lipidIDs
# plot.lipidIDs_sel = plot.lipidIDs[1:50]
# plot.lipidIDs_unsel = plot.lipidIDs[-c(1:50)]#0

p5 <-generate_forest_col(
  data = df_intervals %>% 
    filter(cohort %in% plot.mycohorts) %>%
    filter(lipidID %in% plot.lipidIDs_sel));p5

## 5f. combine plots ----
plot.lims = df_intervals %>% 
  filter(cohort %in% plot.mycohorts) %>%
  select(.xmin,.xmax) %>% range

final_grid <- p1 + p2 + p3 + p4 + p5 +
  plot_layout(ncol = 5, widths = c(1, 1, 1, 1), axes = 'collect', guides='collect') &
  facet_grid(rows = vars(lipidClass), scales = "free_y", space = "free_y") &
  theme(axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size=8),
        axis.text.y = element_text(size = 6.125),
        strip.text.y = element_blank(),
        # strip.text.y = element_text(angle = 90, face = "bold", size=5),
        # 2. Strip away outer plot margins so patchwork columns sit tightly together
        plot.margin = margin(t = 5, r = 2.5, b = 5, l = 2.5)) &
  # limits
  scale_shape_manual(values=c('p.fdr<0.05'=19,'p.fdr>=0.05'=21)) &
  scale_fill_manual(values=c('p.fdr<0.05'='black','p.fdr>=0.05'='white')) &
  coord_cartesian(xlim=plot.lims) &
  xlab("beta (SD units) +/- 95% CI") &
  ylab(NULL);final_grid

ggsave(
  filename = paste0("forestplot_252lipids_PD_and_Controls_nonsevere_gba1eff_unadjBMI_",Sys.Date(),".pdf"),
  plot = final_grid,
  width = 33.15,
  height = 15.34,
  units = "cm",
  device = cairo_pdf
)

## file with legend
final_grid_legend = final_grid & 
  guides(
    color = guide_legend(
      reverse = TRUE,
      ncol = 1,
      override.aes = list(size = 1)
    ),
    shape = guide_legend(
      reverse = FALSE,
      ncol = 1,
      override.aes = list(size = 1)
    ),
    alpha = guide_legend(
      reverse = FALSE,
      ncol = 1,
      override.aes = list(size = 1)
    )
  ) & 
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.height = unit(0.1, "cm"),
    legend.key.spacing.y = unit(0.05, "cm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  );final_grid_legend

ggsave(
  filename = paste0("forestplot_252lipids_PD_and_Controls_nonsevere_gba1eff_unadjBMI_wiLegend_",Sys.Date(),".pdf"),
  plot = final_grid_legend,
  width = 33.8,
  height = 17.8,
  units = "cm",
  device = cairo_pdf
)

## ends ----
