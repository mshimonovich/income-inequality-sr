library(tidyverse)
library(mvmeta)
library(readxl)
library(SciViews)
library(naniar)
library(meta)
library(metafor)
library(ggplot2)
library(grid)
library(rlang)
library(data.table)
library(stringi)
library(reshape2)
library(ggtext)
##use main for list
setwd("C:/Users/mshim/OneDrive/Documents/r_projects/incominequal_sr")

rob_tbl <- readRDS("cleaning/complete.rds")

rob_tbl <- map(rob_tbl, function(x) {x %>% filter(main==1) %>%  mutate(study=str_glue("{author_comb}<sup>{ref_id}</sup>"),
             across(matches('rob'), ~case_when(.=="3"|.=="critical"|.=="Critical"~"Critical",
                                               .=="2"|.=="serious"|.=="Serious"~"Serious",
                                               .=="1"|.=="m"|.=="moderate"|.=="Moderate"~"Moderate",
                                               .=="0"|.=="low"|.=="Low"~"Low",
                                                .=="4"|.=="no information"|.=="No information"~"No information",
                                               TRUE~"."))) %>%  group_by(ref_id) %>% slice(n=1) %>% ungroup()  %>% select(study, starts_with("rob_"), rob)
                                                         })
# converting data to long form for ggplot2 use
rob_long <- map(rob_tbl, function(x) {x %>% melt(., id.var='study') %>% mutate(value=factor(value, levels=c("Low", "Moderate", "Serious", "Critical", "NR"))
                                                                       )})
domain_long <- c("Bias due to confounding", "Bias in selection of participants into the study", "Bias in classification of interventions", "Bias due to deviations from intended interventions", "Bias due to missing data", "Bias in measurement of the outcome", "Bias in selection of the reported result", "Overall Risk of Bias.")
domain_short <- c("Confounding", "Selection", "Intervention classification", "Intervention deviations", "Missing data", "Outcome measurement", "Reporting results", "Overall")
domain_legend <- paste0(domain_short, " = ", domain_long)
domain_caption <- paste(domain_legend, collapse="; ")
rob_val <- c("green3", "yellow1", "orange", "red3", "gray70")                        

robfig_default <- list(labs(caption = paste0("<b>Response options for risk of bias:</b> <br>
<b><span style = 'color:green3;'>Low =</span></b> The study is comparable to a well-performed randomized trial.<br>
<b><span style = 'color:yellow3;'>Moderate =</span></b> The study provides sound evidence for a nonrandomized study but cannot be
considered comparable to a well-performed randomized trial.<br>
<b><span style = 'color:orange2;'>Serious =</span></b> The study has important problem(s).<br>
<b><span style = 'color:red3;'>Critical =</span></b> The study is too problematic to provide any useful evidence.<br>
<b><span style = 'color:gray48;'>No information:</span></b> The study does not provide sufficient information to a make judgement about risk of bias.<br>", domain_caption)), theme_bw(), labs(x="Risk of Bias Domains", y=""),
        theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
              plot.caption.position="plot", 
              plot.caption = element_textbox_simple(size = 10, padding = margin(7, 4, 4, 4),  margin = margin(7, 0, 0, 0)),
              legend.position="none",
              panel.border = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.grid.major = element_line(colour="white", linewidth=1.5),
              axis.line = element_blank(), 
              axis.text.x.top = element_text(vjust = .1, hjust=.2, angle = 49, size=12, colour="black"),
              axis.text.y = element_markdown(colour="black", size=12), 
              axis.ticks = element_blank(),
              axis.title.x.top = element_text(size=14, colour="black", face="bold")
              ),scale_x_discrete(labels = str_wrap(domain_short, width = 1.3), position = "top"), scale_color_manual(values= rob_val))
                   
robfig_srh <- ggplot(data=rob_long$srh,aes(x=variable, y=fct_rev(study))
                ) + geom_point(color='black', shape=21, size=4.3, aes(fill=factor(value))
                                ) + theme(aspect.ratio = 10 /8.3) + scale_fill_manual(values=rob_val) + robfig_default + coord_fixed(ratio = 3/10) 

robfig_srh

dev.copy(pdf,
         file = paste0("outputs/rob_srh.pdf"),
         width = 8.0,
         height = 10)
dev.off()

robfig_acm <- ggplot(data=rob_long$acm,aes(x=variable,  y=fct_rev(study))) + geom_point(aes(color=value), size = 5) + robfig_default
robfig_acm <-ggplot(data=rob_long$acm,aes(x=variable, y=fct_rev(study))  
       ) + geom_point(color='black', shape=21, size=4.3, aes(fill=factor(value))
                      ) + theme(aspect.ratio = 7 /10, axis.title.y=element_text(lineheight = 1.2)
                                ) + scale_fill_manual(values=rob_val) + robfig_default + coord_fixed(ratio = 5.2/10) 

robfig_acm

dev.copy(pdf,
         file = paste0("outputs/rob_acm.pdf"),
         width = 7.8,
         height = 9)
dev.off()