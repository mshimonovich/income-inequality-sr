
analysis <- readRDS("cleaning/complete_main.rds")
meta <- readRDS("analysing/main_meta.rds")

meta<-map(list(meta$srh$main, meta$acm$main), ~update(., subgroup=rob))
names(meta) <- list("srh", "acm")

##function creating dataframes for foresplot, want to make a newone so that main can be accessed for other things (like MR)
main <- map(analysis, ~.x %>% select(matches(c('study_|sample_|id')), setting, yi, si, rob) %>% arrange(rob)  %>% group_by(covidence_id) %>% mutate(
    group_dup=cur_group_id()) %>% ungroup() %>% mutate(
      group_id=rleid(group_dup)) %>% group_by(rob) %>% 
    mutate(robgrp_dup=cur_group_id()) %>% ungroup() %>% mutate(rob_id=rleid(robgrp_dup)
    ) %>% select(-c(group_dup, robgrp_dup)) %>% group_by(ref_id) %>% mutate(subsamp_len=n()
    ) %>% ungroup()%>% mutate(subsamp_dif=c(0, diff(subsamp_len)
    ), subsamp_cum=cumsum(subsamp_dif), subsamp_pos=rleid(subsamp_cum)
    ) %>% select(-c(subsamp_dif, subsamp_len, subsamp_cum)) %>% mutate(studlab=row_number()))

main <- map(main, ~ .x %>% group_split(subsamp_pos) %>%  map_dfr(~ .x %>%  add_row(.before = 1)
) %>% fill(c(covidence_id:study_desc, rob_id, group_id), .direction="up") %>% ungroup()
%>% mutate(condition=ifelse(sample_desc=="whole sample" & is.na(yi), "drop", "keep")) %>% filter(condition=="keep") %>% select(-c(condition)) %>% 
  mutate(sample_desc=str_to_sentence(sample_desc),
         sample_lbl = case_when(
    is.na(yi) & is.na(subsamp_pos) ~ paste0(" ", study_id, ", ", study_year), 
    !is.na(yi) & sample_desc=="Whole sample"~ paste0(" ", study_id, ", ", study_year),
    !is.na(yi) & sample_desc!="Whole sample"~paste0("    ", sample_desc)),
    est=exp(yi), 
    lcl=exp(yi-(1.96*si)), 
    ucl=exp(yi+(1.96*si)),
    bigstudy = case_when(
      is.na(yi) & is.na(subsamp_pos) ~ TRUE, 
      !is.na(yi) & sample_desc=="Whole sample" ~TRUE,
      !is.na(yi) & sample_desc!="Whole sample"~FALSE
    ),
    ref_id= ifelse(bigstudy=="FALSE", "", ref_id)
  ))

colp <- "#B1BAC1"
coll <- "#232734"
    
####make forestplot for SRH and mortality
table_srh <- main$srh %>% fill(subsamp_pos, .direction="up"
                               ) %>% group_by(group_id) %>% ungroup() %>% mutate(row_order=rev(row_number()),
                                                                    row_space=case_when(rob_id==1~row_order+8, #####to make it easy to get lengths of each subgroup of risk of bias
                                                                              rob_id==2~row_order+4,
                                                                              rob_id==3~row_order+0),
                                                                    row_expand=ifelse(bigstudy==TRUE, 0.2, 0), 
                                                                    row_expand=lead(row_expand, default=0), ###zero the start so everything goes around that
                                                                    row_expand=rev(cumsum(rev(row_expand))), 
                                                                    row_fp=row_space+row_expand,
                                                                    row_fp=row_fp-min(row_fp))


###label with references at superscripts
label_srh <- mapply(function(x,y) as.expression(bquote(.(x)^.(y))), table_srh$sample_lbl, table_srh$ref_id)
res_srh <- rma(yi, sei=si, data=table_srh, method="REML", drop00 = TRUE)
options(OutDec=".")
### weights and format them for text
weights_srh <- paste0(formatC(weights(res_srh), format="f", digits=1), "%")
weights_srh[weights_srh == "NA%"|is.na(weights_srh)] <- ""
### Point sizes 
psize_srh <- weights(res_srh)
psize_srh <- 0.6 + (psize_srh - min(psize_srh, na.rm=TRUE)) / (max(psize_srh, na.rm = TRUE) - min(psize_srh, na.rm=TRUE))
psize_srh[is.na(psize_srh)]<-0
#windows(height = 23, width = 20)
### Margins
par(mar=c(bottom=6, left=0.5, top=2, right=1.3), mgp=c(3,0,0), tcl=0.15)
options(na.action = "na.pass") 



fp_srh<-forest(res_srh, xlim=c(-1.7, 2.8), at=c(0.5, 0.75, 1, 1.25, 1.5), ylim=c(min(table_srh$row_fp)-3, max(table_srh$row_fp)),addfit = FALSE,
       transf=exp, refline=NA,lty=c(1,0,0),xlab="", slab=label_srh, ilab=cbind(table_srh$setting, weights_srh), 
       ilab.xpos=c(-0.53, 1.70),ilab.pos=c(4,4), textpos=c(-1.7, 2.7), psize=psize_srh, pch=18, cex=0.84, efac=c(0,0,0), 
       annosym=c(" (", "–", ")"), mlab="", rows = table_srh$row_fp, rowadj=-0.07)

### add vertical reference line at 1
segments(1, fp_srh$ylim[1]-2.6, 1, max(fp_srh$row)+0.5, col=coll)
### add vertical reference line at the pooled estimate
segments(exp(coef(res_srh)), fp_srh$ylim[1]-1.6, exp(coef(res_srh)), max(fp_srh$row)+0.5, col=colp, lty="33", lwd=0.8)

### add horizontal line at the top
abline(h=max(fp_srh$row)+2, col=coll)
### redraw the x-axis in the chosen color
axis(side=1, at=seq(0.5, 1.5,by=0.25), col=coll, labels=FALSE)
axis(side=1, at=seq(0.75, 1.25,by=0.25), col=coll, col.ticks = coll, labels=FALSE)
## add rounded rectangle
#roundedRect(xleft=min(fp_srh$at), ybottom=fp_srh$ylim[1]-2.6, xright=max(fp_srh$at), ytop=max(fp_srh$row)+0.5, rounding=0.35, border=coll)
par(xpd=NA)

### redraw the CI lines and points in the chosen color
table_srh$ucl[which(table_srh$ucl>1.5)]<-1.5
segments(table_srh$lcl, table_srh$row_fp, table_srh$ucl, table_srh$row_fp, col=colp, lwd=1.75)
Arrowhead(x0=c(1.5, 1.5),y0=c(table_srh$row_fp[which(table_srh$ucl==1.5)]), arr.type = "triangle",lcol=colp, arr.col=colp,  arr.length = 0.1,arr.width = 0.04, arr.adj = 1)
points(table_srh$est, table_srh$row_fp, pch=18, cex=psize_srh*1.15, col="white")
points(table_srh$est, table_srh$row_fp, pch=18, cex=psize_srh, col=colp)

### adjust cex as used in the forest plot and use a bold font
par(cex=fp_srh$cex, font=2)
text(c(fp_srh$ilab.xpos[1],fp_srh$ilab.xpos[2]-0.15), fp_srh$ylim[2]+2.8, c("Location", "Weight (%)"), cex=fp_srh$cex*1.2, pos=c(4, 4))
text(c(fp_srh$textpos[1], fp_srh$xlim[2]+0.05), fp_srh$ylim[2]+2.8, c("Study name, year", "Odds ratio (95% CI)"), cex=fp_srh$cex*1.2, pos=c(4, 2))

text(fp_srh$xlim[1], c(max(table_srh$row_fp[which(table_srh$rob_id==1)])+1.2, 
                                        max(table_srh$row_fp[which(table_srh$rob_id==2)])+1.2,
                       max(table_srh$row_fp[which(table_srh$rob_id==3)])+1.2), pos=4, c("Moderate risk of bias",
                                                                               "Serious risk of bias",
                                                                               "Critical risk of bias"), cex=fp_srh$cex*1.1)

### fit random-effects model in the three subgroups
res.m_srh <- rma(yi, sei=si, subset=(rob_id==1), data=table_srh, drop00 = TRUE)
res.s_srh <- rma(yi, sei=si, subset=(rob_id==2), data=table_srh, drop00 = TRUE)
res.c_srh <- rma(yi, sei=si, subset=(rob_id==3), data=table_srh, drop00 = TRUE)

par(font=3, col=coll)
### add summary polygons for  subgroups and total
addpoly(res.m_srh, row= min(table_srh$row_fp[which(table_srh$rob_id==1)])-1.1, mlab="", efac=1, col=colp, border=colp)
text(x=c(fp_srh$xlim[1], fp_srh$ilab.xpos[2]), 
     y=min(table_srh$row_fp[which(table_srh$rob_id==1)])-1.1, 
            labels=c("Random effects model for subgroup", paste0(formatC((meta$srh$w.random.w[1]/sum(meta$srh$w.random.w))*100, digits=1,  format="f"), "%")), pos=4, cex=fp_srh$cex*1.0)
text(fp_srh$xlim[1], min(table_srh$row_fp[which(table_srh$rob_id==1)])-2.1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                               .(formatC(meta$srh$tau2.w[1], digits=2, format="f")), "; ", chi^2, "=",
                                                               .(formatC(meta$srh$Q.w[1], digits=2, format="f")), "; df=", .(meta$srh$df.predict.w[1]),
                                                               ", P=", .(formatC(meta$srh$pval.Q.w[1], digits=2, format="f")), "; ", I^2, "=",
                                                               .(formatC(meta$srh$I2.w[1]*100, digits=1, format="f")), "%")), cex=fp_srh$cex*0.9)

addpoly(res.s_srh, row= min(table_srh$row_fp[which(table_srh$rob_id==2)])-1.1, mlab="", efac=1, col=colp, border=colp,cex=fp_srh$cex*1.0)
text(x=c(fp_srh$xlim[1], fp_srh$ilab.xpos[2]), 
     y=min(table_srh$row_fp[which(table_srh$rob_id==2)])-1.1, 
     labels=c("Random effects model for subgroup", paste0(formatC((meta$srh$w.random.w[2]/sum(meta$srh$w.random.w))*100, digits=1, format="f"), "%")), pos=4,cex=fp_srh$cex*1.0)
text(fp_srh$xlim[1], min(table_srh$row_fp[which(table_srh$rob_id==2)])-2.1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                               .(formatC(meta$srh$tau2.w[2], digits=2, format="f")), "; ", chi^2, "=",
                                                               .(formatC(meta$srh$Q.w[2], digits=2, format="f")), "; df=", .(meta$srh$df.predict.w[2]),
                                                               ", P=", .(formatC(meta$srh$pval.Q.w[2], digits=2, format="f")), "; ", I^2, "=",
                                                               .(formatC(meta$srh$I2.w[2]*100, digits=1, format="f")), "%")), cex=fp_srh$cex*0.9)

addpoly(res.c_srh, row=min(table_srh$row_fp[which(table_srh$rob_id==3)])-1.1, mlab="", efac=1, col=colp, border=colp,cex=fp_srh$cex*1.0)
text(x=c(fp_srh$xlim[1], fp_srh$ilab.xpos[2]),  
     y=min(table_srh$row_fp[which(table_srh$rob_id==3)])-1.1, 
     labels=c("Random effects model for subgroup", paste0(formatC((meta$srh$w.random.w[3]/sum(meta$srh$w.random.w))*100, digits=1, format="f"), "%")), pos=4,cex=fp_srh$cex*1.0)
text(fp_srh$xlim[1], min(table_srh$row_fp[which(table_srh$rob_id==3)])-2.1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                               .(formatC(meta$srh$tau2.w[3], digits=2, format="f")), "; ", chi^2, "=",
                                                               .(formatC(meta$srh$Q.w[3], digits=2, format="f")), "; df=", .(meta$srh$df.predict.w[3]),
                                                               ", P=", .(formatC(meta$srh$pval.Q.w[3], digits=2, format="f")), "; ", I^2, "=",
                                                               .(formatC(meta$srh$I2.w[3]*100, digits=1, format="f")), "%")), cex=fp_srh$cex*0.9)

#######total
overall<-par(font=1, col="black")
overall
addpoly.default(x=meta$srh$TE.random, sei=meta$srh$seTE.random, row= min(fp_srh$rows)-3.6, mlab="", annotate=FALSE, efac=1, col=colp, border=colp)
text(x=c(fp_srh$xlim[1], fp_srh$ilab.xpos[2]-0.03, fp_srh$textpos[2]+0.03), y=min(fp_srh$rows)-3.6, labels=c("Overall", "100.0%", "1.06 (1.03–1.08)"),cex=c(fp_srh$cex*1.2, fp_srh$cex*1.1, fp_srh$cex*1.1), pos=c(4,4, 2), font=2)

### test for heterogeneity
text(fp_srh$xlim[1], y= min(fp_srh$rows)-4.9, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                                 .(formatC(meta$srh$tau2, digits=2, format="f")), "; ", chi^2, "=",
                                                                 .(formatC(meta$srh$Q, digits=2, format="f")),";")), cex=fp_srh$cex*1.1)
text(fp_srh$xlim[1], y= min(fp_srh$rows)-6.2, pos=4, bquote(paste("df=", .(meta$srh$df.Q),
                                                                 ", P=", .(formatC(meta$srh$pval.Q, digits=2, format="f")), "; ", I^2, "=",
                                                                 .(formatC(meta$srh$I2*100, digits=1, format="f")), "%")), cex=fp_srh$cex*1.1)

### add text for test of overall effect
text(fp_srh$xlim[1],min(fp_srh$rows)-7.5, pos=4, bquote(paste("Test for overall effect: Z=",
                                                                .(formatC(meta$srh$zval.random, digits=2, format="f")),
                                                                ", P=", .(formatC(meta$srh$pval.random, digits=2, format="f")))), cex=fp_srh$cex*1.1)
### add text for subgroup differences
text(fp_srh$xlim[1],min(fp_srh$rows)-8.8, pos=4, bquote(paste("Test for subgroup differences: Q=",
                                                                .(formatC(meta$srh$Q.b.random, digits=2, format="f")), ";")))
text(fp_srh$xlim[1],min(fp_srh$rows)-10.1, pos=4, bquote(paste("df=", .(meta$srh$df.Q.b.random),
                                                                ", P=", .(formatC(meta$srh$pval.Q.b.random, digits=2, format="f")))), cex=fp_srh$cex*1.1)

        
        
left<-stringr::str_wrap("Increase in Gini coefficient improves health", 30, exdent=2)
right<- stringr::str_wrap("Increase in Gini coefficient worsens health", 30, exdent=2)

text(c(fp_srh$at[3],fp_srh$at[3]), min(fp_srh$rows)-8.8, c("Good self-rated health","Poor self-rated health"), pos=c(2,4), offset=0.3, font=2, cex=fp_srh$cex*1.1)
text(c(fp_srh$at[3],fp_srh$at[3]), min(fp_srh$rows)-11.3, c(left,right), pos=c(2,4), offset=0.3, font=2, cex=fp_srh$cex*1.0)
Arrows(x0=fp_srh$at[3]+0.02, y0=min(fp_srh$rows)-9.7, x1=fp_srh$at[3]+1, y1=min(fp_srh$rows)-9.7,  arr.type = "curved", arr.length = 0.15,arr.width = 0.05, arr.adj = 1)
Arrows(x0=fp_srh$at[3]-0.02, y0=min(fp_srh$rows)-9.7, x1=fp_srh$at[3]-1, y1=min(fp_srh$rows)-9.7,  arr.type = "curved", arr.length = 0.15,arr.width = 0.05, arr.adj = 1)


dev.copy(pdf,
         file = paste0("outputs/foresplot_srh.pdf"),
         width = 8.0,
         height = 10)
dev.off()


########dpme








