################Mortality

table_acm <- main$acm %>% fill(subsamp_pos, .direction="up"
) %>% group_by(group_id) %>% ungroup() %>% mutate(row_order=rev(row_number()),
                                                  row_space=case_when(rob_id==1~row_order+4, #####to make it easy to get lengths of each subgroup of risk of bias
                                                                      rob_id==2~row_order+0),
                                                  row_expand=ifelse(bigstudy==TRUE, 0.3, 0), 
                                                  row_expand=lead(row_expand, default=0), 
                                                  row_expand=rev(cumsum(rev(row_expand))), 
                                                  row_fp=row_space+row_expand,
                                                  row_fp=row_fp-min(row_fp)) ###zero the start so everything goes around that

### Label with superscript references
label_acm <- mapply(function(x,y) as.expression(bquote(.(x)^.(y))), table_acm$sample_lbl, table_acm$ref_id)
res_acm <- rma(yi, sei=si, data=table_acm, method="REML", drop00 = TRUE)
options(OutDec=".")
### weights and format them for text
weights_acm <- paste0(formatC(weights(res_acm), format="f", digits=1), "%")
weights_acm[weights_acm == "NA%"|is.na(weights_acm)] <- ""
### Point sizes 
psize_acm <- weights(res_acm)
psize_acm[is.na(psize_acm)]<-0
psize_acm <- 0.6 + (psize_acm - min(psize_acm, na.rm=TRUE)) / (max(psize_acm, na.rm = TRUE) - min(psize_acm, na.rm=TRUE))

### Margins

par(mar=c(bottom=6, left=0.5, top=3, right=1.5), mgp=c(3,0,0), tcl=0.15)
options(na.action = "na.pass")

fp_acm<-forest(res_acm, xlim=c(-1.4, 2.85), at=c(0.75, 1.0, 1.25, 1.50, 1.75), ylim=c(min(table_acm$row_fp)-3, max(table_acm$row_fp)),addfit = FALSE,
               transf=exp, refline=NA,lty=c(1,0,0),xlab="", slab=label_acm, ilab=cbind(table_acm$setting, weights_acm), 
               ilab.xpos=c(0.1, 1.9),ilab.pos=c(4,4), textpos=c(-1.4, 2.85), psize=psize_acm, pch=18, cex=0.85, efac=c(0,0,0), 
               annosym=c(" (", "–", ")"), mlab="", rows = table_acm$row_fp, rowadj=-0.05)


### add vertical reference line at 1
segments(1, fp_acm$ylim[1]-2.1, 1, max(fp_acm$row)+0.5, col=coll)
### add vertical reference line at the pooled estimate
segments(exp(coef(res_acm)), fp_acm$ylim[1]-2.1, exp(coef(res_acm)), max(fp_acm$row)+0.5, col=colp, lty="33", lwd=0.8)

### add horizontal line at the top
abline(h=max(fp_acm$row)+2.0, col=coll)
### redraw the x-axis in the chosen color
axis(side=1, at=seq(0.75, 1.75,by=0.25), col=coll, labels=FALSE)
axis(side=1, at=seq(1, 1.5,by=0.25), col=coll, col.ticks = coll, labels=FALSE)

## add rounded rectangle
#roundedRect(xleft=min(fp_acm$at),ybottom=fp_acm$ylim[1]-2.3, xright=max(fp_acm$at), ytop=max(fp_acm$row)+0.45, rounding=0.35, border=coll)
par(xpd=NA)

### redraw the CI lines and points in the chosen color
segments(table_acm$lcl, table_acm$row_fp, table_acm$ucl, table_acm$row_fp, col=colp, lwd=1.75)
points(table_acm$est, table_acm$row_fp, pch=18, cex=psize_acm*1.15, col="white")
points(table_acm$est, table_acm$row_fp, pch=18, cex=psize_acm, col=colp)


### adjust cex as used in the forest plot and use a bold font
par(cex=fp_acm$cex, font=2)
text(c(fp_acm$ilab.xpos[1]-0.05,fp_acm$ilab.xpos[2]-0.15), fp_acm$ylim[2]+2.8, c("Location", "Weight (%)"), cex=fp_acm$cex*1.2, pos=c(4, 4))
text(c(fp_acm$textpos[1], fp_acm$textpos[2]+0.07), fp_acm$ylim[2]+2.8, c("Study name, year", "Risk ratio (95% CI)"), cex=fp_acm$cex*1.2, pos=c(4, 2))

text(fp_acm$xlim[1], c(max(table_acm$row_fp[which(table_acm$rob_id==1)])+1.2, 
                       max(table_acm$row_fp[which(table_acm$rob_id==2)])+1.2), pos=4, c("Moderate risk of bias",
                                                                                        "Serious risk of bias"), cex=fp_acm$cex*1.1)

### fit random-effects model in the three subgroups
res.m_acm <- rma(yi, sei=si, subset=(rob_id==1), data=table_acm, drop00 = TRUE)
res.s_acm <- rma(yi, sei=si, subset=(rob_id==2), data=table_acm, drop00 = TRUE)

par(font=3, col=coll)

### add summary polygons for  subgroups and total
addpoly(res.m_acm, row= min(table_acm$row_fp[which(table_acm$rob_id==1)])-1.1, mlab="Random effects model for subgroup", efac=1, col=colp, border=colp,cex=fp_acm$cex*1.0)
text(x=fp_acm$ilab.xpos[2], y=min(table_acm$row_fp[which(table_acm$rob_id==1)])-1.1, labels=paste0(formatC((meta$acm$w.random.w[1]/sum(meta$acm$w.random.w))*100, digits=1,  format="f"), "%"), pos=4, cex=fp_acm$cex*1.0)
text(fp_acm$xlim[1], min(table_acm$row_fp[which(table_acm$rob_id==1)])-2.1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                                                              .(formatC(meta$acm$tau2.w[1], digits=2, format="f")), "; ", chi^2, "=",
                                                                                              .(formatC(meta$acm$Q.w[1], digits=2, format="f")), "; df=", .(meta$acm$df.predict.w[1]),
                                                                                              ", P=", .(formatC(meta$acm$pval.Q.w[1], digits=2, format="f")), "; ", I^2, "=",
                                                                                              .(formatC(meta$acm$I2.w[1]*100, digits=1, format="f")), "%")), cex=fp_acm$cex*1)

addpoly(res.s_acm, row= min(table_acm$row_fp[which(table_acm$rob_id==2)])-1, mlab="Random effects model for subgroup", efac=1, col=colp, border=colp,cex=fp_acm$cex*1.0)
text(x=fp_acm$ilab.xpos[2], y=min(table_acm$row_fp[which(table_acm$rob_id==2)])-1, labels=paste0(formatC((meta$acm$w.random.w[2]/sum(meta$acm$w.random.w))*100, digits=1, format="f"), "%"), pos=4,cex=fp_acm$cex*1.0)
text(fp_acm$xlim[1], min(table_acm$row_fp[which(table_acm$rob_id==2)])-2, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                                                              .(formatC(meta$acm$tau2.w[2], digits=2, format="f")), "; ", chi^2, "=",
                                                                                              .(formatC(meta$acm$Q.w[2], digits=2, format="f")), "; df=", .(meta$acm$df.predict.w[2]),
                                                                                              ", P=", .(formatC(meta$acm$pval.Q.w[2], digits=2, format="f")), "; ", I^2, "=",
                                                                                              .(formatC(meta$acm$I2.w[2]*100, digits=1, format="f")), "%")), cex=fp_acm$cex*1)

#######total
par(font=1, col="black")

addpoly.default(x=meta$acm$TE.random, sei=meta$acm$seTE.random, row= min(fp_acm$rows)-3.8, mlab="", efac=1, border=colp,col=colp, annotate=FALSE)
text(x=c(fp_acm$xlim[1], fp_acm$ilab.xpos[2]-0.05, fp_acm$textpos[2]+0.03), y=min(fp_acm$rows)-3.8, labels=c("Overall", "100.0%", "1.02 (1.00–1.04)"),cex=fp_acm$cex*1.1, pos=c(4,4, 2), font=2)

### test for heterogeneity
text(fp_acm$xlim[1], min(fp_acm$rows)-5.1, pos=4, bquote(paste("Test for heterogeneity: ", tau^2, "=",
                                                                 .(formatC(meta$acm$tau2, digits=2, format="f")), "; ", chi^2, "=",
                                                                 .(formatC(meta$acm$Q, digits=2, format="f")))), cex=fp_acm$cex*1.1, font=1)

text(fp_acm$xlim[1], min(fp_acm$rows)-6.4, pos=4, bquote(paste("df=", .(meta$acm$df.Q),
                                                                 ", P=", .(formatC(meta$acm$pval.Q, digits=2, format="f")), "; ", I^2, "=",
                                                                 .(formatC(meta$acm$I2*100, digits=1, format="f")), "%")), cex=fp_acm$cex*1.1, font=1)


### add text for test of overall effect
text(fp_acm$xlim[1],min(fp_acm$rows)-7.7, pos=4, bquote(paste("Test for overall effect: Z=",
                                                                  .(formatC(meta$acm$zval.random, digits=2, format="f")),
                                                                  ", P=", .(formatC(meta$acm$pval.random, digits=2, format="f")))), cex=fp_acm$cex*1.1)
### add text for subgroup differences
text(fp_acm$xlim[1],min(fp_acm$rows)-9.0, pos=4, bquote(paste("Test for subgroup differences: Q=",
                                                                .(formatC(meta$acm$Q.b.random, digits=2, format="f")))), cex=fp_acm$cex*1.1)
                                                              
text(fp_acm$xlim[1],min(fp_acm$rows)-10.3, pos=4, bquote(paste("df=", .(meta$acm$df.Q.b.random),
                                                                ", P=", .(formatC(meta$acm$pval.Q.b.random, digits=2, format="f")))), cex=fp_acm$cex*1.1)


left<-stringr::str_wrap("Increase in Gini coefficient improves health", 30, exdent=2)
right<- stringr::str_wrap("Increase in Gini coefficient worsens health", 30, exdent=2)

text(c(fp_acm$at[3],fp_acm$at[3]), min(fp_acm$rows)-8.1, c("Lower mortality","Higher mortality"), pos=c(2,4), offset=0.3, font=2, cex=fp_acm$cex*1.0)

text(c(fp_acm$at[3],fp_acm$at[3]), min(fp_acm$rows)-10.4, c(left,right), pos=c(2,4), offset=0.3, font=2, cex=fp_acm$cex*1.0)
Arrows(x0=fp_acm$at[3]+0.02, y0=min(fp_acm$rows)-8.8, x1=fp_acm$at[3]+0.8, y1=min(fp_acm$rows)-8.8,  arr.type = "curved", arr.length = 0.15,arr.width = 0.05, arr.adj = 1)
Arrows(x0=fp_acm$at[3]-0.02, y0=min(fp_acm$rows)-8.8, x1=fp_acm$at[3]-0.8, y1=min(fp_acm$rows)-8.8,  arr.type = "curved", arr.length = 0.15,arr.width = 0.05, arr.adj = 1)


dev.copy(pdf,
         file = paste0("outputs/foresplot_acm.pdf"),
         width = 8.0,
         height = 10)
dev.off()

