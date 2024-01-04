##create data file for conducting main meta-analysis and meta-regressions
analysis <- readRDS("cleaning/complete_main.rds")

#############Meta-analysis
datasets_or <-list()
datasets_or$main <-analysis$srh
models_or <- map(datasets_or,
                 ~metagen(TE=yi,
                          seTE=si,
                          data = .,
                          comb.fixed = FALSE,
                          comb.random = TRUE,
                          method.tau = "REML",
                          hakn = FALSE,
                          prediction = FALSE,
                          sm = "OR",
                          n.e = sample_size)
)

datasets_rr<-list()
datasets_rr$main <-analysis$acm
models_rr <- map(datasets_rr, 
                 ~metagen(TE=yi,
                          seTE=si,
                          data = .,
                          comb.fixed = FALSE,
                          comb.random = TRUE,
                          method.tau = "REML",
                          hakn = FALSE,
                          prediction = FALSE,
                          sm = "RR",
                          n.e = sample_size)
)
saveRDS(list(srh = models_or, acm=models_rr), file = "analysing/main_meta.rds")

#######create file for meta-regression
metareg.func<- function(model) {
  output <- list()
  output$rob <- metareg(model, ~ rob)
  output$gscale_cat <- metareg(model, ~gscale_cat)
  output$gini_cat <- metareg(model, ~ gini_cat)
  output$region_cat <- metareg(model, ~ region_cat)
  output$usa_cat <- metareg(model, ~ usa_cat)
  output$age_cat <- metareg(model, ~age_cat)
  output$adjust_cat <- metareg(model, ~ adjust_cat)
  output$fu_cat <- metareg(model, ~fu_cat) 
  output$timelag_cat <- metareg(model, ~timelag_cat)
  return(output)
}


modelreg <- list()
modelreg$srh <- map(models_or, metareg.func)
modelreg$acm <- map(models_rr, metareg.func)
saveRDS(object=analysis, file = "analysing/meta_regress.rds")
