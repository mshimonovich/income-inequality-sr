pacman::p_load(ggplot2, tidyverse, dplyr, dosresmeta, readxl, SciViews, naniar, meta, mvmeta, metafor, ggplot2, grid, rlang, data.table, stringi, gt, glue, shape, berryFunctions)
##1. Create list for all-cause mortality (acm) and self-rated health (srh) dataframes and transform estimates to logform 
dflist<-list()
dflist$srh <- read_excel("cleaning/prep_srh.xlsx", na = "NA")
dflist$acm <- read_excel("cleaning/prep_acm.xlsx", na = "NA")
dflist$acm<-dflist$acm %>% mutate(rat_re=ifelse(out_est=="logit" & is.na(rat_re), exp(b_re), rat_re),
                    se_re=ifelse(is.na(se_re) & out_est=="HR", b_re/tscore_re, abs(se_re)),            
               lcl_re=ifelse(out_est=="logit" & is.na(lcl_re), exp(b_re-(qnorm(1-0.05/2)*se_re)), lcl_re),
               ucl_re=ifelse(out_est=="logit" & is.na(ucl_re), exp(b_re+(qnorm(1-0.05/2)*se_re)), ucl_re)) %>% 
            mutate(across(c(rat_re, lcl_re, ucl_re),
                   ~case_when(out_est=="logit" | out_est == "OR"~ (.)/((1-rate_lev)+(rate_lev*.)),##OR to RR
                                 out_est=="HR"~ (1-exp(.*ln(1-rate_lev)))/1, ##HR to RR
               TRUE~.),  .names = "{.col}_new"))
dflist$srh <- dflist$srh %>% mutate(across(c(se_re, b_re), ~ifelse(out_est=="probit" & !is.na(.), .*1.81, .)))


#2. transforming estimates of categorical exposures to continuous
##generating id for each group of estimates (needed for categorical exposures)
dflist1<-map(dflist, function(x) {
  x <- x %>% mutate(
       b=ifelse(!is.na(rat_re), ln(rat_re), b_re), 
            se=if_else(!is.na(ucl_re), ((ln(ucl_re)-ln(lcl_re))/(2 * qnorm(1 - 0.05/2))), se_re, NA_real_))  %>% 
    group_by(covidence_id, author, study_id, sample_desc, gscale_id, timelag_n, adjust_cat, dummy_cat) %>%
    mutate(id_dup=cur_group_id()) %>% ungroup() %>% mutate(id=rleid(id_dup))
}
)
###reverse code if SRH upper and exposure type is not Median share
dflist1$srh<-dflist1$srh %>% 
  mutate(across(c(b, se), 
                ~ ifelse(exp_meas!="median share" & out_meas=="upper", .x*-1, .x),
                .names = "{.col}")) %>% mutate(se=abs(se))

# getting the dose for categorical exposures where a range for each level were report
dflist2<-map(dflist1, function(x){ 
  x %>% mutate (min="\\d+\\.*\\d*\\s?(?=-)|(?<=>)\\s?\\d+\\.*\\d*", 
                   max= "(?<=-)\\s?\\d+\\.*\\d*|(?<=<)\\s?\\d+\\.*\\d*", 
                   dose_min=as.numeric(str_extract(dose_re, min)), 
                   dose_max=as.numeric(str_extract(dose_re, max))) %>% 
  group_by(id) %>% mutate(SD=sd(unlist(select(cur_data(), dose_min:dose_max)), na.rm=TRUE),
                   dose_min=if_else(is.na(dose_min), dose_max-(2*SD), dose_min),
                   dose_max=if_else(is.na(dose_max), dose_min+(2*SD), dose_max),
                   dose=if_else(!is.na(SD), map2_dbl(dose_min, dose_max, ~mean(c(.x, .y))), as.numeric(dose_re)),
                   dose0=case_when(               
                          first(dose)<last(dose)~dose-first(dose),
                          first(dose)>last(dose)~dose-last(dose),
                          TRUE~NA_real_))
                    }                  
                )

##need to impute missing SE for estimates of both continuous and categorical exposures
#will first do it for categorical exposures as we need to find SE for the outcome of each level
dflist3<-map(dflist2, function(x){
  replna.lev <- x %>%
    filter(!is.na(se) & exp_var=="categorical")%>% 
    group_by(id) %>%
    nest() %>% 
    mutate(Slist = map(data, ~ if (any(is.na(.x$cases_lev) | is.na(.x$n_lev))){
      diag(.x$se[.x$se != 0]^2, nrow = sum(.x$se != 0))
    } else {
      covar.logrr(cases = .x$cases_lev, n = .x$n_lev, y = .x$b, v = I(.x$se^2), 
                  type = "cc")
    }),
    v.list=map(Slist, ~diag(.x)), ##as they are variances, use sqrt to find SE
    v.df =map_dfr(v.list, ~as.data.frame(t(.))), 
    se.df = map_df(v.df, ~sqrt(.x)))
  
   x <- x %>% group_by(id) %>%
    mutate(se_lev=rep_len(c(0L, colMeans(replna.lev$se.df, na.rm=TRUE)), length.out=n()), #repeat vector of mean se
      se= ifelse(is.na(se) & exp_var=="categorical", se_lev, se))
    return(x)
 }
)
##transform categorical estimates after missing data was imputed 
dflist4<-map(dflist3, function(x){
  lev<-x %>% ##running again with missing SE filled
    filter(!is.na(se) & exp_var=="categorical")%>% 
    group_by(id) %>%
    nest() %>% 
    mutate(Slist = map(data, ~ if (any(is.na(.x$cases_lev) | is.na(.x$n_lev))){
      diag(.x$se[.x$se != 0]^2, nrow = sum(.x$se != 0))
    } else {
      covar.logrr(cases = .x$cases_lev, n = .x$n_lev, y = .x$b, v = I(.x$se^2), 
                  type = "cc")
    }),
    lin = map2(data, Slist, 
                  ~ dosresmeta(b ~ dose, se = se, data = .x, covariance = "user",
                               Slist = .y)),
    y=map_dbl(lin, coef), s=sqrt((map_dbl(lin, vcov))))  %>%  slice(n()) %>% select(id, y, s) 
  #y = final estimates for dose response
   #s = final SE for dose response, renamed from b(se)to match metafor style
  
#run meta analysis to get weighted SE used to replace missing SE in outcomes for continuous exposures
  meta.lev <- x %>% filter(exp_var=="categorical") %>% group_by(id) %>% slice(n()) %>% left_join(lev, by=c('id'='id')) 
  meta.cont <- x %>%  filter(exp_var=="continuous") %>% mutate(y=b*(1/dose), s=se*(1/dose))
  replna.meta <- bind_rows(meta.cont, meta.lev)
  replna.main<- replna.meta %>% filter(main==1) 
  
  main <- metagen(TE=y, seTE=s, data = replna.main,
        studlab = author, comb.fixed = FALSE, comb.random = TRUE,
        method.tau = "REML", hakn = FALSE,
         approx.seTE="",  prediction = FALSE,
         sm = "OR", n.e = sample_size,
         id = id)
  se_cont<-main$seTE.random 
  x<-replna.meta %>% mutate(s=ifelse(is.na(se), se_cont, s))
  })

meta_list<-map(dflist4, function (x) {
  x <- x %>% select(-matches('min|max|_re|_lev'), -c(id_dup, dose, dose0, SD, b, se)) %>% mutate(
    yi = y*0.05, si = s*0.05, timelag_cat=ifelse(timelag_n>=6.00, 1, 0)) %>% group_by(covidence_id) %>% 
  select(c(covidence_id:sample_size, id, y:si, everything())) %>%  relocate(timelag_cat, .after=timelag_n
              ) 
  })
#yi and si are 0.05 increases in gini
#id is each estimate, group_id is numeric vector for estimates grouped by article/covidence_id


############create databases now
info_srh <- data.frame(read_excel("cleaning/studyinfo_srh.xlsx", na = "NA"))
info_acm <- data.frame(read_excel("cleaning/studyinfo_acm.xlsx", na = "NA"))
info_list <- list(srh=info_srh, acm=info_acm)


ref_srh <- data.frame(read_excel("cleaning/endnote_srh.xlsx"))
ref_acm <- data.frame(read_excel("cleaning/endnote_acm.xlsx"))
ref_list <- list(srh=data.frame(ref_srh), acm=data.frame(ref_acm))

##clean ref files
ref_list<-map(ref_list, function (x) {
  ##to find if need to put et al., or & after first author for study table
  x <- x %>% mutate(author_1=str_extract(reference, "\\w+(?=\\,)|\\w+-\\w+(?=\\,)"),
                    author_2 = str_extract(reference, "(?<=\\.\\,\\s)\\w+(?=\\,)|(?<=\\s\\&\\s)\\w+"),
                    author_3 = str_extract_all(reference, "(?<=\\&\\s)\\w+", simplify=TRUE),
                    author_3 = ifelse(author_2==author_3[,1], NA_character_, author_3[,1]),
                    pub_year_dup=str_extract(reference, "(?<=\\.\\s)\\d{4}"),
                    author_comb=case_when(is.na(author_2)& is.na(author_3)~paste0(author_1, ", ", pub_year_dup),
                                          !is.na(author_2) & is.na(author_3) ~paste0(author_1, " and ", author_2, ", ", pub_year_dup),
                                          !is.na(author_2) & !is.na(author_3)~paste0(author_1, ", et al., ", pub_year_dup))                    
  ) %>% select(-c(reference:pub_year_dup))
})

complete<-  function(x, y, z) {
  x %>% left_join(y, by = c('covidence_id'='covidence_id')) %>%  left_join(z, by=c('covidence_id'='covidence_id'))}
complete<-purrr::pmap(list(meta_list, ref_list, info_list),  ~complete(x = ..1, y = ..2, z= ..3))


complete<-map(complete, function(x) x %>% relocate(rob, .after = main)%>% relocate(auth_dup, .after = gscale_desc)) 
complete <- map(complete, ~.x %>% mutate(rown = row_number()))
complete$acm <- complete$acm %>% mutate(fu_med=ifelse(is.na(fu_med), map2_dbl(fu_min, fu_max, ~mean(c(.x, .y))), fu_med),
                              fu_cat = case_when(
                                fu_med<5.00~ "0",
                                fu_med>=5.00 & fu_med<=10.00 ~ "1",
                                fu_med > 10.00 ~ "2",
                                TRUE~NA_character_),
                              fu_n = ifelse(!is.na(fu_mn), fu_mn, fu_med)
) %>% relocate(c(fu_n, fu_cat), .after=timelag_n)

complete$srh <- complete$srh %>% mutate(fu_cat=sample(0:2, n(), replace = TRUE)) %>% relocate(fu_cat, .before=timelag_cat)
complete$srh[complete$srh$covidence_id=="#11971", "adjust_ind"] <- "age, gender, income, marital status, race, health care coverage, health check up in last year, smoking habit"
complete$srh[complete$srh$covidence_id=="#11971", "adjust_area"] <- "per capita median income, % mistrust in area"


complete$srh <- complete$srh %>% mutate(rown = row_number())
complete<-map(complete, function (x) {
  x %>% mutate(dplyr::across(ends_with("_cat"), forcats::as_factor)) %>% 
    mutate( timelag_cat = dplyr:::recode_factor(timelag_cat, "0"= "< 6 years", "1"= "> 6 years"),
            usa_cat = dplyr:::recode_factor(usa_cat, "0"="Not-USA", "1"="USA"),
            region_cat = dplyr:::recode_factor(region_cat, "0"="Worldwide", "1"="Asia and Pacific", "2"="Europe and Central Asia", "3"="Latin America and the Caribbean", "4"="North America"),
            gscale_cat=dplyr:::recode_factor(gscale_cat, "0"="local areas (within-country)", "1"="regions (within-country)", "2"="national (between-country)"),
            gini_cat=dplyr:::recode_factor(gini_cat, "0"="<0.30", "1"="0.30-0.40", "2"=">0.40"),
            adjust_cat=dplyr:::recode_factor(adjust_cat, "0"="individual-level variables only", "1"="plus area-level variables"),
            age_cat=dplyr:::recode_factor(age_cat, "0"="< 60 years old", "1"="> 60 years old"),
              rob=forcats::as_factor(rob),
            rob=recode_factor(rob, "0" = "Low", "1" = "Moderate", "2" = "Serious", "3" = "Critical"),
                fu_cat=fct_recode(fu_cat, "<5 years" ="0", "5-10 years" = "1", ">10 years" = "2")
            )
  }
)

complete$srh[complete$srh$covidence_id=="#11971", "adjust_ind"] <- "age, gender, income, marital status, race, health care coverage, health check up in last year, smoking habit"
complete$srh[complete$srh$covidence_id=="#11971", "adjust_area"] <- "per capita median income, % mistrust in area"
complete$srh[complete$srh$covidence_id=="#9933", "study_desc"] <- "National Longitudinal Survey of Youth (NLSY)"
complete$srh[complete$srh$covidence_id=="#2117", "study_desc"] <- "Cancer Prevention Study-II (CPS-II)"

character <-map(complete, function(x) 
  x %>% filter(main==1) %>% select(-c(auth_dup))
)

complete_main<-map(complete, function(x){
  x %>% filter(main==1) %>% select(-c(women_prop:ncol(.)))
}
)

saveRDS(object=complete, file = "cleaning/complete.rds")
saveRDS(object=complete_main, file = "cleaning/complete_main.rds")
saveRDS(object=character, file = "outputs/stchar.rds")

####now all missing SEs are filled, can move onto meta analysis file
##pre-analysis cleaning done

