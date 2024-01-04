 
##use main for list

setwd("C:/Users/mshim/OneDrive/Documents/r_projects/incominequal_sr")

#######################################################start with some cleaning
analysis <- readRDS("analysis.rds")

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

pred.func <- function (x){
  x %>% enframe %>% unnest_longer(value, values_to = "list", indices_to = "group") %>% rowwise() %>% mutate(
    broom = list(broom::tidy(list, conf.int = TRUE)), se.tau2=list$se.tau2, tau2=list$tau2.f, len=length(list$b))
}
modelreg_pred <- map(modelreg, pred.func)

pred2.func <- function(x){
  x.2<- x %>% filter(len==2) %>% mutate(predls=list(predict(object=list, newmods=rbind(c(0, 1)), trans=exp, digits=3)))
  x_pred.2<-map2_df(x.2$group, x.2$predls, data.frame)
  x.2 %>% unnest(broom) %>% bind_cols(x_pred.2)
}
metaregtbl.2 <- map(modelreg_pred, pred2.func)

pred3.func <- function(x){
  x.3<- x %>% filter(len==3) %>% mutate(predls=list(predict(object=list, newmods= rbind(c(0, 0), c(1, 0), c(0, 1)), trans=exp, digits=3)))
  x_pred.3<-map2_df(x.3$group, x.3$predls, data.frame)
  x.3 %>% unnest(broom) %>% bind_cols(x_pred.3)
}
metaregtbl.3 <- map(modelreg_pred, pred3.func)

modelreg.4<- modelreg_pred$acm %>% filter(len==4) %>% mutate(predls=list(predict(object=list, newmods= rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0), c(0, 0, 1)), trans=exp, digits=3)))
metareg_pred.4<-map2_df(modelreg.4$group, modelreg.4$predls, data.frame)
metaregtbl.4 <- list()
metaregtbl.4$acm <-modelreg.4 %>% unnest(broom) %>% bind_cols(metareg_pred.4)

modelreg.5<-modelreg_pred$srh %>% filter(len==5) %>% mutate(predls=list(predict(object=list, newmods=rbind(c(0,0, 0, 0), c(1,0, 0, 0), c(0,1, 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1)), trans=exp, digits=3)))
metareg_pred.5<-map2_df(modelreg.5$group, modelreg.5$predls, data.frame)
metaregtbl.5 <- list()
metaregtbl.5$srh <-modelreg.5 %>% unnest(broom) %>% bind_cols(metareg_pred.5)


metaregbind.func <- function (x, y) { 
  bind_rows(x, y)
}
metaregtbl.a <- purrr::pmap(list(metaregtbl.2, metaregtbl.3), ~metaregbind.func(x=..1, y=..2))
metaregtbl.b<- list()
metaregtbl.b$srh <- metaregtbl.a$srh %>% bind_rows(metaregtbl.5$srh)
metaregtbl.b$acm <- metaregtbl.a$acm %>% bind_rows(metaregtbl.4$acm)
###refls_acm <- list("rob", "gscale_cat", "gini_cat", "region_cat", "usa_cat", "age_cat", "adjust_cat", "fu_cat", "timelag_cat")
refls<-list()
refls$srh <- as.list(unique(metaregtbl.b$srh$group))
refls$acm <- as.list(unique(metaregtbl.b$acm$group))
  

reftbl <-list ()
reftbl$srh<- refls$srh %>% 
  map(~{
    analysis$srh %>% group_by_at(.x) %>% 
      summarise(group_ref = paste(ref_id, collapse=","),
                count = n_distinct(ref_id)
      ) %>% 
      rename_at(1, ~"subgroup") %>% 
      mutate(biggroup=.x)
  }) %>% bind_rows() %>% select(biggroup, subgroup, group_ref, count)


reftbl$acm <- refls$acm %>% 
      map(~{
    analysis$acm %>% group_by_at(.x) %>% 
      summarise(group_ref = paste(unique(ref_id), collapse=","),
                count = n_distinct(ref_id)
      ) %>% 
      rename_at(1, ~"subgroup") %>% 
      mutate(biggroup=.x)
  }) %>% bind_rows() %>% select(biggroup, subgroup, group_ref, count)

reftbl$acm <- reftbl$acm %>% mutate(group_ref=stri_replace_all_fixed(group_ref, pattern="s13", replacement="m999"))
reftbl <-map(reftbl, ~ .x %>% mutate(group_ref = str_split(group_ref, ','),
                                     group_ref = map_chr(group_ref, ~toString(unique(.x))),
                                      group_ref = str_squish(group_ref)))
############hyphen list of references
a<-c(1:38)
b<-lag(a)
dash <- paste0("s", b, "-s", a)[-1]
comma <- paste0("s", b, ", s", a)[-1]
reftbl$srh <-reftbl$srh %>% mutate(ref_val = str_remove_all(group_ref, "s"), ref_val=(strsplit(ref_val, ", "))
) %>% mutate(ref_val=map(ref_val, ~as.numeric(unlist(.))),
             ref_val=map(ref_val, ~sort(.)),
             new_vec = map(ref_val, ~unname(tapply(., c(0, cumsum(diff(.) > 1)), function(y) {
               paste0(unique(range(y)), collapse = "-")
             }))),
             new_vec=map_chr(new_vec, ~toString(.)),
             new_vec=paste0("s", new_vec),
             new_vec=stri_replace_all_fixed(new_vec, pattern = c(", ", "-"),
                                            replacement = c(", s","-s"), vectorize_all = FALSE),
             new_vec= stri_replace_all_regex(new_vec, dash, comma, vectorize_all=FALSE),
             group_ref=new_vec
                          )

d<-c(1:14)
e<-lag(d)
dash.acm <- paste0("m", e, "-m", d)[-1]
comma.acm <- paste0("m", e, ", m", d)[-1]
reftbl$acm <-reftbl$acm %>% mutate(ref_val = str_remove_all(group_ref, "m"), ref_val=(strsplit(ref_val, ", "))
) %>% mutate(ref_val=map(ref_val, ~as.numeric(unlist(.))),
             ref_val=map(ref_val, ~sort(.)),
             new_vec = map(ref_val, ~unname(tapply(., c(0, cumsum(diff(.) > 1)), function(y) {
               paste0(unique(range(y)), collapse = "-")
             }))),
             new_vec=map_chr(new_vec, ~toString(.)),
             new_vec=paste0("m", new_vec),
             new_vec=stri_replace_all_fixed(new_vec, pattern = c(", ", "-", "m999"),
                                            replacement = c(", m","-m", "s13"), vectorize_all = FALSE),
             new_vec= stri_replace_all_regex(new_vec, dash.acm, comma.acm, vectorize_all=FALSE),
             group_ref=new_vec
)

refmetaregbind.func <- function (x, y){
  x %>% select(group, term, pred, ci.lb:pi.ub, p.value, se.tau2, tau2) %>%  mutate(dummyvar=str_remove(x$term, x$group)
) %>% relocate(dummyvar, .after=group) %>%  bind_cols(y) %>% group_by(group) %>% select(
  c(group, subgroup, group_ref, count, everything())) %>% select(-c(dummyvar, biggroup, term)
  )  
}

metaregtbl.c<- purrr::pmap(list(metaregtbl.b, reftbl), ~refmetaregbind.func(x=..1, y=..2))

group_ord <- c("Geographical scale", "Mean income inequality", "World region", "Study setting", "Average population age", "Covariate adjustments", "Risk of bias", "Time lag", "Follow-up time")

metaregtbl.c <- map(metaregtbl.c, function(x) {x %>% 
     mutate(group=case_when(
      group == "gscale_cat"~"Geographical scale",
      group == "gini_cat"~"Mean income inequality",
      group == "region_cat"~"World region",
      group == "usa_cat"~"Study setting",
      group == "age_cat"~"Average population age",
      group == "adjust_cat"~"Covariate adjustments",
      group == "rob"~"Risk of bias",
      group == "timelag_cat"~"Time lag",
      group == "fu_cat" ~ "Follow-up time"
    )) %>% group_by(group) %>% mutate(group_n=cur_group_id()) %>% ungroup()})


metareg_tbl.d <-map(metaregtbl.c, function(x) {x %>% mutate(group=factor(group, levels=group_ord)
                                ) %>% arrange(group) %>% group_by(group) %>% mutate(base_ref = paste0("Ref. group: ", subgroup[1], ", n&thinsp;=&thinsp;", count[1], "<sup>", group_ref[1], "</sup>"),
                                                                                                                                                     subgroup_ref= paste0(count, "<sup>", group_ref, "</sup>"),
                                                                                                                                                     subgroup_id=paste0(subgroup, "&thinsp;", "(n&thinsp;=&thinsp;", subgroup_ref, ")")) %>% ungroup() %>% mutate(
                                                                                                                                                       across(c('pred', 'pi.lb', 'pi.ub'), ~format(., digits=2, nsmall=2)),
                                                                                                                                                       p.value=format(p.value, digits=2, nsmall=2, scientific=FALSE),
                                                                                                                                                       tau2=format(tau2, digits=2, nsmall=4),
                                                                                                                                                       est=paste0(pred, " (", pi.lb, " - ", pi.ub, ")") ) %>% select(c(group, subgroup, subgroup_id, base_ref, est, p.value, tau2, group_n)) })

metareg_tbl<-metareg_tbl.d$srh %>% left_join(metareg_tbl.d$acm, by="subgroup") %>% group_by(group.x) %>% dplyr::mutate(pos=1:n()
)%>% arrange(pos, group.x)%>% group_by(group.x) %>% mutate(across(ends_with(".y"), ~ifelse(group_n.x==9, lead(., n=1), .))
) %>% arrange(group.x
)%>% mutate(group=paste0("<b>", group.x, "</b>"),
            base_ref=ifelse(pos==1, paste0(group, " ",base_ref.x, "     ", base_ref.y), NA_real_)) %>% fill(base_ref, .direction="down"
            ) %>% slice(-1) %>% ungroup() %>% select(base_ref,  subgroup_id.x, est.x:tau2.x, subgroup_id.y, est.y:tau2.y
            ) %>% mutate(across(ends_with(".y"), ~ifelse(is.na(.), paste0(" "), .)))

metareg_tbl%>% gt(groupname_col = "base_ref") %>% gt::text_transform(
  locations = gt::cells_row_groups(),
  fn = function(x) {
    purrr::map(x, ~ gt::html(.x))}
) %>% fmt_markdown(columns=c("subgroup_id.x", "subgroup_id.y")) %>% tab_style(
  style = cell_text(align = "left", indent = px(20)),
  locations = cells_body(columns=c("subgroup_id.x", "subgroup_id.y"))) %>% cols_label(
    subgroup_id.x=md("**Covariates (n=38)**"),
    subgroup_id.y=md("**Covariates (n=14)**"),
    est.x = html("<b>OR (95% CI)*</b>"),
    est.y = html("<b>RR (95% CI)<sup>\u2020</sup></b>"),
    p.value.y = html("<b>p value<sup>\u2021</sup></b>"),
    p.value.x = html("<b>p value<sup>\u2021</sup></b>"),
    tau2.x = html("<b>&tau;&sup2;<sup>\u00A7</sup></b>"),
    tau2.y = html("<b>&tau;&sup2;<sup>\u00A7</sup></b>"
    )) %>% tab_source_note(source_note = html("*Odds ratio (OR) (95% confidence intervals (CI)) for each covariate of interest reflects change in poor self-rated health (SRH) per 0.05 increase in Gini coefficient. Calculated using random effects models with restricted maximum likelihood estimate. <br>
                                           <sup>\u2020</sup>Risk ratio (RR) (95% CI) for each covariate of interest reflects change in all-cause mortality per 0.05 increase in Gini coefficient.<br>
                                           <sup>\u2021</sup>For testing difference between null hypothesis (coefficient has no effect) compared to coefficient affects estimate.<br>
                                           <sup>\u00A7</sup>Between study variance explained by the covariates (i.e., residual heterogeneity beyond I<sup>2</sup>). <br>
                                           \u00B6Based on World Bank classifications (World Bank, 2021)"))%>%  gtsave(filename = "mrtab1.html")



metareg_tbl <-map(metaregtbl.c, function(x) {x %>% mutate(group=factor(group, levels=group_ord)
) %>% arrange(group) %>% group_by(group) %>% mutate(base_ref = paste0("ref. group: ", subgroup[1], ", n&thinsp;=&thinsp;", count[1], "<sup>", group_ref[1], "</sup>"),
            subgroup_ref= paste0(count, "<sup>", group_ref, "</sup>"),
            subgroup=paste0("<b>", subgroup, "</b>"),
            subgroup=paste0(subgroup, "&thinsp;", "(n&thinsp;=&thinsp;", subgroup_ref, ")")
                ) %>% slice(-1) %>% ungroup() %>% mutate(
        across(c('pred', 'pi.lb', 'pi.ub'), ~format(., digits=2, nsmall=2)),
            p.value=format(p.value, digits=2, nsmall=2, scientific=FALSE),
             tau2=format(tau2, digits=2, nsmall=4),
             est=paste0(pred, " (", pi.lb, " - ", pi.ub, ")")) %>% select(c(group, subgroup, base_ref, est, p.value, tau2))})

mrtabl_srh <- metareg_tbl$srh %>% mutate(group=paste0("<b>", group, "</b>"),
                                      base_ref=paste0(group, " (", base_ref, ")")
                                      ) %>% select(-c(group)) %>% gt(groupname_col = "base_ref") %>% gt::text_transform(
                                        locations = gt::cells_row_groups(),
                                        fn = function(x) {
                                        purrr::map(x, ~ gt::html(.x))}
                          ) %>% fmt_markdown(columns=c("subgroup")
                                          ) %>% cols_label(
                            subgroup=md("**Covariates**"),
                              est = html("<b>OR (95% CI)*<b>"),
                            p.value = html("<b>p value\u2021<b>"),
                            tau2 = html("<b>&tau;&sup2;\u00A7<b>")
                          ) %>%  cols_align(
                            align = "left",
                            columns = everything()
                          ) %>% tab_source_note(source_note = html("*Odds ratio (OR) (95% confidence intervals (CI)) for each covariate of interest reflects change in poor self-rated health (SRH) per 0.05 increase in Gini coefficient. Calculated using random effects models with restricted maximum likelihood estimate. <br>
                                                          \u2020Risk ratio (RR) (95% CI) for each covariate of interest reflects change in all-cause mortality per 0.05 increase in Gini coefficient.<br>
                                                         \u2021For testing difference between null hypothesis (coefficient has no effect) compared to coefficient affects estimate.<br>
                                                          \u00A7Between study variance explained by the covariates (i.e., residual heterogeneity beyond I<sup>2</sup>). <br>
                                                                    \u00B6Based on World Bank classifications (World Bank, 2021)")) %>% tab_style(
                                                                      style = cell_text(align = "left", indent = px(20)),
                                                                      locations = cells_body(columns=subgroup)
                                                                    ) %>% tab_style(
                                                                      style = cell_text(align = "left"),
                                                                      locations = cells_row_groups(groups=everything())) %>%
                                                                        tab_header(
                                                                          title = "Self-rated health (n=38)"
                                                                        )  %>% tab_style(
                                                                          style = cell_text(align = "center"),
                                                                          locations = cells_title(groups=c("title")))
                                                                    
                        mrtabl_srh  %>%  gtsave(filename = "mrtab1_srh.html")

mrtabl_acm <- metareg_tbl$acm %>% select(-c(group)) %>% gt(groupname_col = "base_ref") %>% gt::text_transform(
  locations = gt::cells_row_groups(),
  fn = function(x) {
    purrr::map(x, ~ gt::html(.x))}
) %>% fmt_markdown(columns=c("subgroup"))  %>%  cols_label(
    subgroup=md("**Covariates**"),
    est = html("<b>RR (95% CI)\u2020<b>"),
    p.value = html("<b>p value\u2021<b>"),
    tau2 = html("<b>&tau;&sup2;\u00A7<b>")
  ) %>%  cols_align(
    align = "left",
    columns = everything()
  ) %>% tab_source_note(source_note = html("*Odds ratio (OR) (95% confidence intervals (CI)) for each covariate of interest reflects change in poor self-rated health (SRH) per 0.05 increase in Gini coefficient. Calculated using random effects models with restricted maximum likelihood estimate. <br>
                                  \u2020Risk ratio (RR) (95% CI) for each covariate of interest reflects change in all-cause mortality per 0.05 increase in Gini coefficient.<br>
                                 \u2021For testing difference between null hypothesis (coefficient has no effect) compared to coefficient affects estimate.<br>
                                  \u00A7Between study variance explained by the covariates (i.e., residual heterogeneity beyond I<sup>2</sup>). <br>
                                            \u00B6Based on World Bank classifications (World Bank, 2021)")) %>% tab_style(
                                              style = cell_text(align = "left", indent = px(20)),
                                              locations = cells_body(columns=subgroup)
                                            ) %>% tab_style(
                                              style = cell_text(align = "left"),
                                              locations = cells_row_groups(groups=everything())) %>%
  tab_header(
    title = "All-cause mortality (n=14)"
  ) %>% tab_style(
    style = cell_text(align = "center"),
    locations = cells_title(groups=c("title")))
mrtabl_acm  %>%  gtsave(filename = "mrtab1_acm.html")

###end of metaregression

