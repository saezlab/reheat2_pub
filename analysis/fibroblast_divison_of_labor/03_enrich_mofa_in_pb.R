
# enrich mofa -------------------------------------------------------------

mat <- read.csv("output/pb_state_patient_profiles.csv")%>%
  as.data.frame()%>%
  column_to_rownames("X") %>% 
  as.matrix()

pb_target<- read.csv("output/pb_state_patient_profiles_target.csv")
g.load<- read_csv("output/mofa/gene_loadings.csv")
col_list<- readRDS("color_list_figures.rds")

mat<-mat[,colnames(mat) %in% pb_target$sample_state]

net.mof <- g.load%>%
  filter(ctype=="Fib")%>%
  dplyr::rename(source= Factor, 
                target= "feature", 
                weight= "value")%>%
  filter(abs(weight)>0.1)%>% 
  mutate(weight= weight * -1)


# enrich mofa factors in each patient x cell state
res<-run_ulm(mat, net.mof, .mor ="weight")

# add meta data
res2 <- res%>% left_join(pb_target , by= c("condition"= "sample_state"))%>%
  group_by(batch)%>%
  mutate(score_sc = scale(score))

res2%>% 
  ggplot(., aes(x= cell_state, y= (score), fill = heart_failure ))+
  facet_wrap(~source, scales ="free_y")+
  geom_boxplot()

factors= res$source%>% unique()

# for each factor we will fit a linear mixed model and extract R2 values for 
#cs and hf variables
excluded_feats <- c()
res_lmm= lapply(factors, function(fact){
  print(fact)
  
  mod.df= res2 %>% filter(source == fact)
  
  #check whether there are samples across all 4 batches
  x= table(mod.df$cell_state, mod.df$heart_failure, mod.df$batch)
  
  # length(attr(x, "dimnames")[[3]])
  # 
  # if(length(attr(x, "dimnames")[[3]]) < 4 ){
  #   excluded_feats<<-c(excluded_feats,fact)
  #   print("skipped")
  #   return(NULL)
  # }
  # 
  fit.clust = lmer(score_sc ~ 1 + cell_state +  heart_failure + (1| batch:sample_id) , data = mod.df)
  
  #as a fixed effect: 
  #fit.clust = lm(score ~ 1 + cell_state +  heart_failure + batch, data = mod.df)
  # fit.clust = lmer(value ~ 1 + cell_state + heart_failure + (1| batch)+ (1| sample_id) , data = mod.df)
  # fit.clust = lmer(value ~ 1 + cell_state + heart_failure + (1| batch) , data = mod.df)
  summary(fit.clust)
  #r2beta(fit.clust)
  #extract random effect variance 
  a =summary(fit.clust)$varcor$batch[1] # random effect batch variance 
  b = attr(summary(fit.clust)$varcor, "sc")^2  # random effect residual variance 
  #calc explained var by random effect (batch)
  r.eff.var= a/(a+b)
  
  # use this to get R2 values
  df= r2beta(fit.clust, partial = TRUE)%>% as.data.frame()
  
  df= df %>% mutate(r.eff.var= r.eff.var, 
                    feature = fact)
  
  return(df)
})%>% do.call(rbind, .)
res_lmm

df.plot= res_lmm %>% select(Effect, Rsq, r.eff.var, feature)%>%
  pivot_wider(names_from = Effect, values_from= Rsq)%>%
  mutate(ratio= heart_failure/cell_state, 
         logratio= log10(ratio))%>%
  mutate(cat = ifelse(heart_failure<0.1 & cell_state>0.1, "identity", "low"),
         cat = ifelse(heart_failure>0.1 & cell_state<0.1, "disease", cat),
         cat = ifelse(heart_failure>0.1 & cell_state>0.1, "both", cat))

df.plot%>% ggplot(aes(x= heart_failure, y= cell_state))+
  geom_point(aes(col= cat))+
  geom_label_repel(aes(label=feature))

df.plot$feature<- str_replace_all(df.plot$feature, "Factor", "Fib_MCP")

# For each 
p.lollipop.R2 <-df.plot%>% 
  pivot_longer(cols= c(heart_failure, cell_state), names_to = "covariate", values_to = "R2")%>%
  filter(feature %in% c("Fib_MCP1", "Fib_MCP2"))%>%
  mutate(covariate= ifelse(covariate=="cell_state", "CS", "HF"))%>%
  ggplot(aes(x= feature, y= R2))+
  geom_col(color= "black",  width = 0.6, show.legend = F, fill ="grey")+
  #geom_point(size= 4, show.legend = F)+
  facet_grid(~covariate)+
  theme_cowplot()+
  #scale_fill_manual(values= unname(col_list$mcp_colors))+
  scale_fill_manual(values= c("grey"))+
  theme(axis.text.x = element_text(angle= 90,vjust= 0.5,  hjust= 1),
        strip.background =element_rect(fill="lavender"),
        strip.text = element_text(colour = 'black'))+
  labs(x= "", y= "partial R²")+
  ylim(c(0,1))+
  geom_hline(yintercept = 0)
  scale_y_continuous(c(0, 0.5, 1))

p.lollipop.R2

pdf("output/figures/fibstates_R2_lolli.pdf", 
    width= 2, height= 2.5)
p.lollipop.R2
dev.off()

p.raw.ulms<- res2%>% 
  filter(source %in% c("Factor1","Factor2"))%>%
  filter(heart_failure=="NF")%>%
  ggplot(., aes(x= cell_state, y= score_sc, fill = batch ))+
  facet_wrap(~source, scales ="free_y")+
  geom_boxplot(width= 0.5)+
  theme_cowplot()+
  theme(axis.text.x = element_text(angle= 45, hjust= 1))+
  labs(x= "", y= "ulm score")

p.raw.ulms

p.raw.ulms<- res2%>% 
  filter(source %in% c("Factor1","Factor2"))%>%
  mutate(source = str_replace_all(source, "Factor", "Fib_MCP"))%>%
  ggplot(., aes(x= cell_state, y= score_sc, fill = heart_failure, col= cell_state))+
  facet_wrap(~source)+
  geom_boxplot(width= 0.5, size = .5)+
  scale_fill_manual(values= c( "white", "darkgrey"))+
  theme_cowplot()+
  scale_color_manual(values= c(unname(col_list$fib_state_colors), "black"))+
  theme(axis.text.x = element_text(angle= 90,vjust= 0.5,  hjust= 1),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5, 
                                        linetype = "solid"),
        strip.background =element_rect(fill="lavender"),
        strip.text = element_text(colour = 'black'))+
  labs(x= "", y= "Expression score", fill = "")+
  guides(color = "none" )
p.raw.ulms

pdf("output/figures/fibstates_R2_rawulm.pdf", 
    width= 5.5, height= 2.5)
p.raw.ulms
dev.off()

plot_grid(p.raw.ulms,)

# associate individual composition with over expression  ------------------


prop_data= read.csv("output/fib_sub_analysis/cluster_composition_sub.csv")
meta_data= read.csv("output/fib_sub_analysis/obs_df_sub.csv")
meta_data= meta_data%>% select(sample_id, disease_code, heart_failure, batch, disease_code)%>% distinct()
prop_data

prop_data<- prop_data %>%
  select(-batch, -heart_failure, -sample_id)%>%
  compositions::acomp() %>%
  compositions::clr()%>% 
  bind_cols(prop_data %>%
              select(batch, heart_failure, sample_id))

prop_data= prop_data %>% pivot_longer(cols = c(paste0("X", 0:5)), 
                                      names_to = "cell_type", 
                                      values_to= "cell_type_prop")%>%
  mutate(cell_type = str_replace_all(cell_type, "X", "_"),
         sample_state = paste0(sample_id, cell_type))
prop_data

g.load_filt = g.load %>% 
  filter(ctype=="Fib",
         abs(value) >0.1, 
         Factor %in% c("Factor1", "Factor2"))

tidy_res <- res2%>% left_join(prop_data%>% select(sample_state, cell_type_prop) , by= c("condition"= "sample_state"))

tidy_res %>% 
  filter(source%in% c("Factor1", "Factor2"))%>%
  ggplot(.,aes(x= -score_sc, y= cell_type_prop, color= source))+
  geom_point(aes(shape= heart_failure))+
  facet_wrap(~cell_state + source, scales = "free", ncol = 4)+
  geom_smooth(method="lm", color="black")+
  theme_cowplot()
