library(data.table)
library(dplyr)
library(tidyr)
data=data.frame(fread('HTT_all_columns_death_dates_included.tsv.gz',header=T),stringsAsFactors=F)
data=data[,!(colnames(data)%in%c('pNA_Date_G10_first_reported_Huntington_disease','p20118_Date_first_reported_lt.10._Home_area_population_density_._urban_or_rural__Instance_0'))]
data_new=data[data$A2>=36,]

data_plot=data.frame()
for(col in c(colnames(data_new)[10:ncol(data_new)]) ){
  print(col)
  
  data_plot_sub=cbind(data_new[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','age_when_cap_was_90','p21022_Age_at_recruitment','CAP_at_recruitment')],as.numeric(data_new[,col])-as.numeric(data_new$year_when_cap_was_90))
  colnames(data_plot_sub)[ncol(data_plot_sub)]='Years_since_CAP_was_90'

  data_plot_sub=cbind(data_plot_sub,as.numeric(data_plot_sub[,'Years_since_CAP_was_90'])+as.numeric(data_plot_sub[,'age_when_cap_was_90']))
  colnames(data_plot_sub)[ncol(data_plot_sub)]='age_at_time_of_event'
  data_plot_sub=cbind(data_plot_sub,(as.numeric(data_plot_sub[,'age_at_time_of_event'])*(as.numeric(data_plot_sub[,'A2'])-30))/6.49)
  colnames(data_plot_sub)[ncol(data_plot_sub)]='CAP_at_time_of_event'
  
  data_plot_sub=cbind(data_plot_sub,col)
  colnames(data_plot_sub)[ncol(data_plot_sub)]='Diagnosis'
  data_plot_sub=data_plot_sub[!is.na(data_plot_sub[,'Years_since_CAP_was_90']),]
  data_plot=rbind(data_plot,data_plot_sub)
}
data_plot=data_plot[!is.na(data_plot$Diagnosis),]
# data_plot=data_plot[(data_plot[,'Years_since_CAP_was_90']>0) & (data_plot[,'Diagnosed_Huntingtons']!='Diagnosed with Huntingtons'),]
# data_plot=data_plot%>%
#   group_by(VCF)%>%
#   mutate(score=n())



scoring_df=fread('scoring_system_latest.tsv',header=T)
data_plot <- data_plot %>%
  left_join(scoring_df, by = c("Diagnosis" = "category")) %>%
  mutate(score = replace_na(score, 0)) %>%
  group_by(VCF) %>%
  mutate(total_score = sum(score[CAP_at_time_of_event>=90])) %>%
  arrange(desc(A2),VCF,Years_since_CAP_was_90)



starters=unique(data_plot[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','age_when_cap_was_90','p21022_Age_at_recruitment','CAP_at_recruitment')])
starters$Years_since_CAP_was_90=0
starters$age_at_time_of_event=starters$age_when_cap_was_90
starters$CAP_at_time_of_event=90
starters$Diagnosis='CAP was 90'
starters$score=NA
starters$total_score=NA
data_plot=rbind(starters,data_plot)
# recruitments=unique(data_plot[,c('VCF','A2','Diagnosed_Huntingtons','age_when_cap_was_90','p21022_Age_at_recruitment','CAP_at_recruitment')])
recruitments=starters
recruitments$Years_since_CAP_was_90=NA
recruitments$age_at_time_of_event=recruitments$p21022_Age_at_recruitment
recruitments$CAP_at_time_of_event=recruitments$CAP_at_recruitment
recruitments$Diagnosis='recruited'
recruitments$score=NA
recruitments$total_score=NA
data_plot=rbind(data_plot,recruitments)

out = data_plot %>%
  arrange(desc(A2)) %>%
  group_by(VCF) %>%
  mutate(
    Years_since_CAP_was_90 = if_else(
      Diagnosis == "recruited",
      p21022_Age_at_recruitment - age_when_cap_was_90,
      Years_since_CAP_was_90
    )
  ) %>%
  arrange(desc(A2), Years_since_CAP_was_90)

out$total_score[(grepl(pattern='Death',out$Diagnosis)) & (out$CAP_at_time_of_event>90) ]=out$total_score[(grepl(pattern='Death',out$Diagnosis)) & (out$CAP_at_time_of_event>90) ]+5


colnames(out)=c('eid','year_when_cap_was_90','repeat_size','Diagnosed_Huntingtons_or_not','age_when_cap_was_90',
                'age_at_recruitment','CAP_at_recruitment','years_since_CAP_was_90','age_at_time_of_event',
                'CAP_at_time_of_event','Event','score','total')


out_with_fillers=data.frame()
for(i in unique(out$eid)){
  out_sub=out[(out$eid)==i,]
  out_sub=rbind(data.frame(out_sub),rep(NA,ncol(out_sub)))
  out_with_fillers=rbind(out_with_fillers,out_sub)
}
write.table(out_with_fillers,'500K_icd10_bigger_than_36_all_events.tsv',sep = '\t',row.names = F,quote=F,na="")
