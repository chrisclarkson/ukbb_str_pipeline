#trial on ppl with greater than 40
library(data.table)
library(tidyr)
setwd('~/Downloads')
data=data.frame(fread('data_500k.tsv.gz',header=T),stringsAsFactors=F)
data$VCF=basename(gsub('.vcf','',data$VCF))
data=data[data$gene=='HTT',]
age_at_recruitment=fread('500K_age_at_recruitment.txt',header=T)
data=merge(data,age_at_recruitment,by.x = 'VCF',by.y = 'eid')

data$CAP_at_recruitment=(data$p21022_Age_at_recruitment*(data$A2-30))/6.49
data$age_when_cap_was_90=round((90*6.49)/(data$A2-30),0)

data=data[data$CAP_at_recruitment>90,]
data=data[!is.na(data$CAP_at_recruitment),]
data=data[,!(colnames(data)%in%c('p21022_Age_at_recruitment'))]
meta_data=data.frame(fread('HTT_project_all_features_meta_file.tsv.gz',sep='\t',header=T),stringsAsFactors=F)

idx=(colnames(meta_data)%in%colnames(data))
idx[1]=FALSE
meta_data=meta_data[,!idx]
data=merge(data,meta_data,by.x='VCF',by.y='eid')
data$year_when_cap_was_90=data$p34_Year_of_birth+data$age_when_cap_was_90
data$Diagnosed_Huntingtons=NA
data$Diagnosed_Huntingtons[!is.na(data$p131012_Date_G10_first_reported_huntingtons_disease)]='Diagnosed with Huntingtons'
data$Diagnosed_Huntingtons[is.na(data$p131012_Date_G10_first_reported_huntingtons_disease)]='Not Diagnosed with Huntingtons'

convert_column_to_event=function(data,col,value_to_detect='notna'){
  print(col)
  if(value_to_detect=='notna'){
    new_col=ifelse(!is.na(data[,col]),1,NA)
  }else{
    new_col=ifelse(data[,col]==value_to_detect,1,NA)
  }
  if(grepl(pattern='_i\\d_',col)){
    print('date detected')
    if(grepl(pattern='_i0_',col)){
      new_col=unlist(lapply(data$p53_i0[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i0_')[[1]][2]
      name=paste0('Date_first_reported_',name)
    }else if(grepl(pattern='_i1_',col)){
      new_col=unlist(lapply(data$p53_i1[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i1_')[[1]][2]
      name=paste0('Date_first_reported_',name)
    }else if(grepl(pattern='_i2_',col)){
      new_col=unlist(lapply(data$p53_i2[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i2_')[[1]][2]
      name=paste0('Date_first_reported_',name)
    }else if(grepl(pattern='_i3_',col)){
      new_col=unlist(lapply(data$p53_i3[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=gsub('_._Instance_3',"",strsplit(col,'_i3_')[[1]][2])
      name=paste0('Date_first_reported_',name)
    }
    # field_id=strsplit(col,'_')[[1]][1]
  }else{
    new_col[is.na(new_col)]=0
    # field_id=strsplit(col,'_')[[1]][1]
    name=strsplit(col,'_i\\d_')[[1]][2]
  }
  data=cbind(data,as.numeric(new_col))
  field_id=strsplit(col,'_')[[1]][1]
  new_col_name=paste0(field_id,'_',name)
  colnames(data)[ncol(data)]=new_col_name
  return(data)
}

data=convert_column_to_event(data,'p3005_i0_Fracture_resulting_from_simple_fall_._Instance_0',1)
data=convert_column_to_event(data,'p3005_i1_Fracture_resulting_from_simple_fall_._Instance_1',1)
data=convert_column_to_event(data,'p3005_i2_Fracture_resulting_from_simple_fall_._Instance_2',1)
data=convert_column_to_event(data,'p3005_i3_Fracture_resulting_from_simple_fall_._Instance_3',1)

icd10_translations=fread('icd10_translations.tsv')
translate_icd10_codes=function(data,icd_translations){
  main_diagnosis=data[,'p40001_i0_Underlying_primary_cause_of_death._ICD10_._Instance_0']
  secondary_diagnosis=data[,'p40002_i0_a0_Contributory_secondary_causes_of_death._ICD10_._Instance_0_._Array_0']
  for(code in icd10_translations$icd10){
    description=icd10_translations[icd10_translations$icd10==code,"Description"]
    main_diagnosis=gsub(code,description,main_diagnosis)
    secondary_diagnosis=gsub(code,description,secondary_diagnosis)
  }
  data$main_cause_of_death=main_diagnosis
  data$secondary_cause_of_death=secondary_diagnosis
  return(data)
}
data$death_year=as.numeric(unlist(lapply(data$p40000_i0_Date_of_death_._Instance_0,function(x){return(strsplit(as.character(x),'-')[[1]][1])})))
# data=translate_icd10_codes(data,icd_translations)
data_new=data[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','death_year',"age_when_cap_was_90",'p21022_Age_at_recruitment','p34_Year_of_birth',"CAP_at_recruitment")]
# cols_of_interest=cols_of_interest[grep(pattern='Date',cols_of_interest)]
cols_of_interest=colnames(data)[grep(pattern='Date',colnames(data))]
cols_of_interest=cols_of_interest[-grep(pattern='main_ICD10',cols_of_interest)]
cols_of_interest=cols_of_interest[-grep(pattern='Date_of_death',cols_of_interest)]
# cols_of_interest=cols_of_interest[-grep(pattern='pNA_',cols_of_interest)]
# informally_recorded_names=c('R270','R268','R296','F023','F059','F051','R410','R418','R441','F060','G319',
#                             'F321','R458','F410','R13','J690','R634','S001','S008','S010',
#                             'S018','S0260','S0650','S0660','S099','S3270','S4220','S5201','S521','S5250',
#                             'S6200','S6260','S720','S7210','S7220','S8201','S8220','S8230','S9200','T0240',
#                             'I635','G459','R568','M8190','R69')
# informally_recorded_cols=c()
# for(col in informally_recorded_names){
#   print(col)
#   col=paste0('Date_',col,'_first_reported')
#   col_found=colnames(data)[grep(pattern = col,colnames(data))]
#   print(col_found)
#   informally_recorded_cols=c(informally_recorded_cols,col_found)
# }
# cols_of_interest=c(cols_of_interest,informally_recorded_cols)

for(col in cols_of_interest){
  print(col)
  if(any(!is.na(data[,col]))){
    data_new=cbind(data_new,unlist(lapply(data[,col],function(x){strsplit(as.character(x),'-')[[1]][1]})))
  }else{
    data_new=cbind(data_new,rep(NA,nrow(data_new)))
  }
  colnames(data_new)[ncol(data_new)]=col
}
data_new=data.frame(data_new,stringsAsFactors = F)

data_main=data_new

data_plot=data.frame()
for(col in c(colnames(data_new)[grep(pattern='_Date_',colnames(data_new))],'death_year') ){
  print(col)
  data_plot_sub=cbind(data_new[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','age_when_cap_was_90','p21022_Age_at_recruitment','CAP_at_recruitment')],as.numeric(data_new[,col])-as.numeric(data_new$year_when_cap_was_90))
  colnames(data_plot_sub)[ncol(data_plot_sub)]='Years_since_CAP_was_90'
  if(col=='p3005_Date_Fracture_resulting_from_simple_fall_._Instance_0'){
    print(data_plot_sub[,'Years_since_CAP_was_90'])
  }
  data_plot_sub=cbind(data_plot_sub,as.numeric(data_plot_sub[,'Years_since_CAP_was_90'])+as.numeric(data_plot_sub[,'age_when_cap_was_90']))
  colnames(data_plot_sub)[ncol(data_plot_sub)]='age_at_time_of_event'
  data_plot_sub=cbind(data_plot_sub,(as.numeric(data_plot_sub[,'age_at_time_of_event'])*(as.numeric(data_plot_sub[,'A2'])-30))/6.49)
  colnames(data_plot_sub)[ncol(data_plot_sub)]='CAP_at_time_of_event'
  if(col=='death_year'){
    label=paste0('Death- Primary: ',data_new$main_cause_of_death,'; Secondary: ',data_new$secondary_cause_of_death)
  }else{
    label=strsplit(col,split='first_reported_')[[1]][2]
    
  }
  data_plot_sub=cbind(data_plot_sub,label)
  colnames(data_plot_sub)[ncol(data_plot_sub)]='Diagnosis'
  data_plot_sub=data_plot_sub[!is.na(data_plot_sub[,'Years_since_CAP_was_90']),]
  data_plot=rbind(data_plot,data_plot_sub)
}
data_plot=data_plot[!is.na(data_plot$Diagnosis),]
# data_plot=data_plot[(data_plot[,'Years_since_CAP_was_90']>0) & (data_plot[,'Diagnosed_Huntingtons']!='Diagnosed with Huntingtons'),]
# data_plot=data_plot%>%
#   group_by(VCF)%>%
#   mutate(score=n())



scoring_df=fread('HTT_scoring_system.tsv',header=T)
data_plot <- data_plot %>%
  left_join(scoring_df, by = c("Diagnosis" = "category")) %>%
  mutate(score = replace_na(score, 0)) %>%
  group_by(VCF) %>%
  mutate(total_score = sum(score[CAP_at_time_of_event>=CAP_at_recruitment])) %>%
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
all_scores=unique(out[,c('eid','Diagnosed_Huntingtons_or_not','total')])
write.table(all_scores,'all_scores.csv',sep=',',row.names = F,quote = F)
# 
# out_with_fillers=data.frame()
# for(i in unique(out$eid)){
#   out_sub=out[(out$eid)==i,]
#   out_sub=rbind(data.frame(out_sub),rep(NA,ncol(out_sub)))
#   out_with_fillers=rbind(out_with_fillers,out_sub)
# }
# write.table(out_with_fillers,'500K_icd10_bigger_than_36_cap_gt90.tsv',sep = '\t',row.names = F,quote=F,na="")





