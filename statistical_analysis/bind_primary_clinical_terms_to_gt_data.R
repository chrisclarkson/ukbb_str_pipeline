library(data.table)
library(tidyr)
data=data.frame(fread('data_500k.tsv.gz',header=T),stringsAsFactors=F)
data$VCF=basename(gsub('.vcf','',data$VCF))
data=data[data$gene=='HTT',]

meta_data=data.frame(fread('HTT_project_all_features_meta_file.tsv.gz',sep='\t',header=T),stringsAsFactors=F)

idx=(colnames(meta_data)%in%colnames(data))
idx[1]=FALSE
meta_data=meta_data[,!idx]
data=merge(data,meta_data,by.x='VCF',by.y='eid')
data$year_when_cap_was_90=data$p34_Year_of_birth+data$age_when_cap_was_90
data$CAP_at_recruitment=(data$p21022_Age_at_recruitment*(data$A2-30))/6.49
data$age_when_cap_was_90=round((90*6.49)/(data$A2-30),0)
data=data[data$CAP_at_recruitment>90,]
data=data[!is.na(data$CAP_at_recruitment),]

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

data_new=data[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','death_year',"age_when_cap_was_90",'p21022_Age_at_recruitment','p34_Year_of_birth',"CAP_at_recruitment")]
cols_of_interest=colnames(data)[grep(pattern='Date',colnames(data))]
cols_of_interest=cols_of_interest[-grep(pattern='main_ICD10',cols_of_interest)]
cols_of_interest=cols_of_interest[-grep(pattern='Date_of_death',cols_of_interest)]

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
