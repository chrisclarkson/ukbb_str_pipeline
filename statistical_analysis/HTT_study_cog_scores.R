#trial on ppl with greater than 40
library(data.table)
library(ggplot2)
setwd('~/Downloads')
data=data.frame(fread('data_500k.tsv.gz',header=T),stringsAsFactors=F)
data$VCF=basename(gsub('.vcf','',data$VCF))
data=data[data$gene=='HTT',]
age_at_recruitment=fread('500K_age_at_recruitment.txt',header=T)
data=merge(data,age_at_recruitment,by.x = 'VCF',by.y = 'eid')

data$CAP_at_recruitment=(data$p21022_Age_at_recruitment*(data$A2-30))/6.49
data$age_when_cap_was_90=round((90*6.49)/(data$A2-30),0)

# data=data[data$CAP_at_recruitment>60,]
data=data[!is.na(data$CAP_at_recruitment),]
data=data[,!(colnames(data)%in%c('p21022_Age_at_recruitment'))]
meta_data=data.frame(fread('HTT_project_all_features_meta_file.tsv',sep='\t',header=T),stringsAsFactors=F)
meta_data=meta_data[,!grepl(pattern = '^pNA_Date',colnames(meta_data))]
idx=(colnames(meta_data)%in%colnames(data))
idx[1]=FALSE
meta_data=meta_data[,!idx]
data=merge(data,meta_data,by.x='VCF',by.y='eid')
data$year_when_cap_was_90=data$p34_Year_of_birth+data$age_when_cap_was_90
data$Diagnosed_Huntingtons=NA
data$Diagnosed_Huntingtons[!is.na(data$p131012_Date_G10_first_reported_huntingtons_disease)]='Diagnosed with Huntingtons'
data$Diagnosed_Huntingtons[is.na(data$p131012_Date_G10_first_reported_huntingtons_disease)]='Not Diagnosed with Huntingtons'

cog_score_cols=colnames(data)[grep(pattern='Instance',colnames(data))]
cog_score_cols=cog_score_cols[-grep(pattern='cancer|Array|Plasma|death|Fracture',cog_score_cols)]


convert_column_to_event=function(data,col,value_to_detect='notna'){
  print(col)
  value_to_detect_old=value_to_detect
  if(value_to_detect=='notna'){
    new_col=ifelse(!is.na(data[,col]),1,NA)
  }else if(grepl(pattern = 'gt',value_to_detect)){
    value_to_detect=strsplit(value_to_detect,':')[[1]][2]
    if(grepl(pattern='%',value_to_detect)){
      vals=as.numeric(na.omit(data[,col]))
      quantiles=quantile(vals,seq(0,1,0.01))
      new_col=ifelse(as.numeric(data[,col])<=quantiles[value_to_detect],1,NA)
    }else{
      new_col=ifelse(as.numeric(data[,col])>=value_to_detect,1,NA)
    }
  }else if(grepl(pattern = 'lt',value_to_detect)){
    value_to_detect=strsplit(value_to_detect,':')[[1]][2]
    if(grepl(pattern='%',value_to_detect)){
      vals=as.numeric(na.omit(data[,col]))
      quantiles=quantile(vals,seq(0,1,0.01))
      new_col=ifelse(as.numeric(data[,col])<=quantiles[value_to_detect],1,NA)
    }else{
      new_col=ifelse(as.numeric(data[,col])>=value_to_detect,1,NA)
    }
    }else{
    new_col=ifelse(data[,col]==value_to_detect,1,NA)
  }
  if(grepl(pattern='_i\\d_',col)){
    if(grepl(pattern='_i0_',col)){
      new_col=unlist(lapply(data$p53_i0[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i0_')[[1]][2]
      name=paste0('Date_first_reported_',value_to_detect_old,'_',name)
    }else if(grepl(pattern='_i1_',col)){
      new_col=unlist(lapply(data$p53_i1[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i1_')[[1]][2]
      name=paste0('Date_first_reported_',value_to_detect_old,'_',name)
    }else if(grepl(pattern='_i2_',col)){
      new_col=unlist(lapply(data$p53_i2[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=strsplit(col,'_i2_')[[1]][2]
      name=paste0('Date_first_reported_',value_to_detect_old,'_',name)
    }else if(grepl(pattern='_i3_',col)){
      new_col=unlist(lapply(data$p53_i3[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=gsub('_._Instance_3',"",strsplit(col,'_i3_')[[1]][2])
      name=paste0('Date_first_reported_',name)
    }else if(grepl(pattern='_i4_',col)){
      new_col=unlist(lapply(data$p53_i3[new_col==1],function(x){strsplit(as.character(x),'-')[[1]][1]}))
      name=gsub('_._Instance_4',"",strsplit(col,'_i4_')[[1]][2])
      name=paste0('Date_first_reported_',name)
    }
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

cog_score_cols_out=c()
for(cog_score_col in cog_score_cols){
  data=convert_column_to_event(data,cog_score_col,'lt:10%')
  cog_score_cols_out=c(cog_score_cols_out,colnames(data)[ncol(data)])
}
data_cogs=data[,c('VCF',cog_score_cols_out)]
# data_new=data
# data_plot=data.frame()
# reg_cols=colnames(data_new)[424:ncol(data_new)]
# reg_cols=reg_cols[reg_cols!='p20118_Date_first_reported_Home_area_population_density_._urban_or_rural__Instance_0']
# for(col in  reg_cols){
#   print(col)
#   label=strsplit(col,'_Date_first_reported_')[[1]][2]
#   data_plot_sub=cbind(data_new[,c('VCF','year_when_cap_was_90','A2','Diagnosed_Huntingtons','age_when_cap_was_90','p21022_Age_at_recruitment','CAP_at_recruitment')],as.numeric(data_new[,col])-as.numeric(data_new$year_when_cap_was_90))
#   colnames(data_plot_sub)[ncol(data_plot_sub)]='Years_since_CAP_was_90'
#   data_plot_sub=cbind(data_plot_sub,as.numeric(data_plot_sub[,'Years_since_CAP_was_90'])+as.numeric(data_plot_sub[,'age_when_cap_was_90']))
#   colnames(data_plot_sub)[ncol(data_plot_sub)]='age_at_time_of_event'
#   data_plot_sub=cbind(data_plot_sub,(as.numeric(data_plot_sub[,'age_at_time_of_event'])*(as.numeric(data_plot_sub[,'A2'])-30))/6.49)
#   colnames(data_plot_sub)[ncol(data_plot_sub)]='CAP_at_time_of_event'
#   data_plot_sub$metric=data_new[,grep(pattern=paste0('i\\d_',label),colnames(data_new))]
#   data_plot_sub=cbind(data_plot_sub,label)
#   colnames(data_plot_sub)[ncol(data_plot_sub)]='metric_label'
#   data_plot_sub=data_plot_sub[!is.na(data_plot_sub[,'Years_since_CAP_was_90']),]
#   data_plot=rbind(data_plot,data_plot_sub)
# }
# data_plot=data_plot[!is.na(data_plot$metric),]
# data_plot$metric_label=gsub('__Instance_\\d','',data_plot$metric_label)
# data_plot=unique(data_plot)
# 
# pdf('test.pdf',height = 10,width = 10)
# print(ggplot(data_plot[!((data_plot$Diagnosed_Huntingtons=='Diagnosed with Huntingtons') & (data_plot$A2<36)),],aes(CAP_at_time_of_event,as.numeric(as.character(metric)),col=Diagnosed_Huntingtons))+geom_point(cex=0.1)+geom_smooth(method='lm')+facet_wrap(~metric_label,scales = 'free_y'))
# dev.off()



