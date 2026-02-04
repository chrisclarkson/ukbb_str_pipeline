library(data.table)
data=data.frame(fread('HTT_all_columns_death_dates_included_mental_health.tsv.gz',header=T),stringsAsFactors=F)
data=data[,!(colnames(data)%in%c('pNA_Date_G10_first_reported_Huntington_disease','p20118_Date_first_reported_lt.10._Home_area_population_density_._urban_or_rural__Instance_0'))]

convert_column_to_event=function(data,col,value_to_detect='notna',new_col_name=NA){
  print(value_to_detect)
  new_col_name_optional=new_col_name
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
  if(!is.na(new_col_name_optional)){
    new_col_name=new_col_name_optional
  }
  colnames(data)[ncol(data)]=new_col_name
  return(data)
}

convert_single_column_to_event=function(data,cols,date_cols,value_to_detect,new_col_name){
  data=cbind(data,substr(as.character(data[,date_cols]),1,4))
  colnames(data)[ncol(data)]=new_col_name
  data[which(data[,cols]!=value_to_detect),new_col_name]=NA
  return(data)
}


# cols_oi=c()
# for(i in 2:7){
#   meta_data=fread(paste0('most_needed_latest',i,'.csv'),sep=',')
#   cols_oi=c(cols_oi,colnames(meta_data)[2:ncol(meta_data)])
#   rm(list = 'meta_data')
# }

encodings=fread('ukbb_encodings.txt',sep=' ',header=T)
bind=data[,c('VCF','p34_Year_of_birth')]
colnames(bind)=c('eid','yob')
for(i in 2:9){
  meta_data=data.frame(fread(paste0('most_needed_latest',i,'.csv'),sep=','),stringsAsFactors = F)
  if(sum(colnames(meta_data)%in%encodings$column)>0){
    bind=merge(bind,meta_data[,c(1,which(colnames(meta_data)%in%na.omit(c(encodings$column,unique(encodings$date_col)))))],by='eid')
  }
}
meta_data=data.frame(fread(paste0('depression.csv'),sep=','),stringsAsFactors = F)
if(sum(colnames(meta_data)%in%encodings$column)>0){
  bind=merge(bind,meta_data[,c(1,which(colnames(meta_data)%in%na.omit(c(encodings$column,unique(encodings$date_col)))))],by='eid')
}
dated_diagnoses = read.csv("dated_diagnoses.csv.gz")[,c(1,8,9,10,11)]
colnames(dated_diagnoses) = c("eid", "p53_i0", "p53_i1", "p53_i2", "p53_i3")
bind=merge(bind,dated_diagnoses,by='eid')
bind=bind[,-grep(pattern='.y$',colnames(bind))]
colnames(bind)=gsub(".x$","",colnames(bind))

bind=bind[,-2]
new_cols=c()
for(i in 1:nrow(encodings)){
  if(encodings$column[i]%in%colnames(bind)){
  if(is.na(encodings$date_col[i])){
    bind=convert_column_to_event(bind,encodings$column[i],
                                 encodings$value_to_detect[i],
                                 encodings$new_col[i])
  }else{
    bind=convert_single_column_to_event(bind,encodings$column[i],
                                        encodings$date_col[i],
                                        encodings$value_to_detect[i],
                                        encodings$new_col[i])
  }
  new_cols=c(new_cols,colnames(bind)[ncol(bind)])
  }else{print(paste(encodings$column[i],'not available'))}
}



data=merge(data,bind[,c('eid',new_cols)],by.x='VCF',by.y='eid')
write.table(data,'ukbb_htt_all_columns_encoded_11_11_25.tsv',sep = '\t',row.names = F,quote = F)

