library(data.table)
setwd('~/Downloads')
# template=read.csv('ARI_for_chris.csv',header=T)[,1:7]

# data=data.frame(fread('HTT_all_columns_death_dates_included_mental_health_28_9_25.tsv.gz',header=T),stringsAsFactors=F)
# data=data[,!(colnames(data)%in%c('pNA_Date_G10_first_reported_Huntington_disease','p20118_Date_first_reported_lt.10._Home_area_population_density_._urban_or_rural__Instance_0'))]

calculate_rates=function(data,cols_oi,metric_stratifier_df=NULL,cap_stratifier_df2=NULL){
  data_out=data[,c('VCF',cols_oi)]
  from_col=2
  percentages=data.frame()
  if(is.null(metric_stratifier_df)){
    for(name in colnames(data_out)[from_col:ncol(data_out)]){
      portion=sum(!is.na(data_out[,name]))
      denominator=nrow(data_out)
      percentages=rbind(percentages,c(name,portion/denominator,portion,denominator))
    }
    colnames(percentages)=c('measurement','rate','n','n_total')
    percentages$rate=as.numeric(percentages$rate)*100
  }else{
    label_strats=function(x){
      metric_stratifier_df$distance=as.numeric(x)-as.numeric(metric_stratifier_df[,1])
      if(all(metric_stratifier_df$distance<0)){
        label=metric_stratifier_df[1,2]
      }else{
        metric_stratifier_df=metric_stratifier_df[metric_stratifier_df$distance>=0,]
        label=metric_stratifier_df[which.min(metric_stratifier_df$distance),2]
      }
      return(label)
    }
    for(name in colnames(data_out)[from_col:ncol(data_out)]){
      data_out$strat=NA
      vals=data_out[!is.na(data_out[,name]),name]
      if(any(!is.na(vals))){
        if(length(vals)>1){
          replacement=unlist(lapply(vals,label_strats))
        }else{
          replacement=label_strats(vals)
        }
        # print(length(replacement))
        # print(sum(!is.na(data_out[,name])))
        # print(!is.na(data_out[,name]))
        data_out$strat[which(!is.na(data_out[,name]))]=replacement
        for(strat in na.omit(unique(data_out$strat))){
          # print(strat)
          data_out_sub=data_out[data_out$strat==strat,]
          if(nrow(data_out_sub)>0){
            # print(head(data_out_sub[,c('VCF','strat',name)]))
            portion=sum(!is.na(data_out_sub[,name]))
            if(is.null(cap_stratifier_df2)){denominator=nrow(data_out_sub)}else{denominator=cap_stratifier_df2[cap_stratifier_df2[,2]==strat,3]}
            percentages=rbind(percentages,c(name,100*(portion/denominator),portion,denominator,strat))
          }
        }
      }
    }
    colnames(percentages)=c('measurement','rate','n','n_total','category')
  }
  percentages$rate=as.numeric(percentages$rate)
  return(percentages)
}


get_or_stats <- function(n1, total1, n2, total2) {
  # 2x2 table
  mat <- matrix(c(n1, total1 - n1, n2, total2 - n2), nrow = 2, byrow = TRUE)
  
  # Add 0.5 correction for zeros
  if (any(mat == 0)) mat <- mat + 0.5
  if(mat[1,2]<0){print()}
  # Odds ratio and CI
  or <- (mat[1,1] * mat[2,2]) / (mat[1,2] * mat[2,1])
  se_log_or <- sqrt(1/mat[1,1] + 1/mat[1,2] + 1/mat[2,1] + 1/mat[2,2])
  ci_low <- exp(log(or) - 1.96 * se_log_or)
  ci_high <- exp(log(or) + 1.96 * se_log_or)
  
  # Fisher's exact p-value
  pval <- fisher.test(matrix(c(n1, total1 - n1, n2, total2 - n2),
                             nrow = 2, byrow = TRUE))$p.value

  return(c(odds_ratio = or, lower_95 = ci_low, upper_95 = ci_high, p_value = pval))
}


fish_test=function(comparisons){
  # Apply function to each row
  results <- t(apply(comparisons, 1, function(row){
    get_or_stats(as.numeric(row["n.x"]),
                 as.numeric(row["n_total.x"]),
                 as.numeric(row["n.y"]),
                 as.numeric(row["n_total.y"]))
  }))
  
  # Combine with data and apply FDR correction
  df_out <- cbind(comparisons, as.data.frame(results))
  df_out$p_value_fdr <- p.adjust(df_out$p_value, method = "fdr")
  
  # Format for readability
  df_out <- df_out %>%
    mutate(
      odds_ratio = round(odds_ratio, 3),
      lower_95 = round(lower_95, 3),
      upper_95 = round(upper_95, 3),
      p_value = signif(p_value, 3),
      p_value_fdr = signif(p_value_fdr, 3)
    )
  return(df_out[order(as.numeric(df_out$p_value)),])
}

construct_stratifier_df=function(increments,column_oi){
  stratifier_df=data.frame(cap_scores=increments,stages=c(paste(increments[seq(1,length(increments)-1,1)],increments[seq(2,length(increments),1)],sep='-'),'>'))
  stratifier_df=stratifier_df[-nrow(stratifier_df),]
  numbers=c()
  intervals=cbind(increments[seq(1,length(increments)-1,1)],increments[seq(2,length(increments),1)])
  to=nrow(intervals)
  for(i in 1:to){
    numbers=c(numbers,sum((as.numeric(column_oi)>=as.numeric(intervals[i,1])) & (as.numeric(column_oi)<as.numeric(intervals[i,2]))))
  }
  print(stratifier_df)
  stratifier_df$normalisers=numbers
  return(stratifier_df)
}

# quantiles=quantile(data$CAP_at_recruitment,c(seq(0,0.8,0.05),seq(0.82,1,0.01)),na.rm = T)
# quantiles=c(quantiles[quantiles<90])
# cap_stratifier_df=data.frame(cap_scores=unname(quantiles),paste('cap quantile',names(quantiles)))
# colnames(cap_stratifier_df)=c('cap_scores','stages')

# cap_stratifier_df=rbind(cap_stratifier_df,data.frame(cap_scores=seq(90,120,10),stages=c('stage 0','stage 1','stage 2','stage 3')))
# cap_stratifier_df=data.frame(cap_scores=seq(50,120,10),stages=c('<60','60-70','70-80','80-90','stage 0','stage 1','stage 2','stage 3'))




data=data.frame(fread('HTT_all_columns_death_dates_included_mental_health_28_9_25.tsv.gz',header=T),stringsAsFactors=F)
data=data[,!(colnames(data)%in%c('pNA_Date_G10_first_reported_Huntington_disease','p20118_Date_first_reported_lt.10._Home_area_population_density_._urban_or_rural__Instance_0'))]

cols_oi=c(colnames(data)[grep(pattern = '_first_reported_',colnames(data))],colnames(data)[689:758])
cols_oi=colnames(data)[10:ncol(data)]
for(col in cols_oi){
  data_sub=data[!is.na(data[,col]),c('p34_Year_of_birth','A2',col)]
  data_sub$age_at_event=as.numeric(data_sub[,col])-data_sub$p34_Year_of_birth
  if(col=='Date_first_reported__Mood_swings'){
    print(data_sub$age_at_event)
  }
  data_sub[,col]=(data_sub$age_at_event*(data_sub$A2-30))/6.49
  if(col=='Date_first_reported__Mood_swings'){
    print(data_sub[,col])
  }
  data[!is.na(data[,col]),col]=data_sub[,col]
}

increments=seq(0,120,10)
cap_stratifier_df=construct_stratifier_df(increments,data$CAP_at_recruitment)

data_lt27=data[data$A2<27,]
data_gt27=data[data$A2>=27,]
data_gt27=data_gt27[data_gt27$CAP_at_recruitment>=0,]
percentage_cap60=calculate_rates(data_gt27,cols_oi=cols_oi,metric_stratifier_df = cap_stratifier_df,cap_stratifier_df2=cap_stratifier_df )
percentage_lt27=calculate_rates(data_lt27,cols_oi=cols_oi)
comparisons_new=merge(percentage_cap60,percentage_lt27,by='measurement')
comparisons_new$n.x[as.numeric(comparisons_new$n.x)>as.numeric(comparisons_new$n_total.x)]=comparisons_new$n_total.x[as.numeric(comparisons_new$n.x)>as.numeric(comparisons_new$n_total.x)]
tests=fish_test(comparisons_new)
View(tests[(tests$p_value<0.05),])

write.table(tests,'HTT_rate_comparison_CAP_stages_vs_lt27_27_oct.tsv',sep='\t',row.names = F,quote = F)





data=data.frame(fread('HTT_all_columns_death_dates_included_mental_health_28_9_25.tsv.gz',header=T),stringsAsFactors=F)
data=data[,!(colnames(data)%in%c('pNA_Date_G10_first_reported_Huntington_disease','p20118_Date_first_reported_lt.10._Home_area_population_density_._urban_or_rural__Instance_0'))]

cols_oi=c(colnames(data)[grep(pattern = '_first_reported_',colnames(data))],colnames(data)[689:ncol(data)])
for(col in cols_oi){
  data_sub=data[!is.na(data[,col]),c('p34_Year_of_birth','A2',col)]
  # data_sub$age_at_event=data_sub[,col]-data_sub$p34_Year_of_birth
  # data_sub[,col]=(data_sub$age_at_event*(data_sub$A2-30))/6.49
  data[!is.na(data[,col]),col]=data_sub$A2
}

a2_increments=seq(27,50,2)
construct_stratifier_df(a2_increments,data$A2)

data_lt27=data[data$A2<27,]
data_gt27=data[data$A2>=27,]
percentage_gt60=calculate_rates(data_gt27,cols_oi=cols_oi,metric_stratifier_df = stratifier_df,cap_stratifier_df2=stratifier_df)
percentage_lt27=calculate_rates(data_lt27,cols_oi=cols_oi)
comparisons_new=merge(percentage_gt60,percentage_lt27,by='measurement')

tests1=fish_test(comparisons_new)
View(tests1[(tests1$p_value<0.05) & (tests1$n.x>=5),])

write.table(tests1,'HTT_rate_comparison_A2_binned_vs_lt27.tsv',sep='\t',row.names = F,quote = F)



prepare_for_comparison=function(tests){
  scoring_system=data.frame(fread('scoring_system_latest.tsv',header = T),stringsAsFactors = F)
  plot_comparisons=data.frame()
  for(bin in unique(tests$category)){
    tests_sub=tests[tests$category==bin,]
    plot_comparisons=rbind(plot_comparisons,c(bin,sum(tests_sub$measurement[(tests_sub$p_value<0.05) & (tests_sub$odds_ratio>1)]%in%scoring_system$category,na.rm = T),
                                              sum(!(tests_sub$measurement[(tests_sub$p_value<0.05) & (tests_sub$odds_ratio>1)]%in%scoring_system$category),na.rm = T),
                                              'pvalue uncorrected'))
    plot_comparisons=rbind(plot_comparisons,c(bin,sum((tests_sub$measurement[(tests_sub$p_value_fdr<0.05) & (tests_sub$odds_ratio>1)]%in%scoring_system$category),na.rm = T),
                                                  sum(!(tests_sub$measurement[(tests_sub$p_value_fdr<0.05) & (tests_sub$odds_ratio>1)]%in%scoring_system$category),na.rm = T),
                                              'pvalue corrected'))
  }
  colnames(plot_comparisons)=c('bin','n_relevant_pheno','n_irrelevant_pheno','pvalue_handling')
  return(plot_comparisons)
}
plot_comparisons=prepare_for_comparison(tests)
plot_comparisons$bin=factor(plot_comparisons$bin,levels = c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100','100-110','110-120'))
plot_comparisons$metric='CAP score'
plot_comparisons_a2=prepare_for_comparison(tests1)
plot_comparisons_a2$metric='A2'
# plot_comparisons=rbind(plot_comparisons,plot_comparisons_a2)



g=ggplot(plot_comparisons,aes(size=as.numeric(n_relevant_pheno),col=pvalue_handling,x=bin,y=100*(as.numeric(n_relevant_pheno)/108)))+
  geom_point()+
  ylab('Percentage of\nstatistically significant associations\nthat are clinically relevant')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle('Allele size')

g1=ggplot(plot_comparisons_a2,aes(size=as.numeric(n_relevant_pheno),col=pvalue_handling,x=bin,y=100*(as.numeric(n_relevant_pheno)/108)))+
  geom_point()+
  ylab('Percentage of\nstatistically significant associations\nthat are clinically relevant')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ggtitle('Allele size')

gall=ggpubr::ggarrange(g1,g,nrow=2)
pdf('new2.pdf')
print(gall)
dev.off()

