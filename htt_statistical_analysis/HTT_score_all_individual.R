library(data.table)
library(ggplot2)
library(dplyr)

# --- Step 1: Load main data
DT <- fread("ukbb_htt_all_columns_encoded_11_11_25.tsv.gz",sep='\t') # download from https://drive.google.com/file/d/15-41jACgnZm9yLKjcObkul5nKNfTz_2b/view?usp=drive_link
new_names=colnames(data.frame(fread('ukbb_htt_all_columns_encoded_11_11_25.tsv.gz',sep='\t',header=T,nrows=0),stringsAsFactors=F))
names(DT)=new_names
colnames(DT)=new_names
# --- Step 2: Load scoring info (with Group column)
scoring <- fread("htt_clinical_scoring_system.csv",sep = ',')[Event %in% names(DT)]

# --- Step 3: Select only necessary columns
meta_cols <- names(DT)[1:10] # assuming first 10 are meta columns like p34_Year_of_birth, A2, etc.
event_cols <- scoring$Event
DT <- DT[, c(meta_cols, event_cols), with = FALSE]

# --- Step 4: Melt wide → long format
DT_long <- melt(
  DT,
  id.vars = meta_cols,
  measure.vars = event_cols,
  variable.name = "Event",
  value.name = "Value",
  variable.factor = FALSE
)

# --- Step 5: Keep only non-NA events
DT_long <- DT_long[!is.na(Value)]

# --- Step 6: Join with scoring info to attach Points and Group
DT_long <- DT_long[scoring, on = "Event", nomatch = 0]

# --- Step 7: Within each VCF + Group, keep *maximum* Points
DT_grouped <- DT_long[
  , .(max_points = max(Points, na.rm = TRUE)),
  by = .(VCF, Group)
  ]
DT_grouped$supportive=ifelse(DT_grouped$max_points==1,1,0)
DT_grouped$suggestive2=ifelse(DT_grouped$max_points==3,3,0)
DT_grouped$
  
  =ifelse(DT_grouped$max_points==5,5,0)

# --- Step 8: Compute cumulative (total) score per VCF
DT_scores <- DT_grouped[
  , .(supportive_total_score = sum(supportive, na.rm = TRUE), 
      suggestive2_total_score = sum(suggestive2, na.rm = TRUE), 
      suggestive_total_score = sum(suggestive, na.rm = TRUE)),
  by = VCF
  ]
DT_scores$supportive_total_score[DT_scores$supportive_total_score>3]=3
DT_scores$suggestive2_total_score[DT_scores$suggestive2_total_score>9]=9

# --- Step 9: Merge total_score back to detailed data if desired
DT_long <- DT_long[
  DT_scores,
  on = "VCF"
  ]
DT_long$total_score=DT_long$supportive_total_score+DT_long$suggestive_total_score+DT_long$suggestive2_total_score

# --- Step 10: (optional) Output detailed tally table
data_for_tally <- DT_long[
  , .(
    VCF,
    p34_Year_of_birth,
    A2,
    Event,
    Group,
    Value,
    Points,
    supportive_total_score,
    suggestive_total_score,
    suggestive2_total_score,
    total_score
  )
  ]

# --- Step 11: (optional) Save outputs
# fwrite(DT_scores, "ukbb_total_scores.tsv", sep = "\t")
# fwrite(data_for_tally, "ukbb_data_for_tally.tsv", sep = "\t")

cap_thresh=70
colnames(data_for_tally)=c('VCF','year_of_birth','A2','category','group','year_at_time_of_event','score','supportive','suggestive','suggestive2','cumulative_score')
data_g=unique(data_for_tally[,c('year_of_birth','A2','year_at_time_of_event','VCF','cumulative_score')])
data_g$CAP_at_time_event=(as.numeric(data_g$year_at_time_of_event-data_g$year_of_birth)*(as.numeric(data_g$A2)-30))/6.49
data_g$cap_greater_than90_at_time_of_event=ifelse(data_g$CAP_at_time_event>cap_thresh,'yes','no')
data_g=data_g%>%group_by(VCF)%>%mutate(cap_greater_than90_at_time_of_event=ifelse(any(CAP_at_time_event>cap_thresh),'yes','no'))
g3=ggplot(data_g,aes(y=cumulative_score,x=cap_greater_than90_at_time_of_event))+geom_violin()+xlab(paste0('CAP > ',cap_thresh))
print(g3)
data_for_tally2=data_for_tally
data_for_tally=data_for_tally2
# data_for_tally=data_for_tally[data_for_tally$cumulative_score>=30,]
data_for_tally$CAP_at_time_event=(as.numeric(data_for_tally$year_at_time_of_event-data_for_tally$year_of_birth)*(as.numeric(data_for_tally$A2)-30))/6.49
data_for_tally$cap_greater_than90_at_time_of_event=ifelse(as.numeric(data_for_tally$CAP_at_time_event)>cap_thresh,'yes','no')

data_for_tally_no=data_for_tally[data_for_tally$cap_greater_than90_at_time_of_event=='no',]
total_points=sum(data_for_tally_no$score)
data_for_tally_no=data_for_tally_no%>%
  group_by(category)%>%
  summarise(total=sum(score)/total_points)%>%
  unique()

data_for_tally_yes=data_for_tally[data_for_tally$cap_greater_than90_at_time_of_event=='yes',]
total_points=sum(data_for_tally_yes$score)
data_for_tally_yes=data_for_tally_yes%>%
  group_by(category)%>%
  summarise(total=sum(score)/total_points)%>%
  unique()

comparison=merge(data_for_tally_yes,data_for_tally_no,by='category')
comparison$OR=comparison$total.x/comparison$total.y
colnames(comparison)=c('field','percentage_contribution_to_points_awarded_in_pop_cap_gt70','percentage_contribution_to_points_awarded_in_pop_cap_lt70','OR')
write.table(comparison[order(comparison$OR,decreasing = T),],'compare_point_contributions.csv',sep=',',row.names = F,quote = F)
data_for_tally_no$cap_greater_than90_at_time_of_event='no'
data_for_tally_yes$cap_greater_than90_at_time_of_event='yes'
data_for_tally=rbind(data_for_tally_no,data_for_tally_yes)



data_for_tally$category=factor(data_for_tally$category)
mycolors=cbind(levels(data_for_tally$category),c(RColorBrewer::brewer.pal(12, "Paired")))
data_for_tally$category_label=as.character(data_for_tally$category)
data_for_tally$category_label[data_for_tally$total<0.01]=NA
g2=ggplot(data_for_tally,aes(x=cap_greater_than90_at_time_of_event,y=total,fill=category,col=category,label=category_label))+
  geom_col()+
  geom_text(size = 5,color='black', position = position_stack(vjust = 0.5))+
  scale_color_manual(breaks = mycolors[,1],
                     values = mycolors[,2])+
  scale_fill_manual(breaks = mycolors[,1],
                    values = mycolors[,2])+
  theme(legend.position = 'none')
pdf('test.pdf',height = 30,width = 20)
print(g2)
dev.off()

d=data_for_tally[data_for_tally$A2>=36]
d$CAP_at_time_event=(as.numeric(d$year_at_time_of_event-d$year_of_birth)*(as.numeric(d$A2)-30))/6.49
d=d%>%group_by(VCF)%>%summarise(VCF=unique(VCF),
                                n_datapoints=n(),
                                n_events=length(unique(na.omit(group))),
                                age=max(year_at_time_of_event-year_of_birth),
                                cumulative_score=unique(cumulative_score),
                                final_CAP_score=max(CAP_at_time_event),
                                allele_size=max(A2),
                                died=ifelse(any(grepl(category,pattern='Death')),1,0),
                                diagnosed_with_huntingtons=ifelse(any(grepl(category,pattern='untington')),'Diagnosed with Huntingtons','Not diagnosed with Huntingtons')
                                )


d$allele_size_from36=d$allele_size-36
write.table(d,'htt_3d_plot.tsv',sep = '\t',row.names = F,quote = F)

library(scatterplot3d)
scatterplot3d(x=d$final_CAP_score,
              y=d$cumulative_score,
              z=d$n_events,
              color=ifelse(d$died,'red','blue'),
              pch = ifelse(d$diagnosed_with_huntingtons,16,17),axis=100)
legend("right", legend=c('Not diagnosed','Diagnosed','Died','Alive'),
       col =  c("blue","red","black","black"), pch = c(15,15,16,17))

