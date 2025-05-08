

require(ggplot2); require(scales); require(reshape2); 
#install.packages("dplyr")
require(dplyr)
#require(Hmisc)
library("readxl")
library(RColorBrewer)
library("ggsci")
#install.packages("ggrepel")
library("ggrepel")
library(ggpubr)


library(stringr)



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd('/Users/admin/Documents/support')
getwd()


library(RColorBrewer)
my_palette = c(brewer.pal(9, "RdBu")[c(1,2, 3, 7, 9)])
my_palette_lst =as.list(strsplit(my_palette, " "))


dark2_colors <-brewer.pal(n = 8, name = "Dark2")
dark2_colors

set1_colors <-brewer.pal(n = 8, name = "Set1")
set1_colors

accent_colors <-brewer.pal(n = 8, name = "Accent")
accent_colors

#pl_per_contig=read.csv('v11.3_PerClade_v2_wol_EXTENDED_pl_vs_contig_len.txt',sep=" ",header=T)
runtime_df=read.csv('runtime_mem_all.csv',sep=" ",header=T)
#classif_df=read.csv('runtime_mem_classifiers.csv' ,sep=",",header=T)
head(runtime_df)
head(classif_df)

combined_df =read.csv('combined_classifier_dist_time_mem.csv' ,sep=",",header=T)
head(combined_df)



# Total training time
ggplot(aes(x=claded, y=time/60, group = chunked, fill = model),
       data=combined_df[combined_df$cap=='no',])+
  #stat_summary(geom="bar")+
  facet_grid(. ~ chunked, scales="free")+
  geom_bar(position = "stack", stat = "identity", colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3))), group = condition), stat = "summary", fun = "sum", vjust = -0.6, colour = "black") +
  #stat_summary(geom="pointrange", )+
  #coord_trans(y = "log2")+
  scale_fill_brewer(palette = "Dark2", na.translate = F, labels = c( "Classifier", "Distance model"))+
  #scale_y_continuous(breaks = c(0, 5,  10, 20, 30, 40, 50, 60))+
  theme_classic()+
  ylim (NA, 40)+
  xlab("Condition")+
  ylab("Training time (hr)")+
theme(legend.title = element_blank(), legend.position = c(.86,.85),legend.direction = "vertical")

ggplot(aes(x=condition, y=mem, fill = chunked, group = chunked),
       data=runtime_df)+
  #stat_summary(geom="bar", )+
  geom_bar(position = "dodge",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3)))), stat = "summary", fun = "mean", vjust = 1.4, colour = "black" ) +
  #stat_summary(geom="pointrange", )+
  #coord_trans(y = "log2")+
  scale_y_continuous(breaks = c(0, 5,  10, 20, 30, 40, 50, 60))+
  theme_classic()+
  xlab("Condition")+
  ylab("Peak memory (G)")


 # scale_color_manual(name="", values = c(  "#0571b0","#ca0020" ), 
#                     labels = c( "Claded (True)", "Uncladed"))+
#  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
#        legend.position = c(.2,.9),legend.direction = "vertical")



# Capping on classifier and distance models for chunked input

ggplot(aes(x=model, y=time/60, fill = cap),
       data=combined_df[combined_df$paired_w_cap =='yes',])+
  #stat_summary(geom="bar", )+
  geom_bar(aes(group = condition), position = "dodge2",stat = "identity",  colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3))), group = condition), position = position_dodge(width = 0.9), stat = "identity",  vjust = -0.6, colour = "black" ) +
  #stat_summary(geom="pointrange", )+
  #coord_trans(y = "log2")+
  scale_fill_brewer(palette = "Dark2")+
  #scale_y_continuous(breaks = c(0, 5,  10, 20, 30, 40, 50, 60))+
  theme_classic()+
  ylim (NA, 39)+
  xlab("Condition")+
  ylab("Training time (hr)")


ggplot(aes(x=model, y=mem, fill = cap),
       data=combined_df[combined_df$paired_w_cap =='yes',])+
  #stat_summary(geom="bar", )+
  geom_bar(aes(group = condition), position = "dodge2",stat = "identity",  colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3))), group = condition), position = position_dodge(width = 0.9), stat = "identity",  vjust = -0.6, colour = "black" ) +
  #stat_summary(geom="pointrange", )+
  #coord_trans(y = "log2")+
  scale_fill_brewer(palette = "Dark2")+
  #scale_y_continuous(breaks = c(0, 5,  10, 20, 30, 40, 50, 60))+
  theme_classic()+
  ylim (NA, 61)+
  xlab("Condition")+
  ylab("Peak memory (G)")





ggplot(aes(x=condition, y=mem, fill = model),
       data=combined_df[combined_df$cap=='no',])+
  #stat_summary(geom="bar", )+
  geom_bar(position = "stack",stat = "summary",fun = "mean", colour="black", alpha=0.9)+
  geom_text(aes(label = after_stat(sprintf("%.1f", round(y, digits = 3))), group = condition), stat = "summary", fun = "mean", vjust = -0.6, colour = "black" ) +
  #stat_summary(geom="pointrange", )+
  #coord_trans(y = "log2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_y_continuous(breaks = c(0, 5,  10, 20, 30, 40, 50, 60))+
  theme_classic()+
  ylim (NA, 40)+
  xlab("Condition")+
  ylab("Time (hr)")




ggplot(aes(x=reorder(clade_num,size), y=time/60, color = condition,  group = condition),
       #data=per_clade_error)+
       data=runtime_df)+
  stat_summary(geom="line")+
  theme_classic()+
  coord_trans(y = "log2")+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 10, 15, 20, 30, 40, 50))+
  xlab("Clade")+
  ylab("Running time (hrs)")

ggplot(aes(x=condition, y=time/60, fill = condition),
       #data=per_clade_error)+
       data=runtime_df)+
  stat_summary(geom="bar")+
  theme_classic()
  coord_trans(y = "log2")+
  scale_y_continuous(breaks = c(0, 1, 2, 4, 10, 15, 20, 30, 40, 50))+
  xlab("Clade")+
  ylab("Running time (hrs)")


  




  geom_rect(xmin = 0.5,xmax = 1.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 2.5,xmax = 3.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 4.5,xmax = 5.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 6.5,xmax = 7.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 8.5,xmax = 9.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 10.5,xmax = 11.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) +
  geom_rect(xmin = 12.5,xmax = 13.5,
            ymin = -Inf,ymax = Inf,alpha = 1.0, fill="#f7f7f7", color = NA) 
  #scale_color_brewer(palette = my_palette,name="")+
  scale_color_manual(name="", values = c(  "#0571b0","#ca0020" ), 
                     labels = c( "Claded (True)", "Uncladed"))
 
  #geom_violin()+
  theme_classic()+
  stat_summary(geom="line")+
  stat_summary(fun.data=function(x) mean_ci(x,ci=0.99),
               position=position_dodge(width=0.3))+
  coord_cartesian(ylim=c(0,7))+
  theme(axis.text.x = element_text(angle = 0, vjust = 0.7, hjust=0.5),
        legend.position = c(.2,.9),legend.direction = "vertical")

ggsave("pl_per_clade_true_violin_v3.pdf",width=6.2,height = 4)
getwd()


  

##############################################################################
##############################################################################
##############################################################################

  geom_vline(xintercept = 0.75, linetype=3,color="gray50")+
  theme_classic()+
  
  
  scale_x_continuous(name="Support",labels=percent)+
  scale_y_continuous(name="ECDF")+
  scale_color_manual(name="", values = c(my_colors[2], my_colors [1], my_colors[8] ), 
                     labels=c("Consensus", "Main","Bin median"))+
  scale_linetype_manual(name="", values = c(1, 2), 
                        labels=c("Incorrect", "Correct"))+
  theme(legend.position = c(.20,.25), 
        legend.margin=margin(t = 0.0, unit='cm') )+
  guides(colour = guide_legend(title = NULL, order = 1, reverse=TRUE, ), 
         linetype = guide_legend(title = NULL, order = 2, reverse=FALSE,))

  
  
# sed '$!N;s/\n/,/' kmer_pl_3layer_lowerLrRate.out | awk 'BEGIN { FS = "[, ]" } ; {print $9, $21}' > kmer_pl_3layer_lowerLrRate.csv
