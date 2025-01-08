setwd("C:/Users/echod/Documents/R_figures")

library(ggplot2)
library(dplyr)
library(plotrix)
library(rstatix)


repro <- read.csv("reproductive_span_briggsae_ggplot.csv")
head(repro)
df_repro <- as.data.frame(repro)


#generating the se's
avgdData <- (df_repro 
             %>% group_by(day) 
             %>% get_summary_stats(offspring, type = "mean_se")
) 


#first make a table with true or false, depending on if values are above 0
reproPercent <- repro |> 
  mutate(day = factor(day), Repro = ifelse(offspring > 0, T, F))  |> 
  group_by(day) |>
  dplyr:::summarize(ReproPercent = (sum(Repro)/length(Repro)*100))

#add a column with the percentage reproducing per day
avgdData$percent <- floor(reproPercent$ReproPercent)

# put the ratio reproducing in a graph
ggplot(data=(reproPercent), aes(x=as.factor(day), y = ReproPercent))+
         geom_col()


#make the plot with avdgData
ggplot(avgdData) +
  geom_bar(aes(x=as.factor(day), y=mean), stat="identity")+
  geom_errorbar(aes(x=as.factor(day), y = mean, ymin=mean-se, ymax=mean+se), width=.2)+
  geom_line(aes(x=day, y=percent),stat="identity")+
  geom_text(aes(label=percent, x=day, y=0.95*percent), colour="black")+
  scale_y_continuous(sec.axis = sec_axis(~., name= "Percentage animals reproducing"))+
  xlab("Day")+ 
  ylab("Mean offspring")+
  theme_minimal() +
  theme (
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
  