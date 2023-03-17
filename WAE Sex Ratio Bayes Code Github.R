#Sex ratio binomial models in brms

library(brms)
library(Matrix)
library(tidyverse)
library(tidybayes)
library(ggpp)

############### NO Substock ##########################################
wae<-read.csv("sex_ratio.csv",header=T)

wae %>% filter(Alive=="", Length>399)%>% 
  group_by(Month, Zone) %>% summarize(All=length(Spp))
wae %>% filter(Alive=="", Sex2=="F", Length>399)%>% 
  group_by(Month, Zone) %>% summarize(All=length(Spp))


wae %>% group_by(Year, Month, Sex2) %>% summarize(All=length(Sex2))
wae %>% filter(Alive=="") %>% group_by(Year, Month, Sex2) %>% summarize(Dead=length(Alive))
wae %>% group_by(Year, Month, Sex2) %>% summarize(All=length(Sex2))

d1<-wae %>% filter(Alive=="") %>% group_by(Year, Month, Zone, Sex2) %>% summarize(Female=length(Sex2)) %>%
  filter(Sex2=="F")

d2<-wae %>% filter(Alive=="") %>% group_by(Year, Month, Zone) %>% summarize(total=length(Sex2))

d<-merge(d1, d2)
d<-d %>% select(-Sex2) %>% mutate(Season=plyr::mapvalues(x = Month, from=c(5,9), to=c("Spring", "Fall")),
                                  ZoneWords=plyr::mapvalues(x=Zone, from=c(1,2,3), to=c("Upper","Middle","Lower")),
                                  Zone=as.factor(Zone),
                                  Year=as.factor(Year))
d
str(d)

#total and count of females each in a column
get_prior(Female|trials(total)~Zone*Season +(1|Year),
          data=d,
          family = binomial(link="logit"))

d_binom<-brm(Female|trials(total)~Zone*Season +(1|Year),
             data=d,
             family = binomial(link="logit"), 
             prior = c(prior(normal(0,.5),class="Intercept"),
                       prior(normal(0,.5), class="b"),
                       prior(normal(-1.5, 0.6), coef="Zone3:SeasonSpring"),
                       prior(exponential(25), class="sd")),
            # sample_prior = "only",
             file="SexRatio.rds",
             chains=4, iter=2000, cores=4)
summary(d_binom)
conditional_effects(d_binom, reformula=NA)
pp_check(d_binom,  type="hist")
plot(d_binom)

newdata<-tibble(Season=rep(c("Spring", "Fall"), 6),
                Year=rep(c("2015","2015","2016", "2016"),3),
                Zone=c(rep("1", 4), rep("2",4), 
                       rep("3", 4)), 
                total=rep(100, 12))

newdata2<-tibble(Season=rep(c("Spring", "Fall"), 3),
                 Year=rep(c("2017"),6),
                 Zone=c(rep("1", 2), rep("2",2), 
                        rep("3", 2)), 
                 total=rep(100, 6))

#####
#sample the posterior
pred2<-add_epred_draws(newdata=newdata, d_binom) %>% 
  ungroup() %>% 
  select(.draw, Season,  Zone, Year, total, .epred) 
pred<-pred2%>% 
  pivot_wider(names_from=c(Season, Zone), values_from = .epred)

#summarize means by zone
mean(pred$Spring_1)
mean(pred$Spring_2)
mean(pred$Spring_3)
mean(pred$Fall_1)
mean(pred$Fall_2)
mean(pred$Fall_3)

#summarize SDs by zone
sd(pred$Spring_1)
sd(pred$Spring_2)
sd(pred$Spring_3)
sd(pred$Fall_1)
sd(pred$Fall_2)
sd(pred$Fall_3)

#probability statements
(8000-(sum(pred$Spring_1<45)+sum(pred$Spring_1>55)))/8000
(8000-(sum(pred$Spring_2<45)+sum(pred$Spring_2>55)))/8000
(8000-(sum(pred$Spring_3<45)+sum(pred$Spring_3>55)))/8000
(8000-(sum(pred$Fall_1<45)+sum(pred$Fall_1>55)))/8000
(8000-(sum(pred$Fall_2<45)+sum(pred$Fall_2>55)))/8000
(8000-(sum(pred$Fall_3<45)+sum(pred$Fall_3>55)))/8000

#Violin plot - Figure 2.
Zones<-c("Riverine", "Transition", "Lacustrine")

pred2 %>% 
  mutate(ReservoirZone = plyr::mapvalues(Zone,from=c(1,2,3), to=c("Riverine", "Transition", "Lacustrine")),
         Season2 = plyr::mapvalues(Season, from=c("Spring","Fall"), to=as.factor(c(1,2)))) %>% 
  ggplot(aes(x=ReservoirZone, y=.epred/100,  fill=Season2))+
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = .45, ymax = .55, 
           alpha = .1)+
  geom_violin(position=position_dodge(width = .7))+
  scale_fill_manual(values=c( "gray10", "gray80"), guide="none", 
                    limits=c("1","2"), breaks=c("1","2"))+
  theme_classic()+
  ylab("Proportion of female Walleye")+
  xlab("Reservoir Zone")+
  scale_x_discrete(limits=Zones)+
  scale_y_continuous(limits=c(0.05,.75), labels=seq(0,0.8,0.1),breaks = seq(0,.80,.10))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.title = element_text(size=16,face="bold"),#Removed from plot with 'guide="none" above'
        legend.text = element_text(size=12,face="bold")) #Removed from plot with 'guide="none" above'



#ggsave("Figure 2. Sex ratio violin plots.jpeg", device='jpeg', dpi=600, height=5, width = 6, units="in")

#########################################################################
##############   only fish >400 mm
#########################################################################
d1_400<-wae %>% filter(Length>400, Alive=="") %>% group_by(Year, Month, Zone, Sex2) %>% summarize(Female=length(Sex2)) %>%
  filter(Sex2=="F")

d2_400<-wae %>% filter(Length>400, Alive=="") %>% group_by(Year, Month, Zone) %>% summarize(total=length(Sex2))

d_400<-merge(d1_400, d2_400)
d_400<-d_400 %>% mutate(Season=plyr::mapvalues(x = Month, from=c(5,9), to=c("Spring", "Fall")),
                                  ZoneWords=plyr::mapvalues(x=Zone, from=c(1,2,3), to=c("Upper","Middle","Lower")),
                                  Zone=as.factor(Zone),
                                  Year=as.factor(Year))
d_400
str(d_400)
d_400 %>% summarise(all=sum(total), females=sum(Female))

#total and count of females each in a column
get_prior(Female|trials(total)~Zone*Season +(1|Year),
          data=d,
          family = binomial(link="logit"))

d_binom400<-brm(Female|trials(total)~Zone*Season +(1|Year),
             data=d_400,
             family = binomial(link="logit"), 
             prior = c(prior(normal(0,.5),class="Intercept"),
                       prior(normal(0,.5), class="b"),
                       prior(normal(-1.5, .6), coef="Zone3:SeasonSpring"),
                       prior(exponential(25), class="sd")),
             # sample_prior = "only",
             file="SexRatio400.rds",
             chains=4, iter=2000, cores=4)
summary(d_binom400)
conditional_effects(d_binom400, reformula=NA)
pp_check(d_binom400,  type="hist")
plot(d_binom400)

newdata<-tibble(Season=rep(c("Spring", "Fall"), 6),
                Year=rep(c("2015","2015","2016", "2016"),3),
                Zone=c(rep("1", 4), rep("2",4), 
                       rep("3", 4)), 
                total=rep(100, 12))
pred2_400<-add_epred_draws(newdata=newdata, d_binom400) %>% 
  ungroup() %>% 
  select(.draw, Season, Year, Zone, total, .epred) 
pred_400<-pred2_400%>% 
  pivot_wider(names_from=c(Season, Zone), values_from = .epred)

(8000-(sum(pred_400$Spring_1<45)+sum(pred_400$Spring_1>55)))/8000
(8000-(sum(pred_400$Spring_2<45)+sum(pred_400$Spring_2>55)))/8000
(8000-(sum(pred_400$Spring_3<45)+sum(pred_400$Spring_3>55)))/8000
(8000-(sum(pred_400$Fall_1<45)+sum(pred_400$Fall_1>55)))/8000
(8000-(sum(pred_400$Fall_2<45)+sum(pred_400$Fall_2>55)))/8000
(8000-(sum(pred_400$Fall_3<45)+sum(pred_400$Fall_3>55)))/8000

mean(pred_400$Spring_1)
mean(pred_400$Spring_2)
mean(pred_400$Spring_3)
mean(pred_400$Fall_1)
mean(pred_400$Fall_2)
mean(pred_400$Fall_3)

###Prior Sensitivity
d_binom_sens<-brm(Female|trials(total)~Zone*Season +(1|Year),
                  data=d,
                  family = binomial(link="logit"), 
                  prior = c(prior(normal(0,5),class="Intercept"),
                            prior(normal(0,5), class="b"),
                            prior(normal(0,5), coef="Zone3:SeasonSpring"),
                            prior(exponential(10), class="sd")),
                  # sample_prior = "only",
                  file="SexRatio_Sensitivity2.rds",
                  chains=4, iter=2000, cores=4)
summary(d_binom_sens)
conditional_effects(d_binom_sens, reformula=NA)
pp_check(d_binom_sens,  type="hist")
plot(d_binom_sens)

newdata<-tibble(Season=rep(c("Spring", "Fall"), 6),
                Year=rep(c("2015","2015","2016", "2016"),3),
                Zone=c(rep("1", 4), rep("2",4), 
                       rep("3", 4)), 
                total=rep(100, 12))
pred2_sens<-add_epred_draws(newdata=newdata, d_binom_sens) %>% 
  ungroup() %>% 
  select(.draw, Season, Year, Zone, total, .epred) 
pred_sens<-pred2_sens%>% 
  pivot_wider(names_from=c(Season, Zone), values_from = .epred)

sum(pred_sens$Spring_1>50)/8000
sum(pred_sens$Spring_2>50)/8000
sum(pred_sens$Spring_3>50)/8000
sum(pred_sens$Fall_1>50)/8000
sum(pred_sens$Fall_2>50)/8000
sum(pred_sens$Fall_3>50)/8000

(8000-(sum(pred_sens$Spring_1<45)+sum(pred_sens$Spring_1>55)))/8000
(8000-(sum(pred_sens$Spring_2<45)+sum(pred_sens$Spring_2>55)))/8000
(8000-(sum(pred_sens$Spring_3<45)+sum(pred_sens$Spring_3>55)))/8000
(8000-(sum(pred_sens$Fall_1<45)+sum(pred_sens$Fall_1>55)))/8000
(8000-(sum(pred_sens$Fall_2<45)+sum(pred_sens$Fall_2>55)))/8000
(8000-(sum(pred_sens$Fall_3<45)+sum(pred_sens$Fall_3>55)))/8000

mean(pred_sens$Spring_1)
mean(pred_sens$Spring_2)
mean(pred_sens$Spring_3)
mean(pred_sens$Fall_1)
mean(pred_sens$Fall_2)
mean(pred_sens$Fall_3)

sd(pred_sens$Spring_1)
sd(pred_sens$Spring_2)
sd(pred_sens$Spring_3)
sd(pred_sens$Fall_1)
sd(pred_sens$Fall_2)
sd(pred_sens$Fall_3)

#Figure 3 - Sensitivity Analysis Violin plot
Base<-pred2 %>% 
  mutate(ReservoirZone = plyr::mapvalues(Zone,from=c(1,2,3), to=c("Riverine", "Transition", "Lacustrine")),
         Model = "Base")

over400<-pred2_400 %>% 
  mutate(ReservoirZone = plyr::mapvalues(Zone,from=c(1,2,3), to=c("Riverine", "Transition", "Lacustrine")),
         Model = "Over 400")

Sens<-pred2_sens %>% 
  mutate(ReservoirZone = plyr::mapvalues(Zone,from=c(1,2,3), to=c("Riverine", "Transition", "Lacustrine")),
         Model = "Prior")

df1<-full_join(Base, over400)
df2<-full_join(df1, Sens)


Zones<-c("Riverine", "Transition", "Lacustrine")
Models<-c("Base", "Over 400", "Prior")

library(gghalves)
df2 %>% ggplot(aes(x=ReservoirZone, y=.epred/100,  fill=Model))+
  annotate("rect", xmin = 0.4, xmax = 3.6, ymin = .45, ymax = .55, 
           alpha = .1)+
  geom_half_violin(data=df2[which(substr(df2$Season,1,144000)=='Spring'),],
               side="l")+
  geom_half_violin(data=df2[which(substr(df2$Season,1,144000)=='Fall'),],
              side="r")+
  #stat_summary(data=df2[which(substr(df2$Season,1,144000)=='Fall'),],fun=mean, geom="point", shape=rep(c(4,4,4),3), col=rep(c("White", "Black","Black"),3),show.legend=FALSE, size=rep(c(3,3,3),3), position=position_nudge(x=-0.1))+
  scale_fill_manual(values=c("gray10", "gray45", "gray80"), 
                    labels=Models, 
                    breaks=Models)+
  theme_classic()+
  ylab("Proportion of female Walleye")+
  xlab("Reservoir Zone")+
  labs(fill = "Model")+
  scale_x_discrete(limits=Zones)+
  scale_y_continuous(limits=c(0.0,.8), labels=seq(0,0.8,0.1),breaks = seq(0,.80,.10))+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"),
        axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)),
        axis.title.x = element_text(margin = margin(t = 5, r = 0, b = 0, l = 0)),
        legend.title = element_text(size=16,face="bold"),#Removed from plot with 'guide="none" above'
        legend.text = element_text(size=12,face="bold")) #Removed from plot with 'guide="none" above'

#ggsave("Figure 3. Sensitivity analysis violin plot.jpeg", device='jpeg', dpi=600, height=5, width=6, units="in")

