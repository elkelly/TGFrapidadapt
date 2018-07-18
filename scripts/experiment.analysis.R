## analysis for paper: Targeted gene flow could facilitate rapid adaptation in an endangered marsupial 
## Kelly, E. & Phillips BL.
## edited by Ella Kelly 08/09/2017

# load packages
library(PropCIs)
library(lme4)
library(ggplot2)

# load data
setwd("~/Dropbox/TeamQuoll/Manuscripts/b-Submitted/Ella/Heritable toad-smarts/analysis/data")
toaddata<-read.csv("exp.data.csv")
pedi<-read.csv("pedigree.csv")
toaddata<-merge(toaddata, pedi[,2:3], "arks")# insert mums

# clean up data
exp1<-subset(toaddata, toaddata$trial==1)
exp1$treatment<-droplevels(exp1$treatment, " ") # remove neophobia level
exp1$attack.type<-droplevels(exp1$attack.type, "") # remove "N" (no appraoch) level

####
## analysis for exp 1. (foraging behaviour) ##
####

## acquistion behaviour ##
# convert to binary responce variable
exp1$Binary <- 0
exp1$Binary[exp1$attack.type == "S"] <- 1

# linear mixed model
model1<-glmer(Binary ~ origin * treatment + (1|mum)+(1|arks), data = exp1, family = binomial)
summary(model1)
model2<-glmer(Binary ~ origin + treatment + (1|mum)+(1|arks), data = exp1, family = binomial)
anova(model1,model2)

#factor effects (take out interaction term)
summary(model2)
model3<-glmer(Binary ~ treatment + (1|mum)+(1|arks), data = exp1, family = binomial)
model4<-glmer(Binary ~ origin + (1|mum)+(1|arks), data = exp1, family = binomial)
anova(model3,model2)
anova(model4,model2)

#plot
con.table<- table(exp1$attack.type, exp1$treatment,exp1$origin)
con.ftable<-ftable(con.table)
frame<-as.data.frame(con.ftable)
names(frame)<-c("attack.type","treatment","origin","count")

frame$origin <- factor(frame$origin, levels = c("NT","HYB","QLD"))
labels<-c(NT="Toad-naïve", HYB="Hybrid",  QLD="Toad-exposed")
ggplot(data = frame, aes(x = treatment, y = count, fill = attack.type)) + 
  geom_bar(stat="identity", position = "fill", colour="black") +
  facet_grid(. ~ origin, labeller=labeller(origin = labels)) + 
  scale_fill_manual(values=c("#000000", "#999999", "#FFFFFF"), labels=c("Bite", "Paw", "Sniff"), guide = guide_legend(title = "Attack type"))+
  theme_bw()+
  theme(axis.title = element_text(size=18))+
  theme(axis.title.y=element_text(margin=margin(0,20,0,0)))+
  theme(legend.title = element_text(size=18))+
  theme(legend.text = element_text(size=13))+
  theme(legend.key.size = unit(1.5, 'lines'))+
  theme(axis.text = element_text(size=13))+
  theme(strip.text.x = element_text(size=16, color="black", face="bold"))+
  theme(strip.background = element_rect(colour="black", fill="white",size=.5, linetype="solid"))+
  labs(x = "Prey type", y= "Proportion exhibiting behaviour")

## invesitgating behaviour ##

#log transform
log.spent<-log(exp1$spent.sec+1)
exp1<-cbind(exp1, log.spent)

#linear mixed model
model1<-lmer(log.spent ~ origin * treatment + (1|mum)+(1|arks), data = exp1)
summary(model1)
model2<-lmer(log.spent ~ origin + treatment + (1|mum)+(1|arks), data = exp1)
summary(model2)
anova(model1,model2)
model3<-lmer(log.spent ~ origin + (1|mum)+(1|arks), data = exp1)
anova(model2,model3)
model4<-lmer(log.spent ~ treatment + (1|mum)+(1|arks), data = exp1)
anova(model2,model4)

## plot
exp1$origin <- factor(exp1$origin, levels = c("NT","HYB","QLD"))
ggplot(exp1, aes(x=origin, y=spent.sec, fill=treatment)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#FFFFFF", "#999999"), labels=c("Mouse", "Toad"), guide = guide_legend(title = "Treatment")) +
  theme_classic()  +
  labs(x = "Origin", y= "Time spent investigating (secs)")+
  scale_x_discrete(breaks=c("HYB", "NT", "QLD"),
                   labels=c( "Hybrid", "Toad-naïve","Toad-exposed"))

####
## analysis for exp 2. (toad consumption) ##
####
exp2<-subset(toaddata, toaddata$leg=="uneaten"|toaddata$leg=="eaten") # select data

#create binary responce variable
exp2$Binary <- 0
exp2$Binary[exp2$leg == "eaten"] <- 1

# linear mixed model
m1<-glmer(Binary~origin+TBZ+ (1|mum)+(1|arks),data=exp2,family=binomial)
summary(m1)
m2<-glmer(Binary~TBZ+(1|mum)+(1|arks),data=exp2,family=binomial)
anova(m1,m2)
m3<-glmer(Binary~origin+(1|mum)+(1|arks),data=exp2,family=binomial)
anova(m1,m3)

#plot
pQ<-6/21 # toad exposed
nQ<-21
seQ<-sqrt(pQ*(1-pQ)/nQ)
pN<-26/42 # toad naive
nN<-42
seN<-sqrt(pN*(1-pN)/nN)
pH<-4/13 # Hybrid
nH<-13
seH<-sqrt(pH*(1-pH)/nH)

props<-c(pN, pH, pQ)
se<-c(seN, seH, seQ)

barx<-barplot(props, ylim=c(0,1), ylab="Proportion eating toad leg (±SE)", names.arg=c("Toad-naïve", "Hybrid", "Toad-exposed"))
segments(barx, props - se, barx, props + se, lwd = 1) #se
arrows(barx, props - se, barx, props + se, lwd = 1, angle=90, length = 0.05)
arrows(barx, props + se, barx, props - se, lwd = 1, angle=90, length = 0.05)

####
## litter numbers ##
####

litters<-read.csv("litter.no.csv")
summary(litters)
aggregate(. ~ pop, litters, function(x) c(mean = mean(x), sd = sd(x)))
anova(lm(litter~pop, data=litters))
