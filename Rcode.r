#Effects of predator cues and pesticide resistance on the toxicity of a (bio)pesticide mixture
#Delnat V., Janssens L., and Stoks R. (2020). 
#Pest Management Science, 76: 1448-1455. https://doi.org/10.1002/ps.5658
#R code tested on 17/09/2020


######Packages######
#install packages
install.packages("car")     
install.packages("lme4")     
install.packages("lsmeans")     
install.packages("effects")     
install.packages("afex")     

#load packages
library(car)     
library(lme4)
library(lsmeans)
library(effects)     
library(afex)

#Versions of the packages used
R.version #v3.6.1
packageVersion("car") #v3.0.8
packageVersion("lme4") #v1.1.23
packageVersion("lsmeans") #v2.30.0
packageVersion("effects") #v4.1.4
packageVersion("afex") #v0.27.2


######Import dataset######

#Binomial mortality data of Culex quinquefasciatus
dataBIN=read.csv("./Delnat-et-al_CulexMortality_Binomial.csv", sep=",", na.strings=c(""))
dataBIN$Strain=factor(dataBIN$Strain,levels=c("S-Lab","Ace-1R")) 
str(dataBIN)
#CPF = Chlorpyrifos; Cues = predator cues; Bti = Bacillus thuringiensis israelensis


######Culex quinquefasciatus - Mortality - GLZ######

#effect coding due to interaction effects and unbalanced number of replicated vials per treatment combination (also use Anova type III instead of type II)
set_sum_contrasts()

#assumption 
glmMortality=glm(Mortality ~ (Strain+Bti+CPF+Cues)^4, data=dataBIN, na.action=na.omit, family=quasibinomial(link=logit))
summary(glmMortality)
#Dispersion parameter for quasibinomial family taken to be 1.004346 --> OK

#Generalized linear mixed model with a binomial error structure and the logit link
#individuals as the unit of replication, yet took into account that animals from the same vial were not independent by adding Vial as a random factor
glmerMortality=glmer(Mortality ~ (Strain+Bti+CPF+Cues)^4 + (1|Vial), data=dataBIN, na.action=na.omit, family=binomial(link=logit),
                     control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=1e5)))
Anova(glmerMortality, type="III") 

#contrasts with false discovery rate (fdr) correction for pairwise posthoc comparisons
interact1<-pairs(lsmeans(glmerMortality, ~Bti|Strain, adjust="none"))
interact2<-pairs(lsmeans(glmerMortality, ~Strain|Bti, adjust="none"))
test(rbind(interact1,interact2), adjust="fdr")

interact3<-pairs(lsmeans(glmerMortality, ~CPF|Strain, adjust="none"))
interact4<-pairs(lsmeans(glmerMortality, ~Strain|CPF, adjust="none"))
test(rbind(interact3,interact4), adjust="fdr")

interact5<-pairs(lsmeans(glmerMortality, ~Bti*CPF|Strain, adjust="none"))
interact6<-pairs(lsmeans(glmerMortality, ~Strain|Bti*CPF, adjust="none"))
test(rbind(interact5,interact6), adjust="fdr")

interact7<-pairs(lsmeans(glmerMortality, ~Cues|Strain, adjust="none"))
interact8<-pairs(lsmeans(glmerMortality, ~Strain|Cues, adjust="none"))
test(rbind(interact7,interact8), adjust="fdr")

interact9<-pairs(lsmeans(glmerMortality, ~Cues|Bti, adjust="none"))
interact10<-pairs(lsmeans(glmerMortality, ~Bti|Cues, adjust="none"))
test(rbind(interact9,interact10), adjust="fdr")

#effect plots of full model and of significant interactions
plot(effect(mod=glmerMortality, term="Strain*Bti*CPF*Cues"),type = "response")
plot(effect(mod=glmerMortality, term="Bti*Cues"),type = "response")
plot(effect(mod=glmerMortality, term="Bti*CPF"),type = "response")
plot(effect(mod=glmerMortality, term="Strain*Bti"),type = "response")
plot(effect(mod=glmerMortality, term="Strain*Cues"),type = "response")
