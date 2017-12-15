#################
# NICU antenatal steroids analysis
# Citation: Goldstein ND, Kenaley KM, Locke R, Paul DA. The Joint Effects of Antenatal Steroids and Gestational Age on Improved Outcomes in Neonates. Matern Child Health J. 2017 Nov 10.
# 10/23/15 -- Neal Goldstein
#################


### FUNCTIONS ###

library(gmodels) #CrossTable
library(psych) #describe, describeBy
library(boot) #bootstrapping
library(mediation) #mediation analysis

#returns additive measure from an interaction model, see: http://www.ncbi.nlm.nih.gov/pubmed/17726040
bootAdditive = function(data,index,outcome,method)
{
  bootdata = data[index,]
  if (outcome=="Death")
  {
    model = glm(Death ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=bootdata, family=poisson()) 
    
  } else if (outcome=="IVH") {
    model = glm(IVH_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=bootdata, family=poisson()) 
    
  } else if (outcome=="BPD") {
    model = glm(BPD_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=bootdata, family=poisson()) 
    
  }
  
  if (method=="RERI") {
    return(as.numeric(exp(coef(model)["Gestational_age_additive"] + coef(model)["Mom_steroids_additive0"] + coef(model)["Mom_steroids_additive0:Gestational_age_additive"]) - exp(coef(model)["Gestational_age_additive"]) - exp(coef(model)["Mom_steroids_additive0"]) + 1))
  } else if (method=="AP") {
    RERI = exp(coef(model)["Gestational_age_additive"] + coef(model)["Mom_steroids_additive0"] + coef(model)["Mom_steroids_additive0:Gestational_age_additive"]) - exp(coef(model)["Gestational_age_additive"]) - exp(coef(model)["Mom_steroids_additive0"]) + 1
    return(as.numeric(RERI / exp(coef(model)["Gestational_age_additive"] + coef(model)["Mom_steroids_additive0"] + coef(model)["Mom_steroids_additive0:Gestational_age_additive"])))
  } else if (method=="Synergy") {
    return(as.numeric((exp(coef(model)["Gestational_age_additive"] + coef(model)["Mom_steroids_additive0"] + coef(model)["Mom_steroids_additive0:Gestational_age_additive"]) - 1) / ((exp(coef(model)["Gestational_age_additive"]) - 1) + (exp(coef(model)["Mom_steroids_additive0"]) - 1))))
  }
}

#returns a*b from indirect model
bootIndirect = function(data,index,outcome)
{
  bootdata = data[index,]
  if (outcome=="Death")
  {
    #with mediator of steroids
    eq2 = glm(Death ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_death), data=bootdata, family=poisson()) 
    b =  as.numeric(coef(eq2)["Mom_steroids"])
    
    #relation of exposure with mediator
    eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_death), data=bootdata, family=binomial(link="logit")) 
    a = as.numeric(coef(eq3)["Gestational_age"])
    
  } else if (outcome=="IVH") {
    
    #with mediator of steroids
    eq2 = glm(IVH_severe ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_IVH), data=bootdata, family=poisson()) 
    b =  as.numeric(coef(eq2)[3])

    #relation of exposure with mediator
    eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_IVH), data=subset(bootdata, !is.na(bootdata$IVH_severe)), family=binomial(link="logit")) 
    a = as.numeric(coef(eq3)[2])
    
  } else if (outcome=="BPD") {
    
    #with mediator of steroids
    eq2 = glm(BPD_severe ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_BPD), data=bootdata, family=poisson()) 
    b =  as.numeric(coef(eq2)[3])
    
    #relation of exposure with mediator
    eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_BPD), data=subset(bootdata, !is.na(bootdata$BPD_severe)), family=binomial(link="logit")) 
    a = as.numeric(coef(eq3)[2])
    
  }
  
  return(a*b)
}


### READ DATA ###

load("NICU.2016-04-28.RData")


### SUBSET and CREATE COHORT ###

#no birthweight means not admitted to NICU
NICU = NICU[!is.na(NICU$Birthweight), ]

#limit data to 1997 to 2015
NICU = NICU[NICU$Admission_year>=1997 & NICU$Admission_year<=2015, ]

#limit to single admissions only
NICU = NICU[NICU$Admission_n==1, ]

#limit to preterm 24-33 weeks only
NICU = NICU[!is.na(NICU$Gestational_age) & NICU$Gestational_age>23 & NICU$Gestational_age<34, ]

#limit to VLBW only
#NICU = NICUNICU$Birthweight<=1500, ]

#exclude 105 major genetic anomalies which may have altered chance of ACS and outcomes
NICU = NICU[NICU$ID!=110, ]
NICU = NICU[NICU$ID!=473, ]
NICU = NICU[NICU$ID!=673, ]
NICU = NICU[NICU$ID!=931, ]
NICU = NICU[NICU$ID!=1014, ]
NICU = NICU[NICU$ID!=1458, ]
NICU = NICU[NICU$ID!=1468, ]
NICU = NICU[NICU$ID!=1534, ]
NICU = NICU[NICU$ID!=1743, ]
NICU = NICU[NICU$ID!=1982, ]
NICU = NICU[NICU$ID!=2068, ]
NICU = NICU[NICU$ID!=2376, ]
NICU = NICU[NICU$ID!=2422, ]
NICU = NICU[NICU$ID!=2491, ]
NICU = NICU[NICU$ID!=2518, ]
NICU = NICU[NICU$ID!=2756, ]
NICU = NICU[NICU$ID!=28423, ]
NICU = NICU[NICU$ID!=9690, ]
NICU = NICU[NICU$ID!=10017, ]
NICU = NICU[NICU$ID!=10666, ]
NICU = NICU[NICU$ID!=10667, ]
NICU = NICU[NICU$ID!=10793, ]
NICU = NICU[NICU$ID!=10874, ]
NICU = NICU[NICU$ID!=10976, ]
NICU = NICU[NICU$ID!=11026, ]
NICU = NICU[NICU$ID!=11107, ]
NICU = NICU[NICU$ID!=11247, ]
NICU = NICU[NICU$ID!=11895, ]
NICU = NICU[NICU$ID!=11995, ]
NICU = NICU[NICU$ID!=12401, ]
NICU = NICU[NICU$ID!=12461, ]
NICU = NICU[NICU$ID!=12953, ]
NICU = NICU[NICU$ID!=13042, ]
NICU = NICU[NICU$ID!=13406, ]
NICU = NICU[NICU$ID!=13489, ]
NICU = NICU[NICU$ID!=14269, ]
NICU = NICU[NICU$ID!=14323, ]
NICU = NICU[NICU$ID!=14576, ]
NICU = NICU[NICU$ID!=14603, ]
NICU = NICU[NICU$ID!=14835, ]
NICU = NICU[NICU$ID!=15349, ]
NICU = NICU[NICU$ID!=15698, ]
NICU = NICU[NICU$ID!=16126, ]
NICU = NICU[NICU$ID!=16309, ]
NICU = NICU[NICU$ID!=16817, ]
NICU = NICU[NICU$ID!=16931, ]
NICU = NICU[NICU$ID!=16994, ]
NICU = NICU[NICU$ID!=17418, ]
NICU = NICU[NICU$ID!=17642, ]
NICU = NICU[NICU$ID!=17758, ]
NICU = NICU[NICU$ID!=17766, ]
NICU = NICU[NICU$ID!=18092, ]
NICU = NICU[NICU$ID!=18292, ]
NICU = NICU[NICU$ID!=18828, ]
NICU = NICU[NICU$ID!=18889, ]
NICU = NICU[NICU$ID!=18907, ]
NICU = NICU[NICU$ID!=19123, ]
NICU = NICU[NICU$ID!=19392, ]
NICU = NICU[NICU$ID!=19519, ]
NICU = NICU[NICU$ID!=19739, ]
NICU = NICU[NICU$ID!=19966, ]
NICU = NICU[NICU$ID!=20068, ]
NICU = NICU[NICU$ID!=20099, ]
NICU = NICU[NICU$ID!=20191, ]
NICU = NICU[NICU$ID!=20209, ]
NICU = NICU[NICU$ID!=20678, ]
NICU = NICU[NICU$ID!=21121, ]
NICU = NICU[NICU$ID!=21385, ]
NICU = NICU[NICU$ID!=21707, ]
NICU = NICU[NICU$ID!=21712, ]
NICU = NICU[NICU$ID!=21785, ]
NICU = NICU[NICU$ID!=21970, ]
NICU = NICU[NICU$ID!=22569, ]
NICU = NICU[NICU$ID!=22730, ]
NICU = NICU[NICU$ID!=23673, ]
NICU = NICU[NICU$ID!=23699, ]
NICU = NICU[NICU$ID!=23940, ]
NICU = NICU[NICU$ID!=24009, ]
NICU = NICU[NICU$ID!=24125, ]
NICU = NICU[NICU$ID!=24146, ]
NICU = NICU[NICU$ID!=24246, ]
NICU = NICU[NICU$ID!=24257, ]
NICU = NICU[NICU$ID!=24410, ]
NICU = NICU[NICU$ID!=24418, ]
NICU = NICU[NICU$ID!=24939, ]
NICU = NICU[NICU$ID!=24967, ]
NICU = NICU[NICU$ID!=25392, ]
NICU = NICU[NICU$ID!=25630, ]
NICU = NICU[NICU$ID!=25988, ]
NICU = NICU[NICU$ID!=25996, ]
NICU = NICU[NICU$ID!=26320, ]
NICU = NICU[NICU$ID!=26359, ]
NICU = NICU[NICU$ID!=26563, ]
NICU = NICU[NICU$ID!=26733, ]
NICU = NICU[NICU$ID!=26808, ]
NICU = NICU[NICU$ID!=26888, ]
NICU = NICU[NICU$ID!=26893, ]
NICU = NICU[NICU$ID!=26894, ]
NICU = NICU[NICU$ID!=26973, ]
NICU = NICU[NICU$ID!=27004, ]
NICU = NICU[NICU$ID!=27425, ]
NICU = NICU[NICU$ID!=27706, ]
NICU = NICU[NICU$ID!=27758, ]
NICU = NICU[NICU$ID!=28130, ]
NICU = NICU[NICU$ID!=28145, ]


### RECODE ###

#center admission year
NICU$Admission_year_centered = scale(NICU$Admission_year, center=T, scale=F)

#center gestational age
NICU$Gestational_age_centered = scale(NICU$Gestational_age, center=T, scale=F)

#center birthweight
NICU$Birthweight_centered = scale(NICU$Birthweight, center=T, scale=F)

#severe IVH (grade III or IV)
NICU$IVH_severe = ifelse((NICU$IVH==1) & (NICU$IVH_grade>=2), 1, 0)

#moderate or severe BPD
NICU$BPD_severe = ifelse((NICU$BPD==1) & (NICU$BPD_grade>=1), 1, 0)

#death
NICU$Death = ifelse(!is.na(NICU$Type_discharge) & NICU$Type_discharge==1, 1, 0)

#composite outcomes
NICU$Death_IVH = ifelse(NICU$Death==1 | NICU$IVH_severe==1, 1, 0)
NICU$Death_BPD = ifelse(NICU$Death==1 | NICU$BPD_severe==1, 1, 0)

#black or hispanic
NICU$Race_ethnicity = ifelse(NICU$Mom_race_ethnicity==1, 1, ifelse(NICU$Mom_race_ethnicity==2, 2, 0))

#ambiguous sex
NICU$Sex = ifelse(!is.na(NICU$Sex) & NICU$Sex=="A", NA, NICU$Sex)

#cohort time to event: outcome and censoring
#this should be adjusted to timing from ACS receipt, not NICU LOS
NICU$Cohort_time_death = NA
NICU$Cohort_time_IVH = NA
NICU$Cohort_time_BPD = NA
#NICU$Cohort_time_CLD = NA

#sensitivity analysis of exposure
set.seed(777)
NICU$Mom_steroids_sensitivity = NA
for (i in 1:nrow(NICU))
{
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #outcome: death
  NICU$Cohort_time_death[i] = as.numeric(NICU$Date_discharge_initial[i] - NICU$Date_admission[i]) + 1

  #outcome: severe IVH
  if (!is.na(NICU$IVH_severe[i]) && NICU$IVH_severe[i]==1) {
    NICU$Cohort_time_IVH[i] = as.numeric(NICU$IVH_onset_date[i] - NICU$Date_admission[i]) + 1
  } else if (!is.na(NICU$IVH_severe[i])) {
    NICU$Cohort_time_IVH[i] = as.numeric(NICU$Date_discharge_initial[i] - NICU$Date_admission[i]) + 1
  }

  #outcome: severe BPD
  if (!is.na(NICU$BPD_severe[i]) && NICU$BPD_severe[i]==1) {
    NICU$Cohort_time_BPD[i] = as.numeric(NICU$BPD_onset_date[i] - NICU$Date_admission[i]) + 1
  } else {
    NICU$Cohort_time_BPD[i] = as.numeric(NICU$Date_discharge_initial[i] - NICU$Date_admission[i]) + 1
  }

#   #outcome: CLD
#   if (NICU$CLD[i]==1) {
#     NICU$Cohort_time_CLD[i] = as.numeric(NICU$CLD_onset_date[i] - NICU$Date_admission[i]) + 1
#   } else {
#     NICU$Cohort_time_CLD[i] = as.numeric(NICU$Date_discharge_initial[i] - NICU$Date_admission[i]) + 1
#   }
  
  #create a probability of imperfect specificity of steroids
  misclas = ifelse(((NICU$Gestational_age[i]-23)*.01)<=0, 0.005, (NICU$Gestational_age[i]-23)*.01)
  
  #re-assign exposure if probability < random
  NICU$Mom_steroids_sensitivity[i] = ifelse((NICU$Mom_steroids[i]==1) & (runif(1, min=0, max=1)<misclas), 0, NICU$Mom_steroids[i])
}
rm(misclas,i)

#load("Steroids.RData")


### DESCRIPTIVES ###

#overall
CrossTable(NICU$Race_ethnicity)
CrossTable(NICU$Chorio_clinical)
CrossTable(NICU$Mom_HTN)
CrossTable(NICU$PROM)
CrossTable(NICU$Delivery_method)
CrossTable(NICU$Outborn)
CrossTable(NICU$Multiple_gestation)
CrossTable(NICU$Sex)
describe(NICU$Gestational_age)
describe(NICU$Birthweight)
describe(NICU$Admission_year)
describe(NICU$LOS)
CrossTable(NICU$Death)
CrossTable(NICU$IVH_severe)
CrossTable(NICU$BPD_severe)

#by steroids
CrossTable(NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=F)
CrossTable(NICU$Race_ethnicity, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Chorio_clinical, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Mom_HTN, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$PROM, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Delivery_method, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Outborn, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Multiple_gestation, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$Sex, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
describeBy(NICU$Gestational_age, NICU$Mom_steroids); t.test(NICU$Gestational_age ~ NICU$Mom_steroids)
describeBy(NICU$Birthweight, NICU$Mom_steroids); t.test(NICU$Birthweight ~ NICU$Mom_steroids)
describeBy(NICU$Admission_year, NICU$Mom_steroids); t.test(NICU$Admission_year ~ NICU$Mom_steroids)
describeBy(NICU$LOS, NICU$Mom_steroids); t.test(NICU$LOS ~ NICU$Mom_steroids)
CrossTable(NICU$Death, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$IVH_severe, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)
CrossTable(NICU$BPD_severe, NICU$Mom_steroids, prop.r=F, prop.t=F, prop.chisq=F, chisq=T)

#changes in ACS over time
CrossTable(NICU$Mom_steroids[NICU$Admission_year==1997], prop.r=F, prop.t=F, prop.chisq=F, chisq=F)
CrossTable(NICU$Mom_steroids[NICU$Admission_year==2015], prop.r=F, prop.t=F, prop.chisq=F, chisq=F)
summary(glm(Mom_steroids ~ Admission_year, data=NICU, family=binomial(link="logit")))


### CONFOUNDER SELECTION ###

#potential confounders
t.test(NICU$Mom_age ~ NICU$IVH_severe)
chisq.test(NICU$Race_ethnicity, NICU$IVH_severe)
chisq.test(NICU$Delivery_method, NICU$IVH_severe)
chisq.test(NICU$Chorio_clinical, NICU$IVH_severe)
chisq.test(NICU$Mom_HTN, NICU$IVH_severe)
chisq.test(NICU$Outborn, NICU$IVH_severe)
chisq.test(NICU$PROM, NICU$IVH_severe)
chisq.test(NICU$Multiple_gestation, NICU$IVH_severe)
t.test(NICU$Birthweight ~ NICU$IVH_severe)
t.test(NICU$Admission_year ~ NICU$IVH_severe)

t.test(NICU$Mom_age ~ NICU$Mom_steroids)
chisq.test(NICU$Race_ethnicity, NICU$Mom_steroids)
chisq.test(NICU$Delivery_method, NICU$Mom_steroids)
chisq.test(NICU$Chorio_clinical, NICU$Mom_steroids)
chisq.test(NICU$PROM, NICU$Mom_steroids)
chisq.test(NICU$Multiple_gestation, NICU$Mom_steroids)
t.test(NICU$Birthweight ~ NICU$Mom_steroids)
t.test(NICU$Admission_year ~ NICU$Mom_steroids)

summary(lm(Gestational_age ~ as.factor(Race_ethnicity), data=NICU))
summary(lm(Gestational_age ~ Delivery_method, data=NICU))
summary(lm(Gestational_age ~ Chorio_clinical, data=NICU))
summary(lm(Gestational_age ~ PROM, data=NICU))
summary(lm(Gestational_age ~ Multiple_gestation, data=NICU))
summary(lm(Gestational_age ~ Birthweight, data=NICU))
summary(lm(Gestational_age ~ Admission_year, data=NICU))


### INTERACTION ANALYSIS ###

#Interaction_model = glm(Death ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
#Interaction_model = glm(IVH_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
#Interaction_model = glm(BPD_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 
Interaction_model = glm(Death ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
Interaction_model = glm(IVH_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
Interaction_model = glm(BPD_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 

#check for overdispersion, where residual deviance >> degrees of freedom
summary(Interaction_model)

#multiplicative interaction
round(coef(Interaction_model),2)
round(exp(coef(Interaction_model)),2)
round(exp(confint(Interaction_model)),2)

#note for additive assessment, we assign category with the lowest risk as the referent group: see http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3115067/
#and transform gestational age relative to term so increase in number is now a risk factor
NICU$Mom_steroids_additive = relevel(as.factor(NICU$Mom_steroids), ref=2)
NICU$Gestational_age_additive = scale(40 - NICU$Gestational_age, center=T, scale=F)
Interaction_model = glm(Death ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
Interaction_model = glm(IVH_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
Interaction_model = glm(BPD_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 

#calculate additive interaction: http://www.ncbi.nlm.nih.gov/pubmed/17726040
additive_left = exp(coef(Interaction_model)["Gestational_age_additive"] + coef(Interaction_model)["Mom_steroids_additive0"] + coef(Interaction_model)["Mom_steroids_additive0:Gestational_age_additive"] - 1)
additive_right = exp(coef(Interaction_model)["Gestational_age_additive"] - 1) + exp(coef(Interaction_model)["Mom_steroids_additive0"] - 1)

#presence of additive interaction?
additive_left == additive_right #none
additive_left < additive_right #negative
additive_left > additive_right #positive

#RERI
RERI = exp(coef(Interaction_model)["Gestational_age_additive"] + coef(Interaction_model)["Mom_steroids_additive0"] + coef(Interaction_model)["Mom_steroids_additive0:Gestational_age_additive"]) - exp(coef(Interaction_model)["Gestational_age_additive"]) - exp(coef(Interaction_model)["Mom_steroids_additive0"]) + 1

#AP
AP = RERI / exp(coef(Interaction_model)["Gestational_age_additive"] + coef(Interaction_model)["Mom_steroids_additive0"] + coef(Interaction_model)["Mom_steroids_additive0:Gestational_age_additive"])

#Synergy index
S = (exp(coef(Interaction_model)["Gestational_age_additive"] + coef(Interaction_model)["Mom_steroids_additive0"] + coef(Interaction_model)["Mom_steroids_additive0:Gestational_age_additive"]) - 1) / ((exp(coef(Interaction_model)["Gestational_age_additive"]) - 1) + (exp(coef(Interaction_model)["Mom_steroids_additive0"]) - 1))

#boot_ci = boot(NICU, bootAdditive, 1000, outcome="Death", method="Synergy", parallel="multicore", ncpus=4)
#boot.ci(boot_ci, type="norm", index=1)

#bootstrap for CI
set.seed(777)
S = NA
for (i in 1:1000)
{
  bootindex = sample(nrow(NICU), nrow(NICU), replace=T)
  bootdata = NICU[bootindex,]
  #Interaction_model = glm(Death ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=bootdata, family=poisson()) 
  #Interaction_model = glm(IVH_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=bootdata, family=poisson()) 
  Interaction_model = glm(BPD_severe ~ Mom_steroids_additive*Gestational_age_additive + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=bootdata, family=poisson()) 
  S = c(S, (exp(coef(Interaction_model)["Gestational_age_additive"] + coef(Interaction_model)["Mom_steroids_additive0"] + coef(Interaction_model)["Mom_steroids_additive0:Gestational_age_additive"]) - 1) / ((exp(coef(Interaction_model)["Gestational_age_additive"]) - 1) + (exp(coef(Interaction_model)["Mom_steroids_additive0"]) - 1)))
}
S = S[-1]
quantile(S, probs=c(0.025, 0.975))


### STRATIFIED PLOTS ###

#plot receipt of ACS by gestation: Figure 1
ga = 24:33
prop_acs = c(sum(NICU$Mom_steroids[NICU$Gestational_age==24]) / sum(NICU$Gestational_age==24),
             sum(NICU$Mom_steroids[NICU$Gestational_age==25]) / sum(NICU$Gestational_age==25),
             sum(NICU$Mom_steroids[NICU$Gestational_age==26]) / sum(NICU$Gestational_age==26),
             sum(NICU$Mom_steroids[NICU$Gestational_age==27]) / sum(NICU$Gestational_age==27),
             sum(NICU$Mom_steroids[NICU$Gestational_age==28]) / sum(NICU$Gestational_age==28),
             sum(NICU$Mom_steroids[NICU$Gestational_age==29]) / sum(NICU$Gestational_age==29),
             sum(NICU$Mom_steroids[NICU$Gestational_age==30]) / sum(NICU$Gestational_age==30),
             sum(NICU$Mom_steroids[NICU$Gestational_age==31]) / sum(NICU$Gestational_age==31),
             sum(NICU$Mom_steroids[NICU$Gestational_age==32]) / sum(NICU$Gestational_age==32),
             sum(NICU$Mom_steroids[NICU$Gestational_age==33]) / sum(NICU$Gestational_age==33))

mean(prop_acs); min(prop_acs); max(prop_acs)
summary(lm(prop_acs ~ ga))

tiff("Figure1.tif",height=4,width=6,units='in',res=1200)
#plot(x=ga, y=prop_acs, ylab="Proportion received antenatal steroids", xlab="Gestational age in weeks", ylim=c(0.4,0.8), pch=16, xaxt="n")
#axis(1, at=ga,labels=ga)
#text(25.5, 0.75, "P for trend: 0.18", cex=0.8)
barplot(height=prop_acs, names.arg=paste(ga), ylab="Proportion received antenatal steroids", xlab="Gestational age in weeks", ylim=c(0,1))
text(2, 0.8, "P for trend: 0.18", cex=0.8)
dev.off()

#note there are too few outcomes in longer gestations for inference

irr_plots = data.frame("outcome"=NA, "gestation"=NA, "irr"=NA, "lo95"=NA, "hi95"=NA)
NICU$Gestational_age_window = ifelse(NICU$Gestational_age==24, 1, ifelse(NICU$Gestational_age==25, 2, ifelse(NICU$Gestational_age==26 | NICU$Gestational_age==27, 3, ifelse(NICU$Gestational_age>=28, 4, NA))))

#death
Interaction_model = glm(BPD_severe ~ Mom_steroids*as.factor(Gestational_age_window) + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=1, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=2, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age_window)2"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=3, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age_window)3"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=4, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age_window)4"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=5, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age_window)5"])), "lo95"=NA, "hi95"=NA))

#IVH
Interaction_model = glm(IVH_severe ~ Mom_steroids*as.factor(Gestational_age) + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=24, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=25, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)25"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=26, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)26"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=27, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)27"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=28, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)28"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=29, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)29"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=30, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)30"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=31, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)31"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=32, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)32"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=33, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)33"])), "lo95"=NA, "hi95"=NA))

#BPD
Interaction_model = glm(BPD_severe ~ Mom_steroids*as.factor(Gestational_age) + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=24, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=25, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)25"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=26, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)26"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=27, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)27"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=28, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)28"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=29, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)29"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=30, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)30"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=31, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)31"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=32, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)32"])), "lo95"=NA, "hi95"=NA))
irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=33, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"] + coef(Interaction_model)["Mom_steroids:as.factor(Gestational_age)33"])), "lo95"=NA, "hi95"=NA))


#obtain irr for stratified GA for each outcome
for (i in 24:33)
{
  Interaction_model = glm(Death ~ Mom_steroids, offset=log(Cohort_time_death), data=subset(NICU, Gestational_age==i), family=poisson())
  irr_plots = rbind(irr_plots, data.frame("outcome"="Death", "gestation"=i, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",1])), "hi95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",2]))))
                          
  Interaction_model = glm(IVH_severe ~ Mom_steroids, offset=log(Cohort_time_IVH), data=subset(NICU, Gestational_age==i), family=poisson())
  irr_plots = rbind(irr_plots, data.frame("outcome"="IVH", "gestation"=i, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",1])), "hi95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",2]))))

  Interaction_model = glm(BPD_severe ~ Mom_steroids, offset=log(Cohort_time_BPD), data=subset(NICU, Gestational_age==i), family=poisson())
  irr_plots = rbind(irr_plots, data.frame("outcome"="BPD", "gestation"=i, "irr"=as.numeric(exp(coef(Interaction_model)["Mom_steroids"])), "lo95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",1])), "hi95"=as.numeric(exp(confint(Interaction_model)["Mom_steroids",2]))))
}
rm(i, Interaction_model)
irr_plots = irr_plots[-1, ]

plot(irr_plots$gestation[irr_plots$outcome=="Death"], irr_plots$irr[irr_plots$outcome=="Death"])


### SENSITIVITY ANALYSIS ###

CrossTable(NICU$Mom_steroids)
CrossTable(NICU$Mom_steroids_sensitivity)

#steroids
Interaction_model = glm(Death ~ Mom_steroids_sensitivity*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
Interaction_model = glm(IVH_severe ~ Mom_steroids_sensitivity*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
Interaction_model = glm(BPD_severe ~ Mom_steroids_sensitivity*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 

#outcomes
Interaction_model = glm(Death_IVH ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
Interaction_model = glm(Death_BPD ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Mom_HTN + PROM + Delivery_method + Outborn + Multiple_gestation + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 

summary(Interaction_model)

#multiplicative interaction
round(coef(Interaction_model),2)
round(exp(coef(Interaction_model)),2)
round(exp(confint(Interaction_model)),2)


### SYSTEMS of EQUATIONS MEDIATION ANALYSIS ###

##Death

#without mediator of steroids
eq1 = glm(Death ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
c = as.numeric(coef(eq1)["Gestational_age"])

#with mediator of steroids
eq2 = glm(Death ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + Delivery_method + PROM + Birthweight + Admission_year, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
c_prime =  as.numeric(coef(eq2)["Gestational_age"])
b =  as.numeric(coef(eq2)["Mom_steroids"])
bSE = summary(eq2)$coefficients["Mom_steroids",2]

#relation of exposure with mediator
eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, data=subset(NICU, !is.na(NICU$Death) & !is.na(Cohort_time_death)), family=binomial(link="logit")) 
a = as.numeric(coef(eq3)["Gestational_age"])
aSE = summary(eq3)$coefficients["Gestational_age",2]

##IVH

#without mediator of steroids
eq1 = glm(IVH_severe ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
c = as.numeric(coef(eq1)[2])

#with mediator of steroids
eq2 = glm(IVH_severe ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
c_prime =  as.numeric(coef(eq2)[2])
b =  as.numeric(coef(eq2)[3])
bSE = summary(eq2)$coefficients["Mom_steroids",2]

#relation of exposure with mediator
eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, data=subset(NICU, !is.na(NICU$IVH_severe) & !is.na(NICU$Cohort_time_IVH)), family=binomial(link="logit")) 
a = as.numeric(coef(eq3)[2])
aSE = summary(eq3)$coefficients["Gestational_age",2]

##BPD

#without mediator of steroids
eq1 = glm(BPD_severe ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 
c = as.numeric(coef(eq1)[2])

#with mediator of steroids
eq2 = glm(BPD_severe ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 
c_prime =  as.numeric(coef(eq2)[2])
b =  as.numeric(coef(eq2)[3])
bSE = summary(eq2)$coefficients["Mom_steroids",2]

#relation of exposure with mediator
eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, data=subset(NICU, !is.na(NICU$BPD_severe) & !is.na(NICU$Cohort_time_BPD)), family=binomial(link="logit")) 
a = as.numeric(coef(eq3)[2])
aSE = summary(eq3)$coefficients["Gestational_age",2]

# ##CLD
# 
# #without mediator of steroids
# eq1 = glm(CLD ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_CLD), data=NICU, family=poisson()) 
# c = as.numeric(coef(eq1)[2])
# 
# #with mediator of steroids
# eq2 = glm(CLD ~ Gestational_age + Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, offset=log(Cohort_time_CLD), data=NICU, family=poisson()) 
# c_prime =  as.numeric(coef(eq2)[2])
# b =  as.numeric(coef(eq2)[3])
# bSE = summary(eq2)$coefficients["Mom_steroids",2]
# 
# #relation of exposure with mediator
# eq3 = glm(Mom_steroids ~ Gestational_age + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight + Admission_year, data=subset(NICU, !is.na(CLD) & !is.na(Cohort_time_CLD)), family=binomial(link="logit")) 
# a = as.numeric(coef(eq3)[2])
# aSE = summary(eq3)$coefficients["Gestational_age",2]

#summary statistics
summary(eq1)
summary(eq2)
summary(eq3)

#a, direct effect
round(a,2)
round(exp(coef(eq3)),2)
round(exp(confint(eq3)),2)

#b, direct effect
round(b,2)
round(exp(coef(eq2)),2)
round(exp(confint(eq2)),2)

#c', direct effect
round(c_prime,2)
round(exp(coef(eq2)),2)
round(exp(confint(eq2)),2)

#mediated/indirect effect
c - c_prime #continuous
a*b #dichotomous, should approximate c - c_prime
round(a*b,2)
round(exp(a*b),2)

#sobel test for significance
pooledSE = sqrt(((a^2)*(bSE^2))+((b^2)*(aSE^2)))
t = (a*b)/pooledSE

#check against normal distribution
2*pnorm(-abs(t))

#confidence
exp((a*b) - (1.96*pooledSE))
exp((a*b) + (1.96*pooledSE))

#bootstrap for CI
boot_ci = boot(NICU, bootIndirect, 1000, outcome="CLD", parallel="multicore", ncpus=4)
boot.ci(boot_ci, type="norm", index=1)

#c, total effect
round(c,2)
round(exp(coef(eq1)),2)
round(exp(confint(eq1)),2)

#percent mediated (since exposure is protective, |c| > |c_prime| indicates mediation)
((c - c_prime) / c) * 100


### COUNTERFACTUAL MEDIATION ANALYSIS ###

#see: http://imai.princeton.edu/research/files/mediationR2.pdf
#vignette("mediation")

##Death
model_outcome = glm(Death ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_death), data=NICU, family=poisson()) 
model_mediator = lm(Gestational_age_centered ~ Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, data=subset(NICU, !is.na(NICU$Death) & !is.na(Cohort_time_death))) 

##IVH
model_outcome = glm(IVH_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_IVH), data=NICU, family=poisson()) 
model_mediator = lm(Gestational_age_centered ~ Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, data=subset(NICU, !is.na(NICU$IVH_severe) & !is.na(NICU$Cohort_time_IVH))) 

##BPD
model_outcome = glm(BPD_severe ~ Mom_steroids*Gestational_age_centered + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, offset=log(Cohort_time_BPD), data=NICU, family=poisson()) 
model_mediator = lm(Gestational_age_centered ~ Mom_steroids + Sex + as.factor(Race_ethnicity) + Chorio_clinical + PROM + Delivery_method + Birthweight_centered + Admission_year_centered, data=subset(NICU, !is.na(NICU$BPD_severe) & !is.na(NICU$Cohort_time_BPD))) 

#run mediation analysis
m.out = mediate(model_mediator, model_outcome, sims=1000, treat="Mom_steroids", mediator="Gestational_age_centered", robustSE=T)
summary(m.out)
plot(m.out)

#proportion mediated for inconsistent mediation
abs(m.out$d.avg) / (abs(m.out$d.avg)+abs(m.out$z.avg)) * 100

