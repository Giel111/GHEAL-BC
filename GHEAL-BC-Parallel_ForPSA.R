# Intro ---
## Outline including all branches,
##    No parameters yet
#setwd('C:/Users/gielv/OneDrive/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel.')
# Libraries ----
#workdir
#setwd('C:/Users/gielv/OneDrive/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel') #<-home
#setwd('Z:/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel')                      #<-work

#library(doParallel)
#registerDoParallel(cores=detectCores(all.tests = FALSE, logical = TRUE)-1)

rm(list=ls())
gc()

StratName = 'DiagDisUtilHigh'
set.seed(111)

library(parallel);
library(doSNOW);
library(pbapply)
library(tidyr)
library(dplyr)


getSingleAttribute <- function(attribute, output, all=F) sapply(unique(output[,"name"]), function(entity) {
  list(
    if(all) {
      output[output[,"name"]==entity & output[,"key"]==attribute, "value"]
    } else {
      tail(output[output[,"name"]==entity & output[,"key"]==attribute, "value"], n=1)
    }
  )
})


RunSimPar <- function(n.patients,n.runs,free.cores=1){
  #set.seed(round(runif(1,1,1000)))
  #set.seed(111)
  #library(parallel)


  cl <- makeCluster(detectCores()-free.cores);
  registerDoSNOW(cl);
  clusterExport(cl, c("getSingleAttribute"));
  
  results <- pbsapply(cl=cl,X=1:n.runs,FUN=function(run){
    
    # PASTE CODE HERE FROM OTHER
    
    library(simmer)
    library(simmer.plot)
    library(simmer.bricks)
    
    library(fitdistrplus)
    library(dplyr)
    library(pracma)
    library(mvtnorm)
    library(ggplot2)
    library(gridExtra)
    #set.seed(9)
    
    
    # STRATEGY SELECTOR ----
    
    # here you can select the screening strategy for this simulation run.
    Start_screen_age <- 50 
    end_screen_age <- 75
    screen_interval <- 2 
    
    manualScreenInput <- F
    manualScreenAges <- c(50,51,52,53,54,55)
    if (manualScreenInput){
      # enter ages at which patient should be screened here:
      ScreenAges <- manualScreenAges
    }
    
    
    
    # Functions ---- 
    
    ## Cancer ----
    LifeExpAndCancerAge <- function(){
      # Generate one random value from the copula
      random_cancer <- 2
      random_life_exp <- 1
      if (runif(1)<(0.15)){ # add a / between 1 and 7 ;)
        while (random_life_exp<random_cancer){
          random_sample <- rmvnorm(1, mean = c(0, 0), sigma = correlation_matrix)
          
          # Transform copula data to distribution data
          random_life_exp <- qweibull(pnorm(random_sample[1,1]), 
                                      shape = LifeExpShape, 
                                      scale = LifeExpScale) *year
          random_cancer <- qnorm(pnorm(random_sample[1,2]), 
                                 mean = CancerExpMean, 
                                 sd = CancerExpSD)  *year
        }
      } else {
        random_life_exp <- rweibull(1,LifeExpShape,LifeExpScale) *year
        random_cancer <- 1000*year
      }
      set <- c(random_life_exp,random_cancer)
      return(set)
    }
    
    
    GompGrow <- function(age,onsetage,GR,IsCured){
      if(IsCured ==1){return(0)} else{
        t <- (age-onsetage ) / month
        if (t<0){
          size<-0
        } else {
          #Determine size of tumour MANC
          Volume <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*t))^4 
          #tumour volume at time t
          size <- 2*(Volume/(4/3*pi))^(1/3)
        }
        return(size)
      }
    }
    
    FindTimeAtSize <- function(size, onsetage, GR, max_size) {
      
      # Solve for the time delta using a binary search algorithm
      low <- 0
      high <- 1000*year  # set an arbitrarily high upper bound for time delta
      while (high - low > 1e-6) {
        mid <- (low + high) / 2
        #size_mid <- GompGrow(onsetage + mid, onsetage, GR, IsCured=0)
        Volume_mid <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*mid))^4 
        #tumour volume at time t
        size_mid <- 2*(Volume_mid/(4/3*pi))^(1/3)
        if (size_mid < size) {
          low <- mid
        } else {
          high <- mid
        }
      }
      
      # Return the estimated time when the given size was reached
      return(onsetage + high)
    }
    
    GompGrow2 <- function(age,onsetage,GR,IsCured,
                          max_size,RegressionSize,StagnateSize){
      age = age/month
      onsetage=onsetage/month
      if (IsCured ==1){return(0)}
      Vm = (4/3)*pi*(max_size/2)**3 
      size <- 0
      delta = age-onsetage
      if (delta<0){
        size <-0
        return (size)
      } 
      #Determine size of tumour MANC
      Volume <- Vm/(1+((Vm/Vc)^0.25-1)*exp(-0.25*GR*delta))^4 
      #tumour volume at time t
      size <- 2*(Volume/(4/3*pi))^(1/3)
      if (size>StagnateSize){
        size <- StagnateSize
      }
      if (size>RegressionSize){
        size <- RegressionSize
        time_reached = FindTimeAtSize(RegressionSize,onsetage,GR,max_size)
        delta2 = delta-time_reached
        size2 = size - (size*pnorm(delta2,50,25))
        return(size2)
      }
      return(size)
    }
    
    RegStagSizes <- function(){
      p <- runif(1)
      RS = 1000*year
      SS=1000*year
      if (p < 0.05){
        RS = rnorm(1,25,3)
        
      } else if (p<0.10){
        SS = rnorm(1,25,3)
        
      }
      
      return(c(RS,SS))
    }
    
    Staging <- function(size, p) {
      size_bin <- SizeTable[SizeTable$TumBinL <= size & SizeTable$TumBinR > size,]
      if (nrow(size_bin) == 0) {
        return(0)
      }
      cumulative_p <- cumsum(size_bin[,-c(1,2)])
      col_index <- which(cumulative_p >= p)[1]
      if (col_index == 1){
        return(0.5)
      } else if(col_index==2){
        return(1)
      } else if(col_index==3){
        return(2)
      } else if(col_index==4){
        return(3)
      } else if(col_index==5){
        return(4)
      }
      #return(colnames(SizeTable)[col_index + 2])
    }
    
    CancerSurvival <- function(stage,lifeExp,currAge){
      p = runif(1)
      if (stage == 0.5){
        NLE = lifeExp #NLE= new life exp
      } else if (stage == 4) {
        NLE = tail(which(SurvivalDf$Stage_IV >= p * 100) - 1,1)*year + currAge
      } else if (stage == 3) {
        NLE = tail(which(SurvivalDf$Stage_III >= p * 100) - 1,1)*year+ currAge
      } else if (stage == 2) {
        NLE = tail(which(SurvivalDf$Stage_II >= p * 100) - 1,1)*year+ currAge
      } else if (stage == 1) {
        NLE = tail(which(SurvivalDf$Stage_I >= p * 100) - 1,1)*year+ currAge
      } else{NLE=lifeExp}
      
      if (NLE  > lifeExp){NLE <- lifeExp}
      
      
      return(NLE)
    }
    
    
    
    ScreenResult <- function(TumorSize){
      min_dia <- rweibull(1,MeanMinDetSize,SDMinDetSize) # from miscan paper
      # define BIRADS based on how much above/below this size
      if (TumorSize ==0){
        birads = ifelse(runif(1)<FP_Perc,3,1) #False Positive
      } else if (TumorSize<0.95*min_dia){
        birads=1
      } else if (TumorSize>0.95*min_dia & TumorSize<1.25*min_dia){
        birads = 3
      } else if (TumorSize>1.25*min_dia & TumorSize<2*min_dia){
        birads = 4
      } else if (TumorSize>2*min_dia){
        birads = 5
      }
      return(birads)
    }
    
    ClinicalCheckAtHome <- function(TumorSize,ClinSizeCheck,PatientAge){
      if (PatientAge %in% ScreenAges){return(0)}
      if ((PatientAge/month) %% 2 == 0){ #round(runif(1,1,2))
        if (TumorSize>ClinSizeCheck){
          return(1)
        } else {return(0)}
      }else {return(0)}
    }
    
    
    GrowRater <- function(startAge){
      GRmultiplier <- 0
      GR <- rgamma(1,shape = Grow_Gamma_shape, rate = Grow_Gamma_rate)
      # print(GR)
      startAge <- startAge/year
      if (startAge <50){
        GRmultiplier = -(62-startAge)/80
      }
      if (startAge >75){
        GRmultiplier <- -(startAge-62)/80
      }
      GR <- GR + GRmultiplier
      if (GR<0){
        GR <- GR - GRmultiplier
      }
      return(GR)
    }
    
    ## Time ----
    WaitAtHome <- function(EntryAge){
      # find next month for now
      div <- EntryAge %% month
      return(month - div)
    }
    
    ScreenTime <- function(){
      return(max(0.1,rnorm(1,MeanScreenTime,SDScreenTime)))
    }
    DiagTime <- function(){
      return(max(0.1,rnorm(1,MeanDiagTime,SDDiagTime)))
    }
    
    
    ## Next Steps ----
    AfterHome <- function(Invite,Screened,Birads,ClinicalCheck,
                          Referral,PatientAge,HealthyLifeExpectancy,TumorSize){
      # Next Steps update
      #All five steps possible:
      #a. Go to screening after invite-> 2
      #b. Go to diagnostics due to tumor size or referral after screening ->3
      #c. Go to hospital due to referral -> 4
      #d. Die due to age or illness-> 5
      #e. Continue in Home trajectory, nothing's needed -> 1 
      #PatientAge <- PatientAge / year
      s <- 0
      if (Invite>0){
        if (PatientAge %in% ScreenAges){ # check if go for screen
          if (Invite==1){
            upt = uptakefirstscreen
          }else if (Invite > 1){
            upt = uptakeotherscreen
          } else {upt = uptakenoscreen}
          if (upt > runif(1)){
            s <- 2
          } 
        }
      }
      
      if (ClinicalCheck ==1 | Birads == 3 | Birads == 4 | Birads ==5){ 
        # check if diagnsostics 
        s <- 3
      }
      
      if (Referral ==1){
        s <-4
      }
      
      if (PatientAge > HealthyLifeExpectancy){
        s <- 5 # die
        if (TumorSize>Vm){print('this should not happen')}
      }
      if (s==0){s<-1}# stay at home
      return(s)
    }
    
    ## Utilities ----
    BaseAgeUtil <- function(Age){
      Age <- Age/year
      if (Age < 31){
        return(1)
      }
      Age <- ceiling((Age-30)/5) # steps of 5
      return(utility_ages[Age,2])
    }
    
    CancerBasedUtil <- function(currentUtil,stage,age,startAge){
      if (stage==0){
        s=1
      } else if (stage==0.5){
        s=2 
      } else if (stage==1){
        s=3
      } else if (stage==2){
        s=4
      } else if (stage ==3){
        s=5
      } else if (stage==4){
        s=6
      }
      delta <- (age-startAge) / year
      utilDec <- utility_decrements[s,2]
      discUtilDec <- utilDec/((1+UtilDiscount)^floor(delta))
      newUtil = currentUtil * (1-discUtilDec)
      return(newUtil)
    }
    
    ## Costs----
    
    ScreenDiagInvitesCost <- function(screened,diaged,invite){
      if (is.na(screened)){screened=0}
      if (is.na(diaged)){diaged=0}
      if (is.na(invite)){invite=0}
      cost<- 0 
      cost<- cost + screened* AvgScreenCost
      cost<- cost + diaged*AvgDiagCost
      cost<- cost + invite*AvgInvCost
      return(cost)
    }
    
    
    TreatmentCost <- function(stage,age,lifeExp){
      # cost is determined by:
      # initial cost for first 12 months post treatment
      # ^always counted
      # costs per month after the first 12, up to the last 6 months
      # counted if this interval exists
      # Costs for the last six months of care
      # ^always counted
      delta <- (lifeExp - age) / month
      if (delta<0){return(0)}
      #delta is how long we have left
      c<-0 # cost variable
      #determine right row from table:
      if (stage==0.5){
        row <- 2
      } else if (stage ==0){
        row <- 1
      } else {
        row<- stage+2
      } 
      InitC <- rnorm(1,TreatCostTable$Initial12Cost[row],
                     TreatCostTable$InitialSD[row]) /12 #round to monthly
      TerminalC <- rnorm(1,TreatCostTable$Terminal6Costs[row],
                         TreatCostTable$TerminalSD[row]) /6 # to monthly
      MidC <- rnorm(1,TreatCostTable$ContinuousMonthlyCost[row],
                    TreatCostTable$ContinuousSD[row]) /2
      
      CInit <- min(delta,12) * InitC #no discount, this is year 0
      CTerminal <- min(delta,6) * TerminalC / 
        (1 + CostDiscount)^((round((lifeExp - age)/year)) - 1) 
      # discount with 4% per year post diagnosis
      CMid <- annuity_cost_monthly(MidC,delta,CostDiscount)
      
      c <- CInit + CTerminal + CMid
      return(c)
      
    }
    
    
    
    VisitationCosts<- function(screened,diaged,hosps){
      costs<- 0
      hosps <- hosps/day
      #travel costs
      costs <- costs + 
        (ifelse(runif(1)<perc_car,1,0) * dist_to_screen * cost_per_km) * screened 
      # for screening visits
      costs <- costs + 
        (ifelse(runif(1)<perc_car,1,0) * dist_to_hosp * cost_per_km) * diaged 
      # for diagnosis in hospital
      costs <- costs + 
        (ifelse(runif(1)<perc_car,1,0) * dist_to_hosp * cost_per_km) * hosps
      
      # productivity loss costs:
      # screening & diag takes average time, hospital takes full working day
      costs <- costs + (screened*MeanScreenTime*prod_cost_hour) + 
        (diaged*MeanDiagTime*prod_cost_hour) + 
        (hosps*8*hour*prod_cost_hour)
      
      return(costs)
    }
    
    PaidCosts <- function(lifeExp,StartAge){
      costs<- 0
      # PAID tool gives insight in expected rest of life costs for if someone 
      #lives longer after cancer treatmenet.
      if (lifeExp>StartAge){
        costlist <- WomenCostList[(round(StartAge)/year):(round(lifeExp)/year)]
        yearlist <- 0:(length(costlist)-1)
        disclist<- costlist / (1+0.04)^yearlist
        costs = sum(disclist)
      }
      return(costs)
    }
    
    annuity_cost_monthly <- function(payment, months, rate) {
      monthly_rate <- (1 + rate)^(1/12) - 1
      total_cost <- 0
      if (months<13){
        return(0)
      }
      months_left = min(0,months-6) #don't count the final 6
      for (i in 13:months_left) {
        discounted_payment <- payment / (1 + monthly_rate)^(i - 1)
        total_cost <- total_cost + discounted_payment
      }
      return(total_cost)
    }
    
    ## Sim ----
    getSingleAttribute <- function(attribute, 
                                   output, 
                                   all=F) sapply(unique(output[,"name"]), 
                                                 function(entity) {
                                                   list(
                                                     if(all) {
                                                       output[output[,"name"]==entity & output[,"key"]==attribute, "value"]
                                                     } else {
                                                       tail(output[output[,"name"]==entity & 
                                                                     output[,"key"]==attribute, "value"], n=1)
                                                     }
                                                   )
                                                 })
    
    getMultipleAttributes <- function(attributes, output) as.data.frame(sapply(attributes, function(attribute) as.numeric(getSingleAttribute(attribute, output))))
    
    
    
    # Parameters ----
    
    ## Patients----
    n.patients=500;
    mon.patients <- ifelse(n.patients<100,2,0);  
    
    LifeExpShape <- 7.937
    LifeExpScale <- 86.788
    # Transform copula data to distribution data
    
    
    ## Data loading ----
    correlation_matrix <- matrix(c(1, 0.5, 0.5, 1), ncol = 2) 
    #for age & cancer incidence
    
    SizeTable <- read.csv('StagingFromSizes.csv')
    ## Time-----
    hour = 1;
    minute = hour/60;
    second=minute/60;
    day=hour*24;
    year = day*365; 
    month = year/12 #not exactly but fine
    
    
    MeanScreenTime <- 4*hour
    SDScreenTime <- 30*minute 
    MeanDiagTime <- 6*hour
    SDDiagTime <- 30*minute
    
    
    ## Utilities ----
    # these are base utilities per age-range, 
    # taken from the MANC model
    utility_ages<-data.frame(c(30,35,40,45,50,
                               55,60,65,70,75,
                               80,85,90,95,100),
                             c(0.9383,0.9145,0.9069,0.8824,0.8639,
                               0.8344,0.8222,0.8072,0.8041,0.779,
                               0.7533,0.6985,0.6497,0.6497,0.6497))
    
    utility_decrements <- data.frame(c('Healthy','DCIS','StageI',
                                       'StageII','StageIII','StageIV'),
                                     c(0,0.1,0.15,0.20,0.25,0.40))
    UtilDiscount <- 0.015
    ScreenDisUtil <- 0.15
    DiagDisUtil <- 0.5
    ## Costs ---- 
    
    TreatCostTable <- read.csv('TreatmentCostsVanLuijt2023Eur.csv')
    CostDiscount <- 0.04
    
    AvgScreenCost <- 100
    AvgInvCost <- 5
    AvgDiagCost <- 212
    
    ## Tumor ----
    #Set tumour growth rate parameters
    CancerExpMean = 61.871
    CancerExpSD = 14.14
    
    Grow_Gamma_shape <- 1.567742
    Grow_Gamma_rate <- 1.933883
    
    
    max_size <- 128 #mm diameter
    start_size <- 0.25 #starting size of tumours, diameter in mm
    Vc = (4/3)*pi*(start_size/2)^3 #Volume at start
    Vm = (4/3)*pi*(max_size/2)^3 #Max volume
    
    #DF to show survival probability after X years at X stage
    SurvivalDf <- data.frame(Years_After_Diagnosis = c(0:10),
                             Stage_I = c(100, 100, 100, 99, 99, 
                                         98, 97, 97, 96, 95, 95),
                             Stage_II = c(100, 99, 97, 95, 93, 
                                          91, 89, 88, 86, 85, 83),
                             Stage_III = c(100, 96, 89, 83, 78, 
                                           73, 69, 65, 62, 60, 58),
                             Stage_IV = c(100, 70, 53, 38, 29, 
                                          22, 17, 13, 11, 9, 7))
    
    ## Detection ----
    MeanMinDetSize <- 0.91 # for screening
    SDMinDetSize <- 0.34 # for screening
    
    MeanLogClin <-  3.22 
    SdLogClin <- 2.25
    AlwaysFoundSize<- 100 # buld in a min! so min(100,clindiagsize)
    
    FP_Perc <- 0.05 #FalsePositive
    
    
    ## Screen Start Params ----
    Start_screen_age <- Start_screen_age *year
    end_screen_age <- end_screen_age * year
    screen_interval <- screen_interval*year
    ScreenAges = seq(Start_screen_age,end_screen_age,by=screen_interval)
    ScreenAges <- append(ScreenAges,1000*year)
    if (manualScreenInput){
      # enter ages at which patient should be screened here:
      ScreenAges <- manualScreenAges * year
    }
    # ScreenAges = c(50,52,54,56,58,60,62,64,66,68,70,72,74)*year
    
    #Screening uptake MANC model
    uptakefirstscreen<- 0.78#0.605
    uptakeotherscreen<-0.9#0.852
    uptakenoscreen<-0.25#0.191
    
    ## Societal Cost Attributes ----
    dist_to_hosp <- 7
    dist_to_screen <- 1.1
    avg_park_cost <- 3
    perc_car <- 0.8 # msot go by car, definetly for hospital visits
    
    cost_per_km <- 0.19 # for Public transport & car
    
    prod_cost_hour <- 31.6 # for women, dutch guidelines
    
    AdditionalLifeYearCosts <- read.csv('PAID_Cost_additional_year_post_BC_Costs_Living_Year_Longer_2023-03-07.csv')
    WomenCostList <- AdditionalLifeYearCosts$Unrelated_Women
    
    # Attributes Recording ----
    allAttributes = c("Index","Name","Alive", "PatientAge", "HealthyLifeExpectancy", 
                      "TumorStartAge", "LastAge", "TotalCosts", 
                      "MedicalCosts", "SocietalCosts","PAIDCosts", "TotalUtility", 
                      "CurrentUtility", "TumorSize", "WorstCancer", 
                      "WorstStage", "CancerStageT", "IsCured", "StagingProb", 
                      "TumorGrowthRate", "max_size", 
                      "RegressionSize", "StagnateSize", "ClinSizeCheck", 
                      "Invite", "Screened", "BIRADS", 
                      "Referral", "NextStep", "EntryAge", "TimeAtHome", 
                      "ClinicalCheck", "TimeAtDiagnostics", 
                      "TotalTimeAtDiagnostics", "DiagnosedThrough", 
                      "DiagnosticVisits", "TimeAtHospital", "TotalTimeAtHospital", 
                      "OriginalLifeExp", "HospitalVisits", "TimeAtScreening",
                      "TotalTimeAtScreening", "FirstStage")
    
    EndDatDf <- data.frame(matrix(rep('test',length(allAttributes)),
                                  ncol=length(allAttributes)))
    colnames(EndDatDf) <- allAttributes
    addAtts <- function(EndDatDf){
      addList <- c(nrow(EndDatDf),
                   get_name(basic_sim),
                   get_attribute(basic_sim,'Alive'),
                   get_attribute(basic_sim,'PatientAge'),
                   get_attribute(basic_sim,'HealthyLifeExpectancy'),
                   get_attribute(basic_sim,'TumorStartAge'),
                   get_attribute(basic_sim,'LastAge'),
                   get_attribute(basic_sim,'TotalCosts'),
                   get_attribute(basic_sim,'MedicalCosts'),
                   get_attribute(basic_sim,'SocietalCosts'),
                   get_attribute(basic_sim,'PAIDCosts'),
                   get_attribute(basic_sim,'TotalUtility'),
                   get_attribute(basic_sim,'CurrentUtility'),
                   get_attribute(basic_sim,'TumorSize'),
                   get_attribute(basic_sim,'WorstCancer'),
                   get_attribute(basic_sim,'WorstStage'),
                   get_attribute(basic_sim,'CancerStageT'),
                   get_attribute(basic_sim,'IsCured'),
                   get_attribute(basic_sim,'StagingProb'),
                   get_attribute(basic_sim,'TumorGrowthRate'),
                   get_attribute(basic_sim,'max_size'),
                   get_attribute(basic_sim,'RegressionSize'),
                   get_attribute(basic_sim,'StagnateSize'),
                   get_attribute(basic_sim,'ClinSizeCheck'),
                   get_attribute(basic_sim,'Invite'),
                   get_attribute(basic_sim,'Screened'),
                   get_attribute(basic_sim,'BIRADS'),
                   get_attribute(basic_sim,'Referral'),
                   get_attribute(basic_sim,'NextStep'),
                   get_attribute(basic_sim,'EntryAge'),
                   get_attribute(basic_sim,'TimeAtHome'),
                   get_attribute(basic_sim,'ClinicalCheck'),
                   get_attribute(basic_sim,'TimeAtDiagnostics'),
                   get_attribute(basic_sim,'TotalTimeAtDiagnostics'),
                   get_attribute(basic_sim,'DiagnosedThrough'),
                   get_attribute(basic_sim,'DiagnosticVisits'),
                   get_attribute(basic_sim,'TimeAtHospital'),
                   get_attribute(basic_sim,'TotalTimeAtHospital'),
                   get_attribute(basic_sim,'OriginalLifeExp'),
                   get_attribute(basic_sim,'HospitalVisits'),
                   get_attribute(basic_sim,'TimeAtScreening'),
                   get_attribute(basic_sim,'TotalTimeAtScreening'),
                   get_attribute(basic_sim,'FirstStage')
      )
      EndDatDf<<- rbind(EndDatDf,addList)
      return(0)
    }
    
    
    # Simulation ----
    
    ## Trajectories ----
    
    ### Init ----
    # trajectory for setting all patient-specific  starting attributes
    
    Initialization <- trajectory()%>%
      set_attribute(key='Alive',value=1) %>%
      set_attribute(key='PatientAge',value=0) %>%
      set_attribute(key=c('HealthyLifeExpectancy','TumorStartAge'),
                    value= function() LifeExpAndCancerAge())%>%
      set_attribute(key='LastAge',value=0)%>%
      set_attribute(key='TotalCosts',value=0)%>%
      set_attribute(key='MedicalCosts',value=0) %>%
      set_attribute(key='SocietalCosts',value=0) %>%
      set_attribute(key='TotalUtility',value=0)%>%
      set_attribute(key='CurrentUtility',value=1) %>%
      set_attribute(key='TumorSize',value=0) %>%
      set_attribute(key='WorstCancer',value=0) %>%
      set_attribute(key='WorstStage',value=0) %>%
      set_attribute(key='CancerStageT',value=0) %>%
      set_attribute(key='IsCured',value=0) %>%
      set_attribute(key='StagingProb',value=function() runif(1)) %>%
      set_attribute(key = 'TumorGrowthRate', 
                    value = function() GrowRater(get_attribute(basic_sim,
                                                               'TumorStartAge')))%>%
      set_attribute(key = 'max_size',
                    value=function() rnorm(1,max_size,10)) %>%
      set_attribute(key = c('RegressionSize','StagnateSize'),
                    value= function() RegStagSizes())%>%
      set_attribute(key='ClinSizeCheck',
                    value=function() min(AlwaysFoundSize,
                                         rlnorm(1,MeanLogClin,SdLogClin))) %>% 
      # set size needed for clinical find for this patient
      set_attribute(key = 'Invite',value=0)%>%
      set_attribute(key = 'Screened',value=0)%>%
      set_attribute(key = 'BIRADS',value=10) %>% # Invalid BIRADS value before checkup
      set_attribute(key = 'Referral',value=0) %>%# No referral
      set_attribute(key='NextStep',value=1)  # go home for first iteration
    #timeout(30*year) 
    #log_(paste0('Going Home'))
    
    ### Home ----
    VisitHome <- trajectory()%>%
      seize(resource='Home')%>%
      set_attribute(key= 'EntryAge',
                    value=function() now(basic_sim)) %>% # record entry age
      set_attribute(key='TimeAtHome',
                    value=function() WaitAtHome(
                      get_attribute(basic_sim,'EntryAge'))) %>% 
      # determine how long at home (either next step or until next month)
      # time update
      timeout(function() get_attribute(basic_sim,'TimeAtHome')) %>%
      set_attribute(key = 'PatientAge',
                    value = function() get_attribute(basic_sim,'TimeAtHome'),
                    mod='+',init=0) %>%
      # Next Steps update
      #All five steps possible:
      #a. Go to screening after invite-> 2
      #b. Go to diagnostics due to tumor size or referral after screening ->3
      #c. Go to hospital due to referral -> 4
      #d. Die due to age or illness-> 5
      #e. Continue in Home trajectory, nothing's needed -> 1 
      set_attribute(key = 'Invite',
                    value= function() ifelse(get_attribute(basic_sim,'PatientAge') 
                                             %in% ScreenAges,1,0),
                    mod='+')%>% # point a  
      set_attribute(key='ClinicalCheck',
                    value=function() ClinicalCheckAtHome(
                      get_attribute(basic_sim,'TumorSize'),
                      get_attribute(basic_sim,'ClinSizeCheck'),
                      get_attribute(basic_sim,'PatientAge'))) %>%
      
      set_attribute('NextStep',
                    value=function() AfterHome(
                      get_attribute(basic_sim,'Invite'),
                      get_attribute(basic_sim,'Screened'), # to check point 
                      get_attribute(basic_sim,'BIRADS'), 
                      get_attribute(basic_sim,'ClinicalCheck'),# to check point b                                            
                      get_attribute(basic_sim,'Referral'),  # point c                                            
                      get_attribute(basic_sim,'PatientAge'), 
                      get_attribute(basic_sim,'HealthyLifeExpectancy'), 
                      get_attribute(basic_sim,'TumorSize') # point D
                    )) %>%  
      
      # Update tumor & cancer parameters:
      set_attribute('TumorSize', 
                    value=function() GompGrow2(
                      get_attribute(basic_sim,"PatientAge"),
                      get_attribute(basic_sim,"TumorStartAge"),
                      get_attribute(basic_sim,'TumorGrowthRate'),
                      get_attribute(basic_sim,'IsCured'),
                      get_attribute(basic_sim,'max_size'),
                      get_attribute(basic_sim,'RegressionSize'),
                      get_attribute(basic_sim,'StagnateSize'))) %>%
      set_attribute('WorstCancer',
                    value = function() ifelse(
                      get_attribute(basic_sim,'TumorSize')>
                        get_attribute(basic_sim,'WorstCancer'),
                      get_attribute(basic_sim,'TumorSize'),
                      get_attribute(basic_sim,'WorstCancer'))) %>%
      
      # Update utilities
      set_attribute(key='CurrentUtility',
                    value = function() BaseAgeUtil(
                      get_attribute(basic_sim,'PatientAge'))) %>%
      set_attribute(key='CurrentUtility',
                    value = function() CancerBasedUtil(
                      get_attribute(basic_sim,'CurrentUtility'),
                      get_attribute(basic_sim,'WorstStage'),
                      get_attribute(basic_sim,'PatientAge'),
                      get_attribute(basic_sim,'TumorStartAge'))) %>%
      set_attribute(key='TotalUtility',
                    value=function() get_attribute(
                      basic_sim,'TimeAtHome')*
                      get_attribute(basic_sim,'CurrentUtility'),
                    mod='+')%>%
      
      
      
      release(resource='Home')
    
    ### Screening ----
    VisitScreening <- trajectory()%>%
      seize('Screening') %>%
      set_attribute(key= 'EntryAge',value=function() now(basic_sim)) %>% 
      # record entry age
      set_attribute(key='TimeAtScreening',value= function() ScreenTime()) %>% 
      # determine how long at home (either next step or until next month)
      # time update
      timeout(function() get_attribute(basic_sim,'TimeAtScreening')) %>%
      set_attribute(key = 'PatientAge',
                    value = function() get_attribute(basic_sim,'TimeAtScreening'),
                    mod='+',
                    init=0) %>%
      set_attribute(key= 'TotalTimeAtScreening',
                    value = function() get_attribute(basic_sim,'TimeAtScreening'),
                    mod='+',
                    init=0) %>%
      
      # add utility
      set_attribute(key = 'CurrentUtility',
                    value= function() get_attribute(
                      basic_sim,'CurrentUtility') - ScreenDisUtil) %>% 
      # disutility for during screening 
      set_attribute(key='TotalUtility',
                    value=function() get_attribute(
                      basic_sim,'TimeAtScreening')*
                      get_attribute(basic_sim,'CurrentUtility'),
                    mod='+')%>%
      # perform screening
      set_attribute(key = 'BIRADS',
                    value =function() ScreenResult(
                      get_attribute(basic_sim,'TumorSize'))) %>% 
      set_attribute(key='Screened',value=1,mod='+',init=0)%>%
      #check next step
      set_attribute(key = 'NextStep',value= 1)%>% #go home after screening
      
      release('Screening')
    
    ### Diagnostics ----
    VisitDiagnostics <- trajectory()%>%
      seize('Diagnostics') %>%
      set_attribute(key= 'EntryAge',value=function() now(basic_sim)) %>% 
      # record entry age
      set_attribute(key='TimeAtDiagnostics',value= function() DiagTime()) %>% 
      # determine how long at home (either next step or until next month)
      # time update
      timeout(function() get_attribute(basic_sim,'TimeAtDiagnostics')) %>%
      set_attribute(key = 'PatientAge',
                    value = function() get_attribute(basic_sim,'TimeAtDiagnostics'),
                    mod='+',
                    init=0) %>%
      set_attribute(key= 'TotalTimeAtDiagnostics',
                    value = function() get_attribute(basic_sim,'TimeAtDiagnostics'),
                    mod='+',init=0)%>%
      # add utility
      set_attribute(key = 'CurrentUtility',
                    value= function() get_attribute(
                      basic_sim,'CurrentUtility') - DiagDisUtil) %>% #disutility for diag 
      set_attribute(key='TotalUtility',
                    value=function() get_attribute(
                      basic_sim,'TimeAtDiagnostics')*
                      get_attribute(basic_sim,'CurrentUtility'),
                    mod='+')%>%
      set_attribute(key = 'DiagnosedThrough',
                    value = function() ifelse(
                      get_attribute(basic_sim,'ClinicalCheck'),1,0)) %>% 
      #through screening or clinical?
      set_attribute(key = "CancerStageT",
                    value= function() Staging(
                      get_attribute(basic_sim,'TumorSize'),
                      get_attribute(basic_sim,'StagingProb'))) %>%
      set_attribute(key = "WorstStage",
                    value= function() ifelse(
                      get_attribute(basic_sim,'CancerStageT')>get_attribute(basic_sim,'WorstStage'),
                      get_attribute(basic_sim,'CancerStageT'),
                      get_attribute(basic_sim,'WorstStage'))) %>%
      set_attribute('DiagnosticVisits',value = 1, mod= '+',init=0) %>%
      set_attribute(key = 'NextStep',
                    value= function() ifelse(
                      get_attribute(basic_sim,'CancerStageT')==0,
                      1,
                      4))%>% 
      release('Diagnostics')
    
    ### Cured ----
    Cured <- trajectory() %>%
      #log_('cured') %>%
      set_attribute(key='TumorSize', value = 0) %>%
      set_attribute(key='IsCured',value=1) %>%
      set_attribute(key='BIRADS',value=10) %>%
      set_attribute(key='CancerStageT',value=0) %>%
      set_attribute(key = 'NextStep',value=1)
    
    ### Residuals ----
    Residual <- trajectory() %>%
      #log_('residual') %>%
      set_attribute('TumorSize', value= 35) %>%
      set_attribute(key='IsCured',value=0) %>%
      set_attribute(key='BIRADS',value=4) %>%
      set_attribute(key='FirstStage',
                    value=function() get_attribute(basic_sim,'WorstStage')) %>%
      set_attribute(key='CancerStageT',
                    value=function() get_attribute(basic_sim,'WorstStage')) %>%
      
      # go back home
      
      set_attribute(key = 'NextStep',value=3)
    
    ### Hospital ----
    VisitHospital <- trajectory()%>%
      seize('Hospital') %>%
      set_attribute(key= 'EntryAge',
                    value=function() now(basic_sim)) %>% # record entry age
      set_attribute(key='TimeAtHospital',
                    value= function() 
                      3^get_attribute(basic_sim,'CancerStageT')*day)%>% 
      # determine how long at home (either next step or until next month)
      # time update
      timeout(function() get_attribute(basic_sim,'TimeAtHospital')) %>%
      set_attribute(key = 'PatientAge',
                    value = function() get_attribute(basic_sim,'TimeAtHospital'),
                    mod='+',
                    init=0) %>%
      set_attribute(key= 'TotalTimeAtHospital',
                    value = function() get_attribute(basic_sim,'TimeAtHospital'),
                    mod='+',
                    init=0) %>%
      
      # update utilities
      set_attribute(key='CurrentUtility',
                    value = function() BaseAgeUtil(
                      get_attribute(basic_sim,'PatientAge'))) %>%
      set_attribute(key='CurrentUtility',
                    value = function() CancerBasedUtil(
                      get_attribute(basic_sim,'CurrentUtility'),
                      get_attribute(basic_sim,'CancerStageT'),
                      get_attribute(basic_sim,'PatientAge'),
                      get_attribute(basic_sim,'TumorStartAge'))) %>%
      set_attribute(key='TotalUtility',
                    value=function() get_attribute(
                      basic_sim,'TimeAtHospital')*
                      get_attribute(basic_sim,'CurrentUtility'),
                    mod='+')%>%
      #instead of cure, go to adjust life exp
      set_attribute(key='OriginalLifeExp', 
                    value=function() get_attribute(
                      basic_sim,'HealthyLifeExpectancy')) %>%
      set_attribute(key = 'HealthyLifeExpectancy',
                    value=function() CancerSurvival(
                      get_attribute(basic_sim,'WorstStage'),
                      get_attribute(basic_sim,'HealthyLifeExpectancy'),
                      get_attribute(basic_sim,'TumorStartAge'))) %>%
      set_attribute(key='MedicalCosts',
                    value=function() TreatmentCost(
                      get_attribute(basic_sim,'WorstStage'),
                      get_attribute(basic_sim,'PatientAge'),
                      get_attribute(basic_sim,'HealthyLifeExpectancy')),
                    mod='+') %>%
      branch(option= function() ifelse(runif(1)<0.88,1,2),continue=c(T,T),
             Cured,
             Residual) %>%
      
      release('Hospital')
    
    
    ### Death ----
    Death <- trajectory()%>%
      #log_('die') %>%
      set_attribute(key='Alive',value=0)%>%
      
      #compute total costs
      set_attribute('MedicalCosts',
                    value = function() ScreenDiagInvitesCost(
                      get_attribute(basic_sim,'Screened'),
                      get_attribute(basic_sim,'DiagnosticVisits'),
                      get_attribute(basic_sim,'Invite')),
                    mod = '+') %>%
      set_attribute(key = 'PAIDCosts',
                    value = function() PaidCosts(
                      get_attribute(basic_sim,'HealthyLifeExpectancy'),
                      get_attribute(basic_sim,'TumorStartAge')),
                    mod = '+') %>%
      set_attribute(key='SocietalCosts',
                    value=function() VisitationCosts(
                      get_attribute(basic_sim,'Screened'),
                      get_attribute(basic_sim,'DiagnosticVisits'),
                      get_attribute(basic_sim,'TotalTimeAtHospital'))) %>%
      
      
      
      # extract all atrributes
      set_attribute(key='Alive',value = function() addAtts(EndDatDf))
    
    ### Base Model ----
    #stitching the blocks together
    basic_model <- trajectory()%>%
      # first initialize
      join(Initialization) %>%
      
      # 'store' patients at home until it's:
      # time to go to screening
      # time to go to diagnostics
      # time to go to hospital
      # time to die
      branch(option=function() get_attribute(basic_sim,'NextStep'), 
             continue=c(T,T,T,T,F),
             VisitHome,
             VisitScreening,
             VisitDiagnostics,
             VisitHospital,
             Death) %>%
      
      
      rollback(amount=1)
    
    
    ## Plot ----
    #plot(basic_model,verbose=TRUE)
    
    # plot(Initialization,verbose=T)
    # plot(VisitHome,verbose=T)
    # plot(VisitScreening,verbose=T)
    # plot(VisitDiagnostics,verbose=T)
    # plot(Cured,verbose=T)
    # plot(Residual,verbose=T)
    # plot(VisitHospital,verbose=T)
    
    ## Run----
    
    basic_sim <- simmer()%>%
      add_resource('Home',capacity=Inf) %>%
      add_resource('Screening',capacity=Inf) %>%
      add_resource('Diagnostics',capacity=Inf) %>%
      add_resource('Hospital',capacity=Inf) %>%
      
      add_generator(name_prefix='patient',
                    trajectory=basic_model,
                    distribution=at(rep(x=0, times=n.patients)),
                    mon=mon.patients)
    
    
    start_time <- Sys.time()   
    basic_sim %>%
      reset() %>%
      run(progress=progress::progress_bar$new()$update,until=110*year,steps=1000)
    
    
    patient_monitor <-
      get_mon_arrivals(basic_sim) %>%
      transform(wait = end_time - start_time - activity_time) %>%
      transform(totalTime = end_time - start_time)
    patient_attributes <- get_mon_attributes(basic_sim)
    #testdf <- getSingleAttribute('SocietalCosts',patient_attributes)
    
    testdf <- get_mon_attributes(basic_sim)
    
    
  #write.csv(testdf,'testdf100000.csv',row.names = F)
  #write.csv(EndDatDf,'Output_DFs/enddatdf_5000.csv',row.names = F)
  return(EndDatDf)
  })
return(results)
}

# Parallel Settings ----

patients = 500
runs=30
run.out <- RunSimPar(n.patients=patients,n.runs=runs)
#run.out[1,]

runOutDF <- as.data.frame(run.out)
runOutDF <- t(runOutDF)
runOutDF <- as.data.frame(runOutDF)


allAttributes = c("Index","Name","Alive", "PatientAge", "HealthyLifeExpectancy", "TumorStartAge", "LastAge", "TotalCosts", 
                  "MedicalCosts", "SocietalCosts","PAIDCosts", "TotalUtility", "CurrentUtility", "TumorSize", "WorstCancer", 
                  "WorstStage", "CancerStageT", "IsCured", "StagingProb", "TumorGrowthRate", "max_size", 
                  "RegressionSize", "StagnateSize", "ClinSizeCheck", "Invite", "Screened", "BIRADS", 
                  "Referral", "NextStep", "EntryAge", "TimeAtHome", "ClinicalCheck", "TimeAtDiagnostics", 
                  "TotalTimeAtDiagnostics", "DiagnosedThrough", "DiagnosticVisits", "TimeAtHospital", "TotalTimeAtHospital", 
                  "OriginalLifeExp", "HospitalVisits", "TimeAtScreening", "TotalTimeAtScreening", "FirstStage")

# newDF <- data.frame(matrix(rep('hi',(pats+1)*runs*length(allAttributes)),ncol=length(allAttributes)))
# colnames(newDF) <- allAttributes
# for (att in seq(length(allAttributes))){
#   print(att)
#   currentList = c()
#   for (j in seq(nrow(runOutDF))){
#     for( i in seq(length(runOutDF[[1,1]]))){
#         currentList <- append(currentList,runOutDF[,att][[j]][i])
#     }
#   }
#   newDF[,att] <- currentList
# }

newDF <- data.frame(matrix(ncol = length(allAttributes), nrow = length(unlist(runOutDF[,1]))))
colnames(newDF) <- allAttributes
for (att in seq(length(allAttributes))){
  print(allAttributes[att])
  newDF[,att] <- unlist(runOutDF[,att])
}

# newDF2 = newDF[!newDF[,1] == "test",]
# newDF2 = newDF %>%
#   filter(newDF[,1]!='test')
newDF2 = newDF[-which(newDF[,1]=='test'),]
newDF3 <- newDF2 %>%
  mutate(Index = row_number())

write.csv(newDF3,paste0('Output_DFs/20230315_',StratName,'_',runs*patients,'EndDatDf.csv'),row.names = F)

EndDatDf <- newDF3


EndDatDf2 = EndDatDf[-1,]

# Reporting ----
library(ggplot2)
library(gridExtra)
hour = 1;
minute = hour/60;
second=minute/60;
day=hour*24;
year = day*365; 
month = year/12

incpercentages <- c()
PercThroughScreens <- c()

df <- EndDatDf

numsdf = df[-1,]
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
meanAge = round(mean(round(as.numeric(numsdf$PatientAge)/year)),2)
medianAge = median(round(as.numeric(numsdf$PatientAge)/year))
modeAge = getmode(round(as.numeric(numsdf$PatientAge)/year))
patients_total <- nrow(df)
cancers_found <- table(as.numeric(df$WorstStage))
inc_perc <-   round((cancers_found[2]+cancers_found[3] + cancers_found[4]+cancers_found[5])/ patients_total * 100,2)
subdf <- df %>%
  dplyr::select(WorstStage,DiagnosedThrough) %>%
  na.omit()


subdf2 <- table(subdf)
subdf3 <- as.data.frame(subdf2)
WorstStage = rep(c(0.5,1,2,3,4),3)
#if no screening, this should be different:

ThroughClin = subdf2[,2][1:5]
ThroughScreen = subdf2[,1][1:5]

ThroughTotal = ThroughClin + ThroughScreen
throughlist = append(ThroughTotal,ThroughClin)
throughlist = append(throughlist,ThroughScreen)
percThroughScreen = round(sum(ThroughScreen) / (sum(ThroughScreen)+ sum(ThroughClin)) *100,2)

incpercentages<-append(incpercentages,inc_perc)
PercThroughScreens<-append(PercThroughScreens,percThroughScreen)

# x<- incpercentages
# #this doesn't work if only one simulation
# mean_x <- mean(x)
# median_x <- median(x)
# mode_x <- density(x)$x[which.max(density(x)$y)]
# ci_x <- t.test(x)$conf.int
# ggplot(data.frame(x), aes(x)) +
#   geom_histogram(bins = 10, color = 'black', fill = 'lightblue') +
#   geom_vline(xintercept = mean_x, color = 'red', linetype = 'dashed') +
#   geom_vline(xintercept = median_x, color = 'green', linetype = 'dashed') +
#   geom_vline(xintercept = mode_x, color = 'purple', linetype = 'dashed') +
#   geom_segment(aes(x = ci_x[1], y = 0, xend = ci_x[2], yend = 0), color = 'blue', size = 1.5) +
#   geom_point(aes(x = mean_x, y = 0), color = 'red', size = 2.5) +
#   geom_point(aes(x = median_x, y = 0), color = 'green', size = 2.5) +
#   geom_point(aes(x = mode_x, y = 0), color = 'purple', size = 2.5) +
#   annotate('text', x = ci_x[1], y = 3, label = round(ci_x[1], 2)) +
#   annotate('text', x = ci_x[2], y = 3, label = round(ci_x[2], 2)) +
#   labs(title = 'Histogram of x with 95% CI and Summary Statistics of incidence percentage',
#        x = 'x', y = 'Frequency')

# x<- PercThroughScreens 
# #this doesn't work if only one simulation
# mean_x <- mean(x)
# median_x <- median(x)
# mode_x <- density(x)$x[which.max(density(x)$y)]
# ci_x <- t.test(x)$conf.int
# ggplot(data.frame(x), aes(x)) +
#   geom_histogram(bins = 10, color = 'black', fill = 'lightblue') +
#   geom_vline(xintercept = mean_x, color = 'red', linetype = 'dashed') +
#   geom_vline(xintercept = median_x, color = 'green', linetype = 'dashed') +
#   geom_vline(xintercept = mode_x, color = 'purple', linetype = 'dashed') +
#   geom_segment(aes(x = ci_x[1], y = 0, xend = ci_x[2], yend = 0), color = 'blue', size = 1.5) +
#   geom_point(aes(x = mean_x, y = 0), color = 'red', size = 2.5) +
#   geom_point(aes(x = median_x, y = 0), color = 'green', size = 2.5) +
#   geom_point(aes(x = mode_x, y = 0), color = 'purple', size = 2.5) +
#   annotate('text', x = ci_x[1], y = 3, label = round(ci_x[1], 2)) +
#   annotate('text', x = ci_x[2], y = 3, label = round(ci_x[2], 2)) +
#   labs(title = 'Histogram of x with 95% CI and Summary Statistics of percentage through screens',
#        x = 'x', y = 'Frequency')

subdf <- df %>%
  dplyr::select('WorstStage','DiagnosedThrough','FirstStage') %>%
  mutate(WorstStage = ifelse(is.na(FirstStage),WorstStage,ifelse(WorstStage>FirstStage,WorstStage,FirstStage))) %>%
  select('WorstStage','DiagnosedThrough') %>%
  na.omit()

subdf2 <- table(subdf)
subdf3 <- as.data.frame(subdf2)
WorstStage = rep(c(0.5,1,2,3,4),3)
#if no screening, this should be different:

ThroughClin = subdf2[,2][1:5]
ThroughScreen = subdf2[,1][1:5]

ThroughTotal = ThroughClin + ThroughScreen
throughlist = append(ThroughTotal,ThroughClin)
throughlist = append(throughlist,ThroughScreen)
through = c('Total','Total','Total','Total','Total',
            'Clinical','Clinical','Clinical','Clinical','Clinical',
            'Screening','Screening','Screening','Screening','Screening')

data3 = data.frame(WorstStage,throughlist,through)

data3$WorstStage <- factor(data3$WorstStage,levels = c("4", "3", "2", "1", "0.5"))

# Create a vector of colors for each level of WorstStage
colors <- c("#cbe9a6", "#f4b26f", "#e576b4", "#1a7c98", "#29bdeb")

# Calculate the total count for each level of through
total_count <- aggregate(data3$throughlist, by = list(data3$through), sum)

# Calculate the percentage for each combination of WorstStage and through
data3$percent <- data3$throughlist / total_count[match(data3$through, total_count$Group.1), "x"] * 100


plot.new()
## Stadiumverdeling plot 2 ----
p2<-ggplot(data3,aes(fill=WorstStage,y=throughlist,x=through))+
  geom_bar(position='fill',stat='identity')+
  labs(title='Stadiumverdeling')+
  scale_fill_manual(values=colors) + 
  geom_text(aes(label = paste0(round(percent), "%")), position = position_fill(vjust = 0.5),size=2)


plot2 <- recordPlot()
plot.new()


subdf4 <- subdf3 %>%
  filter(DiagnosedThrough==0)

# filter subdf3 by DiagnosedThrough == 0 and drop the DiagnosedThrough column
subdf_filtered <- subset(subdf3, DiagnosedThrough == 0, select = -c(DiagnosedThrough))

# calculate the total  count of Freq values for each WorstStage
freq_total <- aggregate(Freq ~ WorstStage, subdf3, sum)$Freq

# divide the Freq column by the total count of Freq values for each WorstStage
subdf_filtered$Freq <- subdf_filtered$Freq / freq_total
# subdf_filtered <- subset(subdf_filtered, row_number() != 6)

# print the new DataFrame
#print(subdf_filtered)



## Percentage plot 1 ----
p1<-ggplot(subdf_filtered,aes(x=WorstStage,y=Freq))+
  geom_bar(stat='identity')+
  labs(title='Percentage tumoren ontdekt via het bevolkingsonderzoek')+
  scale_fill_brewer(palette="Spectral")+
  geom_text(aes(label = paste0(round(Freq*100, 1), "%")),
            position = position_stack(vjust = 0.5),size=2)
plot1 <- recordPlot()
plot.new()

subdf5 <- df %>%
  dplyr::select('WorstStage','TumorStartAge') %>%
  na.omit() %>%
  filter(TumorStartAge<100*year) %>%
  filter(WorstStage>0) %>%
  mutate(TumorStartAge = as.numeric(TumorStartAge)) %>%
  mutate(ages = cut(TumorStartAge,breaks=c(0,50*year,74*year,Inf))) %>%
  dplyr::select(ages,WorstStage)

subdf6 = rbind(as.data.frame(table(subdf5)))

# Calculate the total count for each level of through
total_count2 <- aggregate(subdf6$Freq, by = list(subdf6$WorstStage), sum)

# Calculate the percentage for each combination of WorstStage and through
subdf6$percent <- subdf6$Freq / total_count2[match(subdf6$WorstStage, total_count2$Group.1), "x"] * 100

## leeftijd diagnose plot 3 ----
p3<- ggplot(subdf6,aes(fill=ages,y=WorstStage,x=Freq))+
  geom_bar(position='fill',stat='identity')+
  coord_flip()+
  scale_fill_discrete(labels=c('<50', '50-74', '>74'))+
  scale_y_discrete(labels=c("0.5" = "DCIS", "1" = "Stage I", "2" = "StageII", "3" = "Stage III", "4" = "Stage IV"))+
  labs(title='Leeftijd bij diagnose voor vrouwen')+
  geom_text(aes(label = paste0(round(percent),'%')), position = position_fill(vjust = 0.5),size=2)
plot3 <- recordPlot()
plot.new()

rep(1,min(10,25))

## Survival curve plot 4 ----
subdf7 <- df %>%
  dplyr::select(TumorStartAge, WorstStage, PatientAge,OriginalLifeExp)%>%
  filter(TumorStartAge<100*year) %>%
  filter(WorstStage>0) %>%
  mutate(TumorStartAge = as.numeric(TumorStartAge)) %>%
  mutate(PatientAge = as.numeric(PatientAge)) %>%
  mutate(OriginalLifeExp = as.numeric(OriginalLifeExp)) %>%
  mutate(Delta = (PatientAge - TumorStartAge) / year) %>%
  mutate(Delta = cut(Delta,breaks=c(0,1,2,3,4,5,6,7,8,9,10,Inf))) %>%
  dplyr::select(WorstStage,Delta) %>%
  table() %>%
  as.data.frame(stringsAsFactors=F)
# create a new data frame with cumulative frequencies for each stage

# append(rep(1,min(10,8)),rep(0,(10-length(rep(1,min(10,8))))))

subdf7 <- df%>%
  dplyr::select(WorstStage,PatientAge, TumorStartAge,OriginalLifeExp)%>%
  filter(TumorStartAge<100*year) %>%
  filter(WorstStage>0) %>%
  mutate(WorstStage = as.numeric(WorstStage)) %>%
  mutate(TumorStartAge = as.numeric(TumorStartAge)/year) %>%
  mutate(PatientAge = as.numeric(PatientAge)/year) %>%
  mutate(OriginalLifeExp = as.numeric(OriginalLifeExp)/year) %>%
  mutate(Survival = PatientAge - TumorStartAge) %>%
  mutate(ShouldveSurvived = OriginalLifeExp - TumorStartAge)%>%
  mutate(One = ifelse(round(Survival)>=1,1,
                      ifelse(round(ShouldveSurvived)>=1,0,NA))) %>%
  mutate(Two = ifelse(round(Survival)>=2,1,
                      ifelse(round(ShouldveSurvived)>=2,0,NA))) %>%
  mutate(Three = ifelse(round(Survival)>=3,1,
                        ifelse(round(ShouldveSurvived)>=3,0,NA))) %>%
  mutate(Four = ifelse(round(Survival)>=4,1,
                       ifelse(round(ShouldveSurvived)>=4,0,NA))) %>%
  mutate(Five = ifelse(round(Survival)>=5,1,
                       ifelse(round(ShouldveSurvived)>=5,0,NA))) %>%
  mutate(Six = ifelse(round(Survival)>=6,1,
                      ifelse(round(ShouldveSurvived)>=6,0,NA))) %>%
  mutate(Seven = ifelse(round(Survival)>=7,1,
                        ifelse(round(ShouldveSurvived)>=7,0,NA))) %>%
  mutate(Eight = ifelse(round(Survival)>=8,1,
                        ifelse(round(ShouldveSurvived)>=8,0,NA))) %>%
  mutate(Nine = ifelse(round(Survival)>=9,1,
                       ifelse(round(ShouldveSurvived)>=9,0,NA))) %>%
  mutate(Ten = ifelse(round(Survival)>=10,1,
                      ifelse(round(ShouldveSurvived)>=10,0,NA)))

subdf7 <- subdf7 %>%
  group_by(WorstStage) %>%
  summarise(
    AAZero = c(1,1,1,1,1),
    AOne = mean(One, na.rm = TRUE),
    BTwo = mean(Two, na.rm = TRUE),
    CThree = mean(Three, na.rm = TRUE),
    DFour = mean(Four, na.rm = TRUE),
    EFive = mean(Five, na.rm = TRUE),
    FSix = mean(Six, na.rm = TRUE),
    GSeven = mean(Seven, na.rm = TRUE),
    HEight = mean(Eight, na.rm = TRUE),
    INine = mean(Nine, na.rm = TRUE),
    JTen = mean(Ten, na.rm = TRUE)
  )


subdf7 
# convert from wide to long format
df_long <- tidyr::pivot_longer(subdf7, cols = c(3:11), names_to = "Column", values_to = "Value")

# create the line plot
p4<-ggplot(data = df_long, aes(x = Column, y = Value, group = WorstStage, color = factor(WorstStage))) +
  geom_line() +
  scale_color_manual(values = c("black", "blue", "red", "green", "purple")) +
  labs(x = "Time (years)", y = "Relative survival", color = "Stage",title='Relatieve overleving voor pati?nten gediagnosticeerd met invasieve borstkanker naar stadium bij diagnose') 
plot4 <- recordPlot()
plot.new()


# More plots for finetuning

## Final age distribution ----
p5 <- ggplot(df, aes(x = as.numeric(PatientAge) / year)) +
  geom_histogram(binwidth = 1, color = "black", fill = "lightblue") +
  labs(title = "All patients final age distribution",
       x = "Age (years)",
       y = "Frequency")


## Histogram cancer start age ----
subdf8<- df %>%
  filter(TumorStartAge<1000*year)
p6 <- ggplot(subdf8, aes(x = as.numeric(TumorStartAge) / year)) +
  geom_histogram(bins = 100, color = "black", fill = "lightblue") +
  xlim(0, 100) +
  labs(title = "Histogram of cancer start age, rounded to year",
       x = "Age (years)",
       y = "Frequency")

## actual stage distribtuion ----
p7 <- ggplot(data3,aes(fill=WorstStage,y=throughlist,x=through))+
  geom_bar(position='stack',stat='identity')+
  scale_fill_manual(values=colors)

## Tumor growth rates hist ----
p8 <- ggplot(df,aes(x=as.numeric(TumorGrowthRate)))+
  geom_histogram(binwidth = 0.25 , color='black', fill='lightblue')+
  labs(title='Tumor growth rates',
       x='growth rate',
       y= 'frequency')+
  scale_x_log10()



patients_total <- nrow(df)
cancers_found <- table(as.numeric(df$WorstStage))
inc_perc <-   round((cancers_found[2]+cancers_found[3] + cancers_found[4]+cancers_found[5])/ patients_total * 100,2)
percThroughScreen = round(sum(ThroughScreen) / (sum(ThroughScreen)+ sum(ThroughClin)) *100,2)

numsdf = df[-1,]
meanAge = round(mean(round(as.numeric(numsdf$PatientAge)/year)),2)
medianAge = median(round(as.numeric(numsdf$PatientAge)/year))
modeAge = getmode(round(as.numeric(numsdf$PatientAge)/year))

#finetune stagefromsize diagram plot
## plot stage size hist ----

SizeTestdf<- df %>% 
  dplyr::select(DiagnosedThrough,WorstCancer,WorstStage) %>% 
  filter(WorstStage == 1) %>%
  mutate(WorstCancer = as.numeric(WorstCancer))

ggplot(SizeTestdf, aes(x = WorstCancer, fill = factor(DiagnosedThrough))) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 50) +
  scale_fill_discrete(name = "Diagnosed Through") +
  labs(title = "Histogram of Worst Cancer by Diagnosis Method", 
       x = "Worst Cancer", y = "Count")+
  xlim(c(0,130))


grid.arrange(p5,p6,p7,p8,nrow=2)
grid.arrange(p1,p2,p3,p4,nrow=2)

CostUtilDf <- df %>%
  select(PatientAge,TotalUtility,TotalCosts,MedicalCosts,SocietalCosts,PAIDCosts) %>%
  mutate_all(function(x) as.numeric(x)) %>%
  mutate_all(function(x) ifelse(is.na(x),0,x)) %>%
  mutate(TotalCosts = SocietalCosts + MedicalCosts+PAIDCosts)

write.csv(CostUtilDf,paste0('Output_DFs/20230315_',StratName,'_',runs*patients,'CostUtilDf.csv'))

meanQALY <- mean(CostUtilDf$TotalUtility,na.rm = T) /year
meanTotalCosts <- mean(CostUtilDf$TotalCosts,na.rm = T)
meanMedicalCosts <- mean(CostUtilDf$MedicalCosts,na.rm = T)
meanSocietalCosts <- mean(CostUtilDf$SocietalCosts,na.rm = T)
meanPaidCosts <- mean(CostUtilDf$PAIDCosts,na.rm=T)
# maybe add overhead? 


## KPIs ----

# end_time <- Sys.time()
# RunTime <- end_time - start_time
# print(RunTime)

cat(paste0(strrep("_", 70),'\n',
           'OUTPUT KPIs', '\n',
           strrep("_", 70),'\n',
           'Patients generated: ',patients_total,'\n',
           'incidence percentage: ',inc_perc,'%','\n',
           'Percentage found through screening: ',percThroughScreen,'%','\n',
           'Patient Ages: mean:', meanAge,' mode: ',modeAge, ' median: ',medianAge, '\n',
           'Mean Qulity Adjusted Life Expectancy: ', meanQALY,'\n',
           'Mean Total Costs: ', meanTotalCosts,'\n',
           'Mean Medical Costs: ', meanMedicalCosts,'\n',
           'Mean Societal Costs: ', meanSocietalCosts,'\n',
           'Mean PAID Costs: ', meanPaidCosts,'\n'))







