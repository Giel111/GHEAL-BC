# Intro ---
## Outline including all branches,
##    No parameters yet
#setwd('C:/Users/gielv/OneDrive/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel.')
# Libraries ----
#workdir
setwd('C:/Users/gielv/OneDrive/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel') #<-home
#setwd('Z:/UT/IEM/001 AFSTUDEREN/RProjects/ThesisModel')                      #<-work

#library(doParallel)
#registerDoParallel(cores=detectCores(all.tests = FALSE, logical = TRUE)-1)


KPIdf <- data.frame(kkstrat = numeric(),
                 ICER = numeric(),
                 AvgCosts = numeric(),
                 AvgEffects = numeric(),
                 IncPerc = numeric(),
                 SDCIS = numeric(),
                 SI = numeric(),
                 SII = numeric(),
                 SIII = numeric(),
                 SIV = numeric())
StratDf <- read.csv('2KFactorialDesignOneZeroCols.csv')
getScreenAges <- function(strategy, StratDf) {
  if (strategy==1){
    return(c(200,210)*8760)
  }
  # Get the row for the given strategy
  row <- StratDf[StratDf$ï..ExpNo == strategy,]
  # Get the column names where the value is 1
  colNames <- names(row)[which(row == 1)]
  # Extract the screen ages from the column names
  screenAges <- as.numeric(gsub("X", "", colNames))
  return(screenAges*8760)
}



for (kkstrat in 1:nrow(StratDf)){
  
  
  # basicallly antoher loop for the entire thing
  # first define current strategy to use:
  # STRATEGY SELECTOR ----
  ScreenAges <- getScreenAges(kkstrat,StratDf)
  # here you can select the screening strategy for this simulation run.


  
    
  
  # rm(list=ls())
  gc()
  
  StratName = paste0('1_2KFactorialStrategy_',kkstrat)
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
  
  
  RunSimPar <- function(n.patients,n.runs,free.cores=0,screenages = ScreenAges){
    #set.seed(round(runif(1,1,1000)))
    #set.seed(111)
    #library(parallel)
    ScreenAges <- screenages
  
    cl <- makeCluster(10);
    registerDoSNOW(cl);
    clusterExport(cl, c("getSingleAttribute"));
    
    results <- pbsapply(cl=cl,X=1:n.runs,FUN=function(run){
      set.seed(run)
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
      n.patients=250;
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
      
      FP_Perc <- 0.005 #FalsePositive
      
      
      # ## Screen Start Params ----
      # Start_screen_age <- Start_screen_age *year
      # end_screen_age <- end_screen_age * year
      # screen_interval <- screen_interval*year
      # ScreenAges = seq(Start_screen_age,end_screen_age,by=screen_interval)
      # ScreenAges <- append(ScreenAges,1000*year)
      # if (manualScreenInput){
      #   # enter ages at which patient should be screened here:
      #   ScreenAges <- manualScreenAges * year
      # }
      # # ScreenAges = c(50,52,54,56,58,60,62,64,66,68,70,72,74)*year
      
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
  
  patients = 250
  runs=10
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
  
  write.csv(newDF3,paste0('2KFactorialOutputs/20230316_',StratName,'_',runs*patients,'EndDatDf.csv'),row.names = F)
  
  
  df<- newDF3
  # get only the KPIS i'm interested in:
  
  hour <- 1
  day=hour*24;
  year = day*365;
  

  
  df['TotalUtility'][is.na(df['TotalUtility'])] <- 786867.1
  df['TotalCosts'][is.na(df['TotalCosts'])] <- 0
  df['MedicalCosts'][is.na(df['MedicalCosts'])] <- 0
  df['SocietalCosts'][is.na(df['SocietalCosts'])] <- 0
  df['PAIDCosts'][is.na(df['PAIDCosts'])] <- 0
  df['TotalCosts'] = as.numeric(df$MedicalCosts) + as.numeric(df$SocietalCosts) + as.numeric(df$PAIDCosts)
  

  
  TotalCosts <- sum(as.numeric(df$TotalCosts))
  TotalEffects <- sum(as.numeric(df$TotalUtility)) /year
  TotalLY <- sum(as.numeric(df$PatientAge)) / year
  TotalMedical <-sum(as.numeric(df$MedicalCosts))
  TotalSocietal <- sum(as.numeric(df$SocietalCosts))
  TotalPAID <- sum(as.numeric(df$PAIDCosts))
  
  
  testedPatients <- nrow(df)
  
  AvgCosts <- TotalCosts / testedPatients
  AvgEffects <- TotalEffects / testedPatients
  AvgLY <- TotalLY / testedPatients
  AvgMedical <- TotalMedical / testedPatients
  AvgSocietal <- TotalSocietal / testedPatients
  AvgPAID <- TotalPAID / testedPatients
  
  
  StageTable<-table(df$WorstStage)
  IncPerc <- (sum(StageTable)- StageTable[1]) / testedPatients
  Cancers <- sum(StageTable)- StageTable[1]
  SDCIS <- StageTable[2]/Cancers 
  SI <-  StageTable[3]/Cancers 
  SII <-  StageTable[4]/Cancers 
  SIII <- StageTable[5]/Cancers 
  SIV <-  StageTable[6]/Cancers 
  
  
  
  BaselineCostsPP <- 17367.21
  BaselineEffectsPP <- 72.29
  BaselineLYPP <- 80.13
  BaselineMedical <- 6700.64
  BaselineSocietal <- 517.22
  BaselinePAID <- 10149.34
  
  ICER <- (BaselineCostsPP - AvgCosts) / (BaselineEffectsPP - AvgEffects)
  
  
  KPIList <-c(kkstrat,
               ICER,
               AvgCosts,
               AvgEffects,
               IncPerc,
               SDCIS,
               SI,
               SII,
               SIII,
               SIV)
  KPIdf <- rbind(KPIdf,KPIList)
  
  write.csv(KPIdf,'1_2kFactorialKPIs.csv')

}






