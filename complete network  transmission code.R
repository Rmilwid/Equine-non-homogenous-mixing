library("EpiModel")
library(networkDynamicData)
library(igraph)
library(intergraph)
library(xlsx)
setwd("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/complete network")

nw1<-make_full_graph(24,directed = F,loops=F)

new1<-asNetwork(nw1,bipartite=F)


#this function extends the collected data to span a chosen time range
#dayn= network for day n
#total= number of times you want it repeated.e.g. to get two weeks worth of data, list 7 days and total=2
extend.data.fun<-function(day1,day2=NULL,day3=NULL,day4=NULL,day5=NULL,day6=NULL,day7=NULL,total){
  a<-list(day1)
  if(!is.null(c(day2))) a<-list(day1,day2)
  if(!is.null(day3)) a<-list(day1,day2,day3)
  if(!is.null(day4)) a <-list(day1,day2,day3,day4)
  if(!is.null(day5)) a<-list(day1,day2,day3,day4,day5)
  if(!is.null(day6)) a<-list(day1,day2,day3,day4,day5,day6)
  if(!is.null(day7)) a<-list(day1,day2,day3,day4,day5,day6,day7)
  newcombined<-rep(a,total)
  return(newcombined)
}
newcombined<-extend.data.fun(day1 = new1,day2=new1,day3=new1,day4 = new1,day5=new1,day6 = new1,day7=new1, total=3)





newcombinedDyn<-networkDynamic(network.list = newcombined)

newnetwork<-newcombinedDyn

nw<-delete.vertex.attribute(newnetwork,"status.active")

init.net.mod <- function(x, param, init, control, s) {
  
  # Master Data List
  dat <- list()
  dat$param <- param
  dat$init <- init
  dat$control <- control
  
  dat$attr <- list()
  dat$stats <- list()
  dat$temp <- list()
  
  # Network Parameters
  dat$nw <- x
  dat$param$modes <- 1
  
  # Initialization
  
  ## Infection Status and Time Modules
  n <- network.size(dat$nw)
  dat$attr$status <- rep("s", n)
  dat$attr$status[sample(1:n, init$i.num)] <- "i"
  
  dat$attr$active <- rep(1, n)
  dat$attr$entrTime <- rep(1, n)
  dat$attr$exitTime <- rep(NA, n)
  
  dat$attr$Quar <- rep(0, n)
  dat$attr$vac <- rep(0, n)
  dat$attr$vac[sample(1:n, init$v.num)] <- 1
  
  dat$attr$infTime <- rep(NA, n)
  dat$attr$infTime[dat$attr$status == "i"] <- 1
  dat$attr$infTime[dat$attr$vac == 1] <- 1
  ## Get initial prevalence
  dat <- get_prev.net(dat, at = 1)
  
  return(dat)
}


#############################3



vaccinate<-function(dat,at){
  active <- dat$attr$active
  status <- dat$attr$status
  vaccinated <- dat$attr$vac
  
  
  
  
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  
  #set attributes in attr sublist
  vaccine.rate <- dat$param$vaccine.rate
  
  nVac<-0
  
  
  #decide which horses are eligible for vaccination and vaccinate them at rate vaccine.rate
  
  
  
  
  idsEligVac<-which(active==1 & status =="s" & vaccinated == 0) #only susceptible horses are eligible for vaccination
  nEligVac<-length(idsEligVac) # number of eligible horses to be vaccinated
  if (nEligVac > 0) {
    vecVac <- which(rbinom(nEligVac, 1, vaccine.rate) == 1) #vaccinate a horse with prob vaccine rate. use a            binomial distribution to decide who to vaccinate 
    
    if (length(vecVac) > 0) {
      
      idsVac <- idsEligVac[vecVac] #id of the vaccinated horses
      nVac <- length(idsVac)   #calculate number of vaccinated horses
      
      vaccinated[idsVac] <- 1   #assign vaccinated horses a status of 1
      
    }
    
  }
  
  
  if (at == 2) {
    
    dat$epi$nVac<-c(0,nVac)
  } else {
    
    dat$epi$nVac[at]<-nVac
  }
  
  dat$nw <- nw
  # dat$attr$active<-active 
  dat$attr$status<- status
  dat$attr$vac<-vaccinated
  return(dat)
}



infect <- function(dat, at) {
  
  active <- dat$attr$active
  status <- dat$attr$status
  vaccinated <- dat$attr$vac
  quarantine <- dat$attr$Quar
  infTime<-dat$attr$infTime
  
  nw <- dat$nw
  
  idsSus <- which(active == 1 & status == "s")
  idsInf <- which(active == 1 & status == "i")
  
  
  
  
  
  ###############decide which horses will become infected based on vaccination status
  
  nActive <- sum(active == 1)
  
  nElig <- length(idsInf)
  
  nInf <- 0
  
  if (nElig > 0 && nElig < nActive) {
    
    del <- discord_edgelist(dat, at)
    if (!(is.null(del))) {
      
      del$vaccinatedSus<-dat$attr$vac[del$sus]
      del$vaccinatedInf<-dat$attr$vac[del$inf]
      
      nVac<-length(which((del$vaccinatedSus == 1| del$vaccinatedInf == 1 ) &  del$quarantine ==0))  
      nNoVac<-nrow(del)-nVac
      del$quarantine<-dat$attr$Quar[del$inf]
      del$transProbVac <-dat$param$inf.prob.vac
      del$transProb <- dat$param$inf.prob
      del$actRate <- dat$param$act.rate
      del$finalProbVac <- 1 - (1 - del$transProbVac)^del$actRate
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      
      delnovac <- del[ which(del$vaccinatedSus == 0 & del$vaccinatedInf == 0 & del$quarantine ==0), ]
      delvac <- del[which((del$vaccinatedSus == 1| del$vaccinatedInf == 1) & del$quarantine ==0 ), ]
      transmit <- rbinom(nrow(delnovac), 1, del$finalProb)
      transmitVac <- rbinom(nrow(delvac), 1, del$finalProbVac)
      
      delnovac <- delnovac[which(transmit == 1 ), ]
      delvac <- delvac[which(transmitVac == 1 ), ]
      
      idsNewInf <- unique(delnovac$sus)
      idsNewInfVac <- unique(delvac$sus)
      nInf <- sum(length(idsNewInf),(length(idsNewInfVac)))
      
      if (nInf > 0) {
        status[idsNewInf] <- "e"
        infTime[idsNewInf] <- at
        status[idsNewInfVac] <- "e"
        infTime[idsNewInfVac] <- at
        
        
      }
    }
  }
  
  
  if (at == 2) {
    dat$epi$se.flow <- c(0, nInf)
    
  } else {
    dat$epi$se.flow[at] <- nInf
    
  }
  
  dat$nw <- nw
  dat$attr$active<-active 
  dat$attr$status<- status
  dat$attr$vac<-vaccinated
  dat$attr$Quar <-quarantine
  
  
  return(dat)
}

progress <- function(dat, at) {
  active <- dat$attr$active
  status <- dat$attr$status
  vaccinated <-dat$attr$vac
  quarantine <- dat$attr$Quar
  evi.rate <- dat$param$evi.rate #transmssion rate for E->I for vaccinated horses
  ei.rate <- dat$param$ei.rate #transmssion rate for E->I for unvaccinated horses
  ivr.rate <- dat$param$ivr.rate #transmssion rate for I->R for vaccinated horses
  ir.rate <- dat$param$ir.rate #transmssion rate for I->R for unvaccinated horses
  qr.rate <- dat$param$qr.rate #transmssion rate for Q->R for unvaccinated horses (that's the only possible route for this model)
  quarantine.rate <- dat$param$quarantine.rate
  nQuar<-0
  
  ## E to I progression
  nInf <- 0
  nInfVac <- 0
  idsEligInf <- which(active == 1 & status == "e" & vaccinated == 0)
  idsEligInfVac <- which(active == 1 & status == "e" & vaccinated == 1)
  
  nEligInf <- length(idsEligInf)
  nEligInfVac <- length(idsEligInfVac)
  
  if (nEligInf > 0) {
    vecInf <- which(rbinom(nEligInf, 1, ei.rate) == 1)
    if (length(vecInf) > 0) {
      idsInf <- idsEligInf[vecInf]
      nInf <- length(idsInf)
      status[idsInf] <- "i"
    }
  }
  
  if (nEligInfVac > 0) {
    vecInfVac <- which(rbinom(nEligInfVac, 1, evi.rate) == 1)
    if (length(vecInfVac) > 0) {
      idsInfVac <- idsEligInfVac[vecInfVac]
      nInfVac <- length(idsInfVac)
      status[idsInfVac] <- "i"
    }
  }
  ################decide which horses are eligible for quarantine and quaratine them at rate quarrantine.rate
  nQuar<-0
  
  idsEligQuar<-which(active==1 & status =="i" & vaccinated ==0 & quarantine == 0) #only infectious horses are eligible for quarantine
  nEligQuar<-length(idsEligQuar)
  
  if (nEligQuar > 0) {
    vecQuar <- which(rbinom(nEligQuar, 1, quarantine.rate) == 1) #quarantine a horse with prob quarantine rate. use a binomial distribution to decide who to quarantine
    
    if (length(vecQuar) > 0) {
      
      idsQuar <- idsEligQuar[vecQuar] #id of the quarantined horses
      nQuar <- length(idsQuar)   #calculate number of quarantined horses
      quarantine[idsQuar] <- 1    #assign quarantined horses a status of 1
    }
    
  }
  
  
  
  ## I to R progression
  nRec <- 0
  nRecQuar<-0
  nRecVac<-0
  idsEligRecQuar <- which(active == 1 & status == "i" & quarantine ==1 & vaccinated ==0)
  idsEligRecVac <- which(active == 1 & status == "i" & vaccinated ==1 & quarantine ==0)
  idsEligRec <- which(active == 1 & status == "i" & quarantine ==0 & vaccinated == 0)
  
  
  nEligRec <- length(idsEligRec)
  nEligRecVac <- length(idsEligRecVac)
  nEligRecQuar <- length(idsEligRecQuar)
  
  
  if (nEligRec > 0) {
    
    vecRec <- which(rbinom(nEligRec, 1, ir.rate) == 1 )
    
    if (length(vecRec) > 0) {
      idsRec <- idsEligRec[vecRec]
      nRec <- length(idsRec)
      status[idsRec] <- "r"
      quarantine[idsRec]<-0
      
    }
  }
  if (nEligRecQuar > 0) {
    
    vecRecQuar <- which(rbinom(nEligRecQuar, 1, qr.rate) == 1)
    
    if (length(vecRecQuar) > 0) {
      
      idsRecQuar <- idsEligRecQuar[vecRecQuar]
      nRecQuar <- length(idsRecQuar)
      status[idsRecQuar] <- "r"
      quarantine[idsRecQuar]<-0
      
    }
  }
  
  if (nEligRecVac > 0) {
    
    vecRecVac <- which(rbinom(nEligRecVac, 1, ivr.rate) == 1)
    
    if (length(vecRecVac) > 0) {
      idsRecVac <- idsEligRecVac[vecRecVac]
      nRecVac <- length(idsRecVac)
      status[idsRecVac] <- "r"
      quarantine[idsRecVac]<-0
    }
  }
  
  
  dat$attr$status <- status
  dat$attr$Quar <- quarantine
  dat$attr$vac <- vaccinated
  
  if (at == 2) {
    dat$epi$ei.flow <- c(0, sum(nInf,nInfVac))
    dat$epi$ir.flow <- c(0, sum(nRec,nRecQuar,nRecVac))
    dat$epi$e.num <- c(0, sum(active == 1 & status == "e"))
    dat$epi$r.num <- c(0, sum(active == 1 & status == "r"))
    dat$epi$nQuar <- c(0, nQuar)
  } else {
    dat$epi$ei.flow[at] <- sum(nInf,nInfVac)
    dat$epi$ir.flow[at] <- sum(nRec,nRecQuar,nRecVac)
    dat$epi$e.num[at] <- sum(active == 1 & status == "e")
    dat$epi$r.num[at] <- sum(active == 1 & status == "r")
    dat$epi$nQuar[at] <- nQuar
    
  }
  
  return(dat)
}

param <- param.net(vaccine.rate=0, #vaccination rate
                   # vaccine.efficacy=.5, #vaccine efficacy
                   inf.prob.vac=.5,  #infection prob. for vaccinated horses
                   inf.prob = 1,  #infection prob. for unvaccinated horses
                   act.rate = 1,   #number of acts
                   evi.rate=1/2.52,     #transition rate from exposed to infectious compartments for vaccinated individuals
                   ei.rate = 1/1.75, #transition rate from exposed to infectious compartments for unvacciated individuals
                   ivr.rate=1/2.5,     #transition rate from infectious to recovered compartments for vacciated individuals
                   qr.rate=1/7,       #transition rate from infectious to recovered compartments for quarantined individuals
                   ir.rate = 1/4.8,   #transition rate from infectious to recovered compartments for unvacciated individuals
                   quarantine.rate=0)

init <- init.net(i.num = 1, v.num =.25*24 ,status.rand = FALSE)

control <- control.net(type = "SI", nsteps = 42, nsims = 10000,
                       initialize.FUN = init.net.mod,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       # vaccinate.FUN=vaccinate,
                       recovery.FUN = NULL, 
                       skip.check = TRUE,
                       depend = TRUE,
                       verbose.int = 1,save.nwstats = F,save.network = F,
                       module.order = c(
                         # "vaccinate.FUN",
                         "infection.FUN",
                         "progress.FUN", 
                         "get_prev.FUN"))

sim <- netsim(nw, param, init, control)


##############
##############3table of values

s1<-as.data.frame(sim) #mean across sims
##export mean values to csv
write.csv(s1,file="quar_0_vac_25_mean_values.csv")

