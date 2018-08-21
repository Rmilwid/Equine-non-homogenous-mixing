require(EpiModel)
#assume average population size and average degree from data:
#average population size= (20+28+26+24)/4=24.5=24
setwd("C:/Users/a/Documents/PhD/thesis/thesis doc/Thesis articles/comparing mutiple farms influenza/network with degree 4/")
initialize.net.A <- function(x, param, init, control, s) {
  
  if (control$start == 1) {
    # Master Data List --------------------------------------------------------
    dat <- list()
    dat$param <- param
    dat$init <- init
    dat$control <- control
    
    dat$attr <- list()
    dat$stats <- list()
    dat$temp <- list()
    
    
    # Network Simulation ------------------------------------------------------
    if (class(x$fit) == "network") {
      nw <- simulate(x$formation,
                     basis = x$fit,
                     coef = x$coef.form.crude,
                     constraints = x$constraints)
    } else {
      nw <- simulate(x$fit, control = control$set.control.ergm)
    }
    modes <- ifelse(nw %n% "bipartite", 2, 1)
    if (control$depend == TRUE) {
      if (class(x$fit) == "stergm") {
        nw <- network.collapse(nw, at = 1)
      }
      nw <- sim_nets(x, nw, nsteps = 1, control)
    }
    if (control$depend == FALSE) {
      nw <- sim_nets(x, nw, nsteps = control$nsteps, control)
    }
    nw <- activate.vertices(nw, onset = 1, terminus = Inf)
    
    
    # Network Parameters ------------------------------------------------------
    dat$nw <- nw
    dat$nwparam <- list(x[-which(names(x) == "fit")])
    dat$param$modes <- modes
    
    
    # Initialization ----------------------------------------------------------
    
    ## Infection Status and Time Modules
    dat <- init_status.net(dat)
    
    
    ## Initialize persistent IDs
    if (control$use.pids == TRUE) {
      dat$nw <- init_pids(dat$nw, dat$control$pid.prefix)
    }
    
    
    ## Pull network val to attr
    form <- get_nwparam(dat)$formation
    fterms <- get_formula_terms(form)
    dat <- copy_toall_attr(dat, at = 1, fterms)
    
    
    ## Store current proportions of attr
    dat$temp$t1.tab <- get_attr_prop(dat$nw, fterms)
    
    
    ## Get initial prevalence
    dat <- get_prev.net(dat, at = 1)
  } else {
    dat <- list()
    
    dat$nw <- x$network[[s]]
    dat$param <- x$param
    dat$control <- control
    dat$nwparam <- x$nwparam
    dat$epi <- sapply(x$epi, function(var) var[s])
    names(dat$epi) <- names(x$epi)
    dat$attr <- x$attr[[s]]
    dat$stats <- sapply(x$stats, function(var) var[[s]])
    dat$temp <- list()
  }
  
  return(dat)
}



init_status.net <- function(dat) {
  
  # Variables ---------------------------------------------------------------
  tea.status <- dat$control$tea.status
  i.num <- dat$init$i.num
  i.num.m2 <- dat$init$i.num.m2
  r.num <- dat$init$r.num
  r.num.m2 <- dat$init$r.num.m2
  v.num <- dat$init$v.num
  
  status.vector <- dat$init$status.vector
  num <- network.size(dat$nw)
  form <- get_nwparam(dat)$form
  statOnNw <- "status" %in% get_formula_terms(form)
  
  modes <- dat$param$modes
  if (modes == 1) {
    mode <- rep(1, num)
  } else {
    mode <- idmode(dat$nw)
  }
  
  type <- dat$control$type
  
  
  # Status ------------------------------------------------------------------
  
  ## Status passed on input network
  if (statOnNw == TRUE) {
    status <- get.vertex.attribute(dat$nw, "status")
  } else {
    if (!is.null(status.vector)) {
      status <- status.vector
    } else {
      status <- rep("s", num)
      status[sample(which(mode == 1), size = i.num)] <- "i"
      if (modes == 2) {
        status[sample(which(mode == 2), size = i.num.m2)] <- "i"
      }
      if (type == "SIR") {
        status[sample(which(mode == 1 & status == "s"), size = r.num)] <- "r"
        if (modes == 2) {
          status[sample(which(mode == 2 & status == "s"), size = r.num.m2)] <- "r"
        }
      }
    }
  }
  dat$attr$status <- status
  dat$attr$vac<-rep(0,num)
  dat$attr$vac[sample(1:num, init$v.num)] <- 1
  
  
  ## Save out other attr
  dat$attr$Quar<-rep(0,length(status))
  dat$attr$active <- rep(1, length(status))
  dat$attr$entrTime <- rep(1, length(status))
  dat$attr$exitTime <- rep(NA, length(status))
  if (tea.status == TRUE) {
    dat$nw <- activate.vertex.attribute(dat$nw,
                                        prefix = "testatus",
                                        value = status,
                                        onset = 1,
                                        terminus = Inf)
  }
  
  
  # Infection Time ----------------------------------------------------------
  ## Set up inf.time vector
  idsInf <- which(status == "i")
  infTime <- rep(NA, length(status))
  
  if (!is.null(dat$init$infTime.vector)) {
    infTime <- dat$init$infTime.vector
  } else {
    # If vital dynamics, infTime is a geometric draw over the duration of infection
    if (dat$param$vital == TRUE && dat$param$di.rate > 0) {
      if (dat$control$type == "SI") {
        infTime[idsInf] <- -rgeom(n = length(idsInf), prob = dat$param$di.rate) + 2
      } else {
        infTime[idsInf] <- -rgeom(n = length(idsInf),
                                  prob = dat$param$di.rate +
                                    (1 - dat$param$di.rate)*mean(dat$param$rec.rate)) + 2
      }
    } else {
      if (dat$control$type == "SI" || mean(dat$param$rec.rate) == 0) {
        # if no recovery, infTime a uniform draw over the number of sim time steps
        infTime[idsInf] <- ssample(1:(-dat$control$nsteps + 2),
                                   length(idsInf), replace = TRUE)
      } else {
        infTime[idsInf] <- -rgeom(n = length(idsInf), prob = mean(dat$param$rec.rate)) + 2
      }
    }
  }
  
  dat$attr$infTime <- infTime
  
  return(dat)
}

# 
# init_pids <- function(nw, prefixes=c("F", "M")) {
#   
#   if (is.null(nw$gal$vertex.pid)) {
#     if (nw$gal$bipartite == FALSE) {
#       nw <- initialize.pids(nw)
#     } else {
#       t0.pids <- c(paste0(prefixes[1], 1:length(modeids(nw, 1))),
#                    paste0(prefixes[2], 1:length(modeids(nw, 2))))
#       
#       nw <- set.network.attribute(nw, "vertex.pid", "vertex.names")
#       nw <- set.vertex.attribute(nw, "vertex.names", t0.pids)
#     }
#   }
#   
#   return(nw)
# }








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
                   quarantine.rate=1)

init <- init.net(i.num = 1, v.num =1*24 ,status.rand = FALSE)


control <- control.net(type = "SI", nsteps = 42, nsims = 10000,
                       # control <- control.net(type = "SI", nsteps = 20, nsims = 5,
                       initialize.FUN = initialize.net.A,
                       infection.FUN = infect,
                       progress.FUN = progress,
                       # vaccinate.FUN=vaccinate,
                       recovery.FUN = NULL, 
                       skip.check = TRUE,
                       depend = T, verbose.int = 1,save.nwstats = F,save.network = F,
                       module.order = c(
                         # "vaccinate.FUN",
                         "infection.FUN",
                         "progress.FUN", 
                         "get_prev.FUN"))


nw <- network.initialize(n = 24, directed = FALSE)
## ----netEst1, results = "hide"-------------------------------------------
#assume each horse comes in contact with his stall neighboours (2 horses), one horse in the cross ties, and one horsein the arena
#thereofre each horse has, on averaeg degree=4
#therefore n_edges=4*24/2
#               =48
est1 <- netest(nw, formation=~edges, target.stats=48,coef.diss=dissolution_coefs(~offset(edges), duration=1))




sim <- netsim(est1, param, init, control)



s1<-as.data.frame(sim) #mean across sims
##export mean values to csv
write.csv(s1,file="quar_100_vac_100_mean_values.csv")

