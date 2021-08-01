######################################
### TWO-STEP VAX LITTLE MODEL ########
######################################
# This file contains the parameters
# about the vaccines and vaccination
# campaigns in order to simulate
# different scenarios.

#####################################
######## General Notation ###########
#####################################
# About the units: .DAYS in days
#                  .RATE in /day
#                  .FRAC in [0, 1]
#                  .NUM  integer
# About the ages:
# 3 age classes: juvenile (J) <20yrs
#                adults (A) 20-60yrs
#                elderly (I)   >60yrs
# Parameters which depend on age are
# represented as vectors in order: J A I
source('./functions/parms_epi.R')
source('./functions/opt_vax_rate.R')

############################
#### VAX MANUFACTURE ######
############################

VAX.INITIAL.STORAGE.NUM = 6e6 # Number of READY vaccines on day one
VAX.PRODUCTION.RATE = 5e5 # Number of vaccines produced per day
MAX.VAC.RATE = 1e6 # Max number of vaccine applications per day

# above numbers are for whole country, correcting (roughly) by population
VAX.INITIAL.STORAGE.NUM = POP.ESTADO.REL.FRAC * VAX.INITIAL.STORAGE.NUM
VAX.PRODUCTION.RATE = POP.ESTADO.REL.FRAC * VAX.PRODUCTION.RATE
MAX.VAC.RATE = POP.ESTADO.REL.FRAC * MAX.VAC.RATE

# distributes current (first dose) vaccines among age groups
VAX.DISTR.RATE <- function(V, POP, tol=100) {
  if (POP[3] > tol)
      return(c(0, V - min(V, POP[3]), min(V, POP[3])))
  if (POP[2] > tol)
      return(c(V - min(V, POP[2]), min(V, POP[2]), 0))
  if (POP[1] > tol)
      return(c(min(V, POP[1]), 0, 0))
  return(c(0, 0, 0))
}

OBSERVED.EFFICACY.TO.PARAMETERS <- function(HOSP.PROP,MORT.PROP,ASYMP.PROP,SUSCEP.EFF,HOSP.EFF,MORT.EFF,SYMP.EFF){
  SUSCEP.EFF.PAR <- SUSCEP.EFF
  HOSP.EFF.PAR <- 1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)
  MORT.EFF.PAR <- 1 - (1-MORT.EFF)/(1-HOSP.EFF)
  SYMP.EFF.PAR <- 1 - ((1-SYMP.EFF)*(HOSP.PROP + (1 - HOSP.PROP)*(1-ASYMP.PROP))  - (1 - HOSP.EFF)*HOSP.PROP)/
                  ((1-SUSCEP.EFF)*(1 - (1-HOSP.EFF)/(1-SUSCEP.EFF)*HOSP.PROP)*(1-ASYMP.PROP))
  return(c(SUSCEP.EFF.PAR,HOSP.EFF.PAR,MORT.EFF.PAR,SYMP.EFF.PAR))
}
# Time window between first and second vaccines
VAX.WINDOW.DAYS = 21
#VAX.WINDOW.DAYS = 84 # 12 weeks

# Fraction of people that take the first but not the second dose \theta
SECOND.VAX.LOSS.FRAC = 0.1

# time until end of vaccination schedule, in days
MAX.TIME.DAYS = 300

# plot vaccination schedule
#plot_vac_schedule(OPT.VAX.RATE, VAX.INITIAL.STORAGE.NUM,
#                  VAX.PRODUCTION.RATE, MAX.VAC.RATE, VAX.WINDOW.DAYS,
#                  SECOND.VAX.LOSS.FRAC, MAX.TIME.DAYS)

#  Relative average contact rate and success between vaccinated susceptibles and infectious individuals \beta_v \beta_w
OBS.VAX2.EFFIC.BETA = 0.0
OBS.VAX2.EFFIC.CLIN = 0.5
OBS.VAX2.EFFIC.SEVERE = 0.83
OBS.VAX2.EFFIC.DEATH = 0.95 ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.8
OBS.VAX1.EFFIC.BETA   = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.BETA
OBS.VAX1.EFFIC.CLIN   = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.CLIN
OBS.VAX1.EFFIC.SEVERE = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.SEVERE
OBS.VAX1.EFFIC.DEATH  = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.DEATH

#transforming observed efficacy in efficacy parameters
VAX1.PARS <- OBSERVED.EFFICACY.TO.PARAMETERS(SEVERITY.FRAC,DEATH.FRAC,ASYMPTOMATIC.FRAC,
                                             OBS.VAX1.EFFIC.BETA,OBS.VAX1.EFFIC.SEVERE,
                                             OBS.VAX1.EFFIC.DEATH,OBS.VAX1.EFFIC.CLIN)
VAX2.PARS <- OBSERVED.EFFICACY.TO.PARAMETERS(SEVERITY.FRAC,DEATH.FRAC,ASYMPTOMATIC.FRAC,
                                             OBS.VAX2.EFFIC.BETA,OBS.VAX2.EFFIC.SEVERE,
                                             OBS.VAX2.EFFIC.DEATH,OBS.VAX2.EFFIC.CLIN)
VAX1.EFFIC.BETA   <- VAX1.PARS[1]
VAX1.EFFIC.SEVERE <- VAX1.PARS[2]
VAX1.EFFIC.DEATH  <- VAX1.PARS[3]
VAX1.EFFIC.CLIN   <- VAX1.PARS[4:6] #now the efficacy depends on the epidemiologic parameters
VAX2.EFFIC.BETA   <- VAX2.PARS[1]
VAX2.EFFIC.SEVERE <- VAX2.PARS[2]
VAX2.EFFIC.DEATH  <- VAX2.PARS[3]
VAX2.EFFIC.CLIN   <- VAX2.PARS[4:6] #now the efficacy depends on the epidemiologic parameters
############################
### FIRST DOSE INFORMATION #
############################

#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_v
VAX1.ASYMPTOMATIC.FRAC = 1 - (1 - ASYMPTOMATIC.FRAC) * (1-VAX1.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
VAX1.SEVERITY.FRAC = SEVERITY.FRAC * (1-VAX1.EFFIC.SEVERE) 
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
VAX1.DEATH.FRAC = DEATH.FRAC * (1 - VAX1.EFFIC.DEATH)
VAX1.BETA.RATE = (1- VAX1.EFFIC.BETA)

##############################
### SECOND DOSE INFORMATION ##
##############################
#### Classification of cases related to the severity of cases of once-vaccinated individuals
# Fraction of asymptomatic cases in total cases (pclin) \alpha_w
VAX2.ASYMPTOMATIC.FRAC = 1 - (1 - ASYMPTOMATIC.FRAC) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
VAX2.SEVERITY.FRAC = SEVERITY.FRAC * (1-VAX2.EFFIC.SEVERE)
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
VAX2.DEATH.FRAC = DEATH.FRAC * (1 - VAX2.EFFIC.DEATH)
VAX2.BETA.RATE = (1- VAX2.EFFIC.BETA)

