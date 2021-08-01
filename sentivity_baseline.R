source('SEIR_run.R')
source('./functions/Reff_NGM.R')
library(doParallel)


VAX.INITIAL.STORAGE.NUM = 0# Number of READY vaccines on day one
VAX.PRODUCTION.RATE = 0 # Number of vaccines produced per day
MAX.VAC.RATE = 0 # Max number of vaccine applications per day
VAX.WINDOW.DAYS <- 7 * c(3)
LIST.FIRST.DOSE.REL.EFFIC = c(0)
LIST.REFF0 <- seq(0.9, 1.4, 0.1)
# Integration failure: did not find time in history
REFF.0 <- reff_ngm(sum(POP0[1:21]), POP0[1:3], CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
                   SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                   SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                   SEVERITY.FRAC, DEATH.FRAC)
REFF.0 <- Re(REFF.0$values)
OBS.VAX2.EFFIC.BETA = 0.0
OBS.VAX2.EFFIC.CLIN = 0.0
OBS.VAX2.EFFIC.SEVERE = 0.0
OBS.VAX2.EFFIC.DEATH = 0.0 ## GUESS

# efficacy of first dose is smaller
FIRST.DOSE.REL.EFFIC = 0.0
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

######substitui tudo pq eu tenho preguiça
diffEqs<- func.factory( ASYMPTOMATIC_FRAC = ASYMPTOMATIC.FRAC, BETA_RATE = BETA.RATE,
                           CONT_REDUC_FRAC = CONT.REDUC.FRAC, CONTACT_M = CONTACT.M, DEATH_FRAC =
                             DEATH.FRAC, EXPOSURE_PERIOD_DAYS = EXPOSURE.PERIOD.DAYS,
                           FIRST_DOSE_REL_EFFIC = FIRST.DOSE.REL.EFFIC, MAX_TIME_DAYS = MAX.TIME.DAYS,
                           MAX_VAC_RATE = MAX.VAC.RATE, POP_DISTR = POP.DISTR, POP_ESTADO_REL_FRAC =
                             POP.ESTADO.REL.FRAC, POP_TOTAL_NUM = POP.TOTAL.NUM, REL_INFEC_PRESYMP =
                             REL.INFEC.PRESYMP, SECOND_VAX_LOSS_FRAC = SECOND.VAX.LOSS.FRAC,
                           SEVERE_CONT_REDUC_FRAC = SEVERE.CONT.REDUC.FRAC, SEVERE_PERIOD_DAYS =
                             SEVERE.PERIOD.DAYS, SEVERITY_FRAC = SEVERITY.FRAC, SICKNESS_PERIOD_DAYS
                           = SICKNESS.PERIOD.DAYS, VAX_INITIAL_STORAGE_NUM = VAX.INITIAL.STORAGE.NUM,
                           VAX_PRODUCTION_RATE = VAX.PRODUCTION.RATE, VAX_WINDOW_DAYS =
                             VAX.WINDOW.DAYS, VAX1_ASYMPTOMATIC_FRAC = VAX1.ASYMPTOMATIC.FRAC,
                           VAX1_BETA_RATE = VAX1.BETA.RATE, VAX1_DEATH_FRAC = VAX1.DEATH.FRAC,
                           VAX1_EFFIC_CLIN = VAX1.EFFIC.CLIN, VAX1_EFFIC_DEATH = VAX1.EFFIC.DEATH,
                           VAX1_EFFIC_SEVERE = VAX1.EFFIC.SEVERE, VAX1_SEVERITY_FRAC =
                             VAX1.SEVERITY.FRAC, VAX2_ASYMPTOMATIC_FRAC = VAX2.ASYMPTOMATIC.FRAC,
                           VAX2_BETA_RATE = VAX2.BETA.RATE, VAX2_DEATH_FRAC = VAX2.DEATH.FRAC,
                           VAX2_EFFIC_CLIN = VAX2.EFFIC.CLIN, VAX2_EFFIC_DEATH = VAX2.EFFIC.DEATH,
                           VAX2_EFFIC_SEVERE = VAX2.EFFIC.SEVERE, VAX2_SEVERITY_FRAC =
                             VAX2.SEVERITY.FRAC)$diffEqs


n_cores <- detectCores() - 1
registerDoParallel(cores=6)
# result <- foreach(VAX.WINDOW.DAYS = LIST.VAX.WINDOW.DAYS, .combine = 'rbind') %:%
result <- foreach(REFF = LIST.REFF0, .combine = 'rbind') %dopar% {
    n.history <- 1e5
    step.max <- 1e6
    done <- FALSE
    # update BETA
    BETA.RATE  <- REFF / REFF.0
    # update dependent parameters
    diffEqs <- func.factory(BETA_RATE = BETA.RATE)$diffEqs
    while(done == FALSE){
      print(paste("Rodando Baseline REFF =",REFF," n.history =",n.history, " step.max = ",step.max))
      SOLUTION <- try({dopri(y = POP0, times = TIME.VECTOR, func = diffEqs,
                             parms = c(), n_history = n.history,
                             return_history = FALSE,
                             step_max_n = step.max)})
      if(class(SOLUTION) == "try-error") {####deu pau
        if(grepl("Integration failure: did not find time in history",SOLUTION,fixed=TRUE)){### tenta achar erro de histórico
          if(n.history < 2e6){####se não tiver grande demais, aumenta
            n.history = 2*n.history
          }
          else{###se estiver, finalizaa
            final.sol <- NA * POP0
            done <- TRUE
          }
        }
        else{###é erro de steps
          if(step.max < 1e8){
            step.max = 2*step.max
          }
          else{###se tiver finaliza
            final.sol <- NA * POP0
            done <- TRUE
          }
        }
        
      } else {#### não deu pau
        done <- TRUE
        write.table(SOLUTION, paste0('results/sol_baseline_', REFF, '.csv'))
        final.sol <- SOLUTION[length(TIME.VECTOR),-1]
      }
    }
    c(reff0 = REFF,
      final.sol)
  }
results <- as.data.frame(result)

write.csv(results, "results/sol_baseline2.csv", row.names=FALSE)
