source('SEIR_run.R')
source('./functions/Reff_NGM.R')
library(doParallel)

#LIST.VAX.WINDOW.DAYS <- 7 * c(3, 12)
LIST.VAX.WINDOW.DAYS <- 7 * c(3, 7, 12)
#LIST.VAX.WINDOW.DAYS <- 7 *c(3) 
# efficacy of first dose is smaller
LIST.FIRST.DOSE.REL.EFFIC = seq(0.0, 1.0, 0.1)
#LIST.FIRST.DOSE.REL.EFFIC = c(0.5)
LIST.REFF0 <- seq(0.9, 1.4, 0.1)
#LIST.REFF0 <- c(0.9)
################ UNCOMMENT THIR TO RUN CORONAVAC ########################
 VACCINE <- "CoronaVac"
 OBS.VAX2.EFFIC.BETA = 0.0
 OBS.VAX2.EFFIC.CLIN = 0.5
 OBS.VAX2.EFFIC.SEVERE = 0.83
 OBS.VAX2.EFFIC.DEATH = 0.95 ## GUESS
################ UNCOMMENT THIS TO RUN PFIZER VACCINE ###################
#VACCINE <- "Pfizer"
#OBS.VAX2.EFFIC.BETA = 0.90
#OBS.VAX2.EFFIC.CLIN = 0.94
#OBS.VAX2.EFFIC.SEVERE = 0.87
#OBS.VAX2.EFFIC.DEATH = 0.98 ## GUESS
################ UNCOMMENT THIS TO RUN ASTRAZENECA ######################
# VACCINE <- "AstraZeneca"
# OBS.VAX2.EFFIC.BETA = 0.599
# OBS.VAX2.EFFIC.CLIN = 0.813
# OBS.VAX2.EFFIC.SEVERE = 0.900 ## GUESS
# OBS.VAX2.EFFIC.DEATH = 0.950 ## GUESS
#########################################################################
VAX2.PARS <- OBSERVED.EFFICACY.TO.PARAMETERS(SEVERITY.FRAC,DEATH.FRAC,ASYMPTOMATIC.FRAC,
                                             OBS.VAX2.EFFIC.BETA,OBS.VAX2.EFFIC.SEVERE,
                                             OBS.VAX2.EFFIC.DEATH,OBS.VAX2.EFFIC.CLIN)
VAX2.EFFIC.BETA   <- VAX2.PARS[1]
VAX2.EFFIC.SEVERE <- VAX2.PARS[2]
VAX2.EFFIC.DEATH  <- VAX2.PARS[3]
VAX2.EFFIC.CLIN   <- VAX2.PARS[4:6] #now the efficacy depends on the epidemiologic parameters
# Fraction of asymptomatic cases in total cases (pclin) \alpha_w
VAX2.ASYMPTOMATIC.FRAC = 1 - (1 - ASYMPTOMATIC.FRAC) * (1-VAX2.EFFIC.CLIN)
# Fraction of severe cases in symptomatic cases (IHR) \sigma_v
VAX2.SEVERITY.FRAC = SEVERITY.FRAC * (1-VAX2.EFFIC.SEVERE)
# Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_w
VAX2.DEATH.FRAC = DEATH.FRAC * (1 - VAX2.EFFIC.DEATH)
VAX2.BETA.RATE = (1- VAX2.EFFIC.BETA)
diffEqs <- func.factory(VAX2_BETA_RATE = VAX2.BETA.RATE,
                        VAX2_ASYMPTOMATIC_FRAC = VAX2.ASYMPTOMATIC.FRAC,
                        VAX2_SEVERITY_FRAC = VAX2.SEVERITY.FRAC,
                        VAX2_DEATH_FRAC = VAX2.DEATH.FRAC)$diffEqs
#########################################################################
VAX.INITIAL.STORAGE.NUM = 0 # Number of READY vaccines on day one
diffEqs <- func.factory(VAX_INITIAL_STORAGE_NUM = VAX.INITIAL.STORAGE.NUM)$diffEqs
#########################################################################
REFF.0 <- reff_ngm(sum(POP0[1:21]), POP0[1:3], CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
                   SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                   SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                   SEVERITY.FRAC, DEATH.FRAC)
REFF.0 <- Re(REFF.0$values)
#########################################################################
args <- commandArgs(trailingOnly = TRUE)
REFF <- as.numeric(args[1])
REFF <- 0.1*REFF
VAX.WINDOW.DAYS <- as.numeric(args[2])
FIRST.DOSE.REL.EFFIC <- as.numeric(args[3])/100
#########################################################################
funcionou <- FALSE

            print(paste("Rodando VAX.WINDOW.DAYS =", VAX.WINDOW.DAYS, "REFF =",
                  REFF, "FIRST.DOSE.REL.EFFIC =", FIRST.DOSE.REL.EFFIC))
            # update BETA
            BETA.RATE  <- REFF / REFF.0
            # update dependent parameters
            OBS.VAX1.EFFIC.BETA   = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.BETA
            OBS.VAX1.EFFIC.CLIN   = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.CLIN
            OBS.VAX1.EFFIC.SEVERE = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.SEVERE
            OBS.VAX1.EFFIC.DEATH  = FIRST.DOSE.REL.EFFIC * OBS.VAX2.EFFIC.DEATH
            VAX1.PARS <- OBSERVED.EFFICACY.TO.PARAMETERS(SEVERITY.FRAC,DEATH.FRAC,ASYMPTOMATIC.FRAC,
                                                         OBS.VAX1.EFFIC.BETA,OBS.VAX1.EFFIC.SEVERE,
                                                         OBS.VAX1.EFFIC.DEATH,OBS.VAX1.EFFIC.CLIN)
            VAX1.EFFIC.BETA   <- VAX1.PARS[1]
            VAX1.EFFIC.SEVERE <- VAX1.PARS[2]
            VAX1.EFFIC.DEATH  <- VAX1.PARS[3]
            VAX1.EFFIC.CLIN   <- VAX1.PARS[4:6] #now the efficacy depends on the epidemiologic parameters
            # Fraction of asymptomatic cases in total cases (pclin) \alpha_v
            VAX1.ASYMPTOMATIC.FRAC = 1 - (1 - ASYMPTOMATIC.FRAC) * (1-VAX1.EFFIC.CLIN)
            # Fraction of severe cases in symptomatic cases (IHR) \sigma_v
            VAX1.SEVERITY.FRAC = SEVERITY.FRAC * (1-VAX1.EFFIC.SEVERE) 
            # Fraction of deaths in severe cases/hospitalizations of vaccinated population (IHFR) \mu_v
            VAX1.DEATH.FRAC = DEATH.FRAC * (1 - VAX1.EFFIC.DEATH)
            VAX1.BETA.RATE = (1- VAX1.EFFIC.BETA)


            n.history <- 4e5
            step.max <- 2e6
            done <- FALSE
            diffEqs <- func.factory(VAX_WINDOW_DAYS = VAX.WINDOW.DAYS,
                                    BETA_RATE = BETA.RATE,
                                    VAX1_ASYMPTOMATIC_FRAC = VAX1.ASYMPTOMATIC.FRAC,
                                    VAX1_SEVERITY_FRAC = VAX1.SEVERITY.FRAC,
                                    VAX1_DEATH_FRAC = VAX1.DEATH.FRAC,
                                    )$diffEqs
            while(done == FALSE){
              print(paste("Rodando VAX.WINDOW.DAYS =", VAX.WINDOW.DAYS, "REL.EFFIC =",
                          FIRST.DOSE.REL.EFFIC, "REFF =", REFF," n.history = ",n.history,
                          " step.max = ",step.max))
              SOLUTION <- try({dopri(y = POP0, times = TIME.VECTOR, func = diffEqs,
                                     parms = c(), n_history = n.history,
                                     return_history = FALSE,
                                     step_max_n = step.max)})
              if(class(SOLUTION) == "try-error") {####deu pau
                if(grepl("Integration failure: did not find time in history",SOLUTION,fixed=TRUE)){### tenta achar erro de histórico
                  if(n.history < 2e6){####se não tiver grande demais, aumenta
                    n.history = 2*n.history
                  }
                  else{###se estiver, finaliza
                    final.sol <- NA * POP0
                    print(paste("Resolver primeira tentativa deu pau por falta de histórico. VAX.WINDOW.DAYS =", VAX.WINDOW.DAYS, "REL.EFFIC =",
                                FIRST.DOSE.REL.EFFIC, "REFF =", REFF))
                    done <- TRUE
                  }
                }
                else{###é erro de steps
                  if(step.max < 1e8){
                    step.max = 2*step.max
                  }
                  else{###se tiver finaliza
                    final.sol <- NA * POP0
                    print(paste("Resolver primeira tentativa deu pau por excesso de steps. VAX.WINDOW.DAYS =", VAX.WINDOW.DAYS, "REL.EFFIC =",
                                FIRST.DOSE.REL.EFFIC, "REFF =", REFF))
                    done <- TRUE
                  }
                }
                
              } else {#### não deu pau
                done <- TRUE
                write.table(SOLUTION, paste0('results/sol_reff_',VACCINE,'_', VAX.WINDOW.DAYS, '_', REFF, '_',
                                             FIRST.DOSE.REL.EFFIC, '.csv'))
                final.sol <- SOLUTION[length(TIME.VECTOR),-1]
		funcionou <- TRUE
              }
            }###fim do while
if(funcionou){print(c(paste("funcionou para reff",REFF,"eficacia",FIRST.DOSE.REL.EFFIC,"e janela",VAX.WINDOW.DAYS)))}
