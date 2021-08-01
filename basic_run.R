source("SEIR_run.R")

SIM.DURATION.DAYS = 120 # Days of simulation
#SIM.DURATION.DAYS = 200
TIME.VECTOR <- seq(0, SIM.DURATION.DAYS) 
REFF.0 <- reff_ngm(sum(POP0[1:21]), POP0[1:3], CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
                   SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                   SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                   SEVERITY.FRAC, DEATH.FRAC)
REFF.0 <- Re(REFF.0$values)
BETA.RATE <- 1 / REFF.0
diffEqs <- func.factory(BETA_RATE = BETA.RATE, MAX_TIME_DAYS = MAX.TIME.DAYS)$diffEqs
SOLUTION <- dopri(y = POP0, times = TIME.VECTOR, parms = c(), func = diffEqs,
                  n_history = 4e6, return_history = FALSE,step_max_n = 1e6)

source("functions/plots.R")

plot.simulation(SOLUTION)

plot.simulation(SOLUTION,c("Sv","Sw","Rv","Rw"))

plot.simulation(SOLUTION,c("Du","Dv","Dw"))
