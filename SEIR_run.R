######################################
### TWO-STEP VAX LITTLE MODEL ########
######################################
# This file contain the definition of the differential equations and the
# numerical integrator.

library(dde)
library(zeallot)

# # Reads the epidemiological parameters for COVID-19 and Brazilian demographics
# source('./functions/parms_epi.R')

# Reads the parametrizations of the vaccine and vaccination campaign, calls parms_epi now
source('./functions/parms_vaccine.R')

# Parameters of the numerical integration
SIM.DURATION.DAYS = MAX.TIME.DAYS # Days of simulation
#SIM.DURATION.DAYS = 200
TIME.VECTOR <- seq(0, SIM.DURATION.DAYS) 

# Population is a named vector of tri-vectors
POP0 = c(# Never-vaccinated populations:
  POP0,
        # POP.S  = c(0.7, 0.71, 0.71)*POP.DISTR, # Susceptible
        # POP.E  = c(0.02, 0.02, 0.015)*POP.DISTR, # Exposed (infected but not infectious)
        # POP.A  = c(0.02, 0.01, 0.015)*POP.DISTR, # Asymptomatic infectious
        # POP.I  = c(0.01, 0.01, 0.01)*POP.DISTR, # Mild cases infectious
        # POP.H  = c(0,0,0), # Severe cases/hospitalization infectious
        # POP.R  = c(0.25, 0.25, 0.25)*POP.DISTR, # Recovered
        # POP.D  = c(0,0,0), # Deaths
        ## Once-vaccinated populations:
        # POP.Sv  = c(0,0,0), # Susceptible
        # POP.Ev  = c(0,0,0), # Exposed (infected but not infectious)
        # POP.Av  = c(0,0,0), # Asymptomatic infectious
        # POP.Iv  = c(0,0,0), # Mild cases infectious
        # POP.Hv  = c(0,0,0), # Severe cases/hospitalization infectious
        # POP.Rv  = c(0,0,0), # Recovered
        # POP.Dv  = c(0,0,0), # Deaths
        # ## Twice vaccinated populations:
        # POP.Sw  = c(0,0,0), # Susceptible
        # POP.Ew  = c(0,0,0), # Exposed (infected but not infectious)
        # POP.Aw  = c(0,0,0), # Asymptomatic infectious
        # POP.Iw  = c(0,0,0), # Mild cases infectious
        # POP.Hw  = c(0,0,0), # Severe cases/hospitalization infectious
        # POP.Rw  = c(0,0,0), # Recovered
        # POP.Dw  = c(0,0,0), # Deaths
        ## Total number of vaccines
        VAX = VAX.INITIAL.STORAGE.NUM
        )

func.factory <- function(
    ASYMPTOMATIC_FRAC = ASYMPTOMATIC.FRAC, BETA_RATE = BETA.RATE,
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
        VAX2.SEVERITY.FRAC) {
    # code that generates this beauty:
    # source('./parms_epi.R')
    # source('./parms_vaccine.R')
    # objs = setdiff(ls(), lsf.str())
    # paste0(lapply(objs, function(x) paste0(gsub('\\.', '_', x), " = ", x)), collapse=", ")
    #
    # WARNING: several of those are not actual parameters! They are only used
    # in parms_* to calculate other parameters

    # calculates current rate of vaccination
    OPT.VAX.RATE <- opt_vax_rate(VAX_INITIAL_STORAGE_NUM, VAX_PRODUCTION_RATE,
                                 MAX_VAC_RATE, VAX_WINDOW_DAYS,
                                 SECOND_VAX_LOSS_FRAC, MAX_TIME_DAYS)
    # interpolation for continuous time
    VAX.RATE <- interpolate.VAX.RATE(OPT.VAX.RATE)
    
    diffEqs = function(t, POP, parms) {
        c(POP.S, POP.E, POP.A, POP.I, POP.H, POP.R, POP.D,
          POP.Sv, POP.Ev, POP.Av, POP.Iv, POP.Hv, POP.Rv, POP.Dv,
          POP.Sw, POP.Ew, POP.Aw, POP.Iw, POP.Hw, POP.Rw, POP.Dw, VAX) %<-%
        split(POP, ceiling(seq_along(POP)/3))
        debug = FALSE
        if(debug && (any(c(POP.S, POP.E, POP.A, POP.I, POP.H, POP.R, POP.D,
          POP.Sv, POP.Ev, POP.Av, POP.Iv, POP.Hv, POP.Rv, POP.Dv,
          POP.Sw, POP.Ew, POP.Aw, POP.Iw, POP.Hw, POP.Rw, POP.Dw, VAX) < -1e-8) ||
           ! all(is.finite(c(POP.S, POP.E, POP.A, POP.I, POP.H, POP.R, POP.D,
          POP.Sv, POP.Ev, POP.Av, POP.Iv, POP.Hv, POP.Rv, POP.Dv,
          POP.Sw, POP.Ew, POP.Aw, POP.Iw, POP.Hw, POP.Rw, POP.Dw, VAX))))){
            cat(paste("Pop. negativa!",
                      paste("tempo t = ", t),
                      "Unvaccinated: ",
                      paste(c(POP.S, POP.E), collapse = ', '),
                      paste(c(POP.A, POP.I, POP.H), collapse = ', '),
                      paste(c(POP.R, POP.D), collapse = ', '),
                      "Vaccinated first dose: ",
                      paste(c(POP.Sv, POP.Ev), collapse = ', '),
                      paste(c(POP.Av, POP.Iv, POP.Hv), collapse = ', '),
                      paste(c(POP.Rv, POP.Dv), collapse = ', '),
                      "Vaccinated second dose: ",
                      paste(c(POP.Sw, POP.Ew), collapse = ', '),
                      paste(c(POP.Aw, POP.Iw, POP.Hw), collapse = ', '),
                      paste(c(POP.Rw, POP.Dw), collapse = ', '),
                      "Vaccines: ", VAX,
                      sep="\n"))
            stop()
        }
        # Some repeated factors:
        # Total infectious
        lambda = (BETA_RATE / POP_TOTAL_NUM) * CONTACT_M %*%
            ((POP.A + REL_INFEC_PRESYMP*POP.E  + (1-CONT_REDUC_FRAC)*POP.I  +(1-SEVERE_CONT_REDUC_FRAC)*POP.H) +
            (POP.Av + REL_INFEC_PRESYMP*POP.Ev + (1-CONT_REDUC_FRAC)*POP.Iv +(1-SEVERE_CONT_REDUC_FRAC)*POP.Hv) +
            (POP.Aw + REL_INFEC_PRESYMP*POP.Ew + (1-CONT_REDUC_FRAC)*POP.Iw +(1-SEVERE_CONT_REDUC_FRAC)*POP.Hw))
        
        SAR.NORM =  1.0/(POP.S+POP.R) # Normalization of S, E, A and R
        if(t <=  MAX_TIME_DAYS){
        V.TOTAL = VAX.RATE(t)
        }
        else{
          V.TOTAL <- 0
        }
        Vt = VAX.DISTR.RATE(V.TOTAL, 1/SAR.NORM)
        # delayed variables
        a = t - VAX_WINDOW_DAYS # Delayed time
        
        if (a <= 0){
            Vt.a = rep(0,each=3)
            POP.S.a = c(POP0['POP.S1'], POP0['POP.S2'], POP0['POP.S3'])
            POP.A.a = c(POP0['POP.A1'], POP0['POP.A2'], POP0['POP.A3'])
            POP.R.a = c(POP0['POP.R1'], POP0['POP.R2'], POP0['POP.R3'])
            SAR.NORM.a = 1.0/(POP.S.a+POP.R.a)
        } else {
            POP.S.a = ylag(a, 1:3)
            POP.A.a = ylag(a, 7:9)
            POP.R.a = ylag(a, 16:18)
            V.TOTAL.a = VAX.RATE(a)
            SAR.NORM.a = 1.0/(POP.S.a+POP.R.a)
            Vt.a = VAX.DISTR.RATE(V.TOTAL.a, 1/SAR.NORM.a)
        }
        # print(c(t,Vt,Vt.a))  
        # NEVER-VACCINATED POPULATIONS
        # Susceptible
        dS    = - POP.S * lambda - # Getting infected
                  Vt * POP.S * SAR.NORM  # Getting vaccinated
        
        # Pre-symptomatic
        dE   =   POP.S * lambda - # Getting infected
                 POP.E / EXPOSURE_PERIOD_DAYS # Becoming asymptomatic or mild case
    
        # Assymptomatic
        dA   =   ASYMPTOMATIC_FRAC*(1 - SEVERITY_FRAC)*POP.E / EXPOSURE_PERIOD_DAYS - # Becoming asymptomatic
                 POP.A / SICKNESS_PERIOD_DAYS # Recovering
                 #Vt * POP.A * SAR.NORM  # Getting vaccinated
        
        # Mild Infectious
        dI    =   (1 - ASYMPTOMATIC_FRAC) * (1 - SEVERITY_FRAC) * POP.E / EXPOSURE_PERIOD_DAYS - # Becoming mild case
                  POP.I / SICKNESS_PERIOD_DAYS # Recovering
    
        # Severe Case/Hospitalization
        dH   =   SEVERITY_FRAC * POP.E / EXPOSURE_PERIOD_DAYS - # Becoming severe case
                  POP.H / SEVERE_PERIOD_DAYS
        
        # Recovered
        dR    =   POP.A / SICKNESS_PERIOD_DAYS + # Asymptomatic recovering
                  POP.I / SICKNESS_PERIOD_DAYS + # Mild case recovering
                  (1-DEATH_FRAC)*POP.H / SEVERE_PERIOD_DAYS - # Severe case recovering
                  Vt * POP.R * SAR.NORM  # Getting vaccinated
        
        # Dead
        dD    = DEATH_FRAC *POP.H / SEVERE_PERIOD_DAYS # Severe case dying
        
        # FIRST VACCINATED POPULATIONS
        # Susceptible
        dSv   =   - POP.Sv * VAX1_BETA_RATE * lambda + # Getting infected
                  Vt * POP.S * SAR.NORM  - # Getting vaccinated once
                  (1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.S.a * SAR.NORM.a  # Getting vaccinated twice
        a1 <- POP.Sv * VAX1_BETA_RATE * lambda
        a1 <- a1[3]
        a2<- Vt * POP.S * SAR.NORM 
        a2 <- a2[3]
        # print(c(t,POP.Sv[3],a1,a2,VAX))
        # Pre-symptomatic
        dEv   =   POP.Sv * VAX1_BETA_RATE * lambda - # Getting infected
                  POP.Ev / EXPOSURE_PERIOD_DAYS # Becoming asymptomatic or mild case
    
        # Assymptomatic
        dAv   =   VAX1_ASYMPTOMATIC_FRAC*(1 - VAX1_SEVERITY_FRAC)*POP.Ev / EXPOSURE_PERIOD_DAYS - # Becoming asymptomatic
                  POP.Av / SICKNESS_PERIOD_DAYS # Recovering
                  #Vt * POP.A * SAR.NORM  - # Getting vaccinated
                  #(1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.A.a * SAR.NORM.a  # Getting vaccinated twice
        
        # Mild Infectious
        dIv    =   (1 - VAX1_ASYMPTOMATIC_FRAC) * (1 - VAX1_SEVERITY_FRAC) * POP.Ev / EXPOSURE_PERIOD_DAYS - # Becoming mild case
                   POP.Iv / SICKNESS_PERIOD_DAYS # Recovering
        
        # Severe Case/Hospitalization
        dHv   =   VAX1_SEVERITY_FRAC * POP.Ev / EXPOSURE_PERIOD_DAYS - # Becoming severe case
                  POP.Hv / SEVERE_PERIOD_DAYS
        
        # Recovered
        dRv    =   POP.Av / SICKNESS_PERIOD_DAYS + # Asymptomatic recovering
                   POP.Iv / SICKNESS_PERIOD_DAYS + # Mild case recovering
                   (1 - VAX1_DEATH_FRAC)*POP.Hv / SEVERE_PERIOD_DAYS + # Severe case recovering
                   Vt * POP.R * SAR.NORM - # Getting vaccinated
                   (1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.R.a * SAR.NORM.a  # Getting vaccinated twice
        
        # Dead
        dDv    = VAX1_DEATH_FRAC *POP.Hv / SEVERE_PERIOD_DAYS # Severe case dying
        
        # TWICE VACCINATED POPULATIONS
        # Susceptible
        dSw   =   - POP.Sw * VAX2_BETA_RATE * lambda + # Getting infected
          (1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.S.a * SAR.NORM.a  # Getting vaccinated twice
        
        # Pre-symptomatic
        dEw   =   POP.Sw * VAX2_BETA_RATE * lambda - # Getting infected
                  POP.Ew / EXPOSURE_PERIOD_DAYS # Becoming assymthomatic or mild case
    
        # Assymptomatic
        dAw   =   VAX2_ASYMPTOMATIC_FRAC *(1 - VAX2_SEVERITY_FRAC)* POP.Ew / EXPOSURE_PERIOD_DAYS - # Becoming assymtomathic
                  POP.Aw / SICKNESS_PERIOD_DAYS # Recovering
                  #(1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.A.a * SAR.NORM.a  # Getting vaccinated twice
        
        # Mild Infectious
        dIw    =   (1 - VAX2_ASYMPTOMATIC_FRAC) * (1 - VAX2_SEVERITY_FRAC) * POP.Ew / EXPOSURE_PERIOD_DAYS - # Becoming mild case
                  POP.Iw / SICKNESS_PERIOD_DAYS # Recovering
        
        # Severe Case/Hospitalization
        dHw   =    VAX2_SEVERITY_FRAC * POP.Ew / EXPOSURE_PERIOD_DAYS - # Becoming severe case
                  POP.Hw / SEVERE_PERIOD_DAYS
        
        # Recovered
        dRw    =   POP.Aw / SICKNESS_PERIOD_DAYS + # Asymptomatic recovering
          POP.Iw / SICKNESS_PERIOD_DAYS + # Mild case recovering
          (1 - VAX2_DEATH_FRAC)*POP.Hw / SEVERE_PERIOD_DAYS + # Severe case recovering
          (1 - SECOND_VAX_LOSS_FRAC) * Vt.a * POP.R.a * SAR.NORM.a  # Getting vaccinated twice
        
        # Dead
        dDw    = VAX2_DEATH_FRAC *POP.Hw / SEVERE_PERIOD_DAYS # Severe case dying
        
        # VACCINE
        dVAX  = VAX_PRODUCTION_RATE - sum(Vt) - (1 - SECOND_VAX_LOSS_FRAC) * sum(Vt.a)
        # debug = TRUE
        
        if(debug && ! all(is.finite(c(dS, dE, dA, dI, dH, dR, dD,
          dSv, dEv, dAv, dIv, dHv, dRv, dDv,
          dSw, dEw, dAw, dIw, dHw, dRw, dDw, dVAX)))){
            cat(paste(" *** Derivadas infinitas/NaN! ***",
                      paste("tempo t = ", t),
                      "Unvaccinated: ",
                      paste(c(dS, dE), collapse = ', '),
                      paste(c(dA, dI, dH), collapse = ', '),
                      paste(c(dR, dD), collapse = ', '),
                      "Vaccinated first dose: ",
                      paste(c(dSv, dEv), collapse = ', '),
                      paste(c(dAv, dIv, dHv), collapse = ', '),
                      paste(c(dRv, dDv), collapse = ', '),
                      "Vaccinated second dose: ",
                      paste(c(dSw, dEw), collapse = ', '),
                      paste(c(dAw, dIw, dHw), collapse = ', '),
                      paste(c(dRw, dDw), collapse = ', '),
                      "Vaccines: ", dVAX,
                      SAR.NORM, SAR.NORM.a,"",
                      sep="\n"))
            cat(paste(" *** Populações: ***",
                      "Unvaccinated: ",
                      paste(c(POP.S, POP.E), collapse = ', '),
                      paste(c(POP.A, POP.I, POP.H), collapse = ', '),
                      paste(c(POP.R, POP.D), collapse = ', '),
                      "Vaccinated first dose: ",
                      paste(c(POP.Sv, POP.Ev), collapse = ', '),
                      paste(c(POP.Av, POP.Iv, POP.Hv), collapse = ', '),
                      paste(c(POP.Rv, POP.Dv), collapse = ', '),
                      "Vaccinated second dose: ",
                      paste(c(POP.Sw, POP.Ew), collapse = ', '),
                      paste(c(POP.Aw, POP.Iw, POP.Hw), collapse = ', '),
                      paste(c(POP.Rw, POP.Dw), collapse = ', '),
                      "Vaccines: ", VAX,
                      sep="\n"))
    
            stop()
        }
    
        return(c(dS, dE, dA, dI, dH, dR, dD, dSv, dEv, dAv, dIv, dHv, dRv, dDv, dSw, dEw, dAw, dIw, dHw, dRw, dDw, dVAX))
    }

    return(list(OPT.VAX.RATE = OPT.VAX.RATE,
                VAX.RATE = VAX.RATE,
                diffEqs = diffEqs))
}


source('./functions/Reff_NGM.R')
# SIM.DURATION.DAYS = 120 # Days of simulation
# #SIM.DURATION.DAYS = 200
# TIME.VECTOR <- seq(0, SIM.DURATION.DAYS) 
# REFF.0 <- reff_ngm(sum(POP0[1:21]), POP0[1:3], CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
#               SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
#               SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
#               SEVERITY.FRAC, DEATH.FRAC)
# REFF.0 <- Re(REFF.0$values)
# BETA.RATE <- 1 / REFF.0
# diffEqs <- func.factory(BETA_RATE = BETA.RATE, MAX_TIME_DAYS = MAX.TIME.DAYS)$diffEqs
# SOLUTION <- dopri(y = POP0, times = TIME.VECTOR, parms = c(), func = diffEqs,
#                  n_history = 4e6, return_history = FALSE,step_max_n = 1e6)
