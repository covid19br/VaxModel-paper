#' beta_init_condits.
#' 
#' A function that calculates the initial conditions and the average contact rate and success between unvaccinated susceptibles and infectious individuals
#' given the age distributed new hospitalizations by week. The method assumes that the time derivatives of infected variables are zero.
#'
#' @param NEW.HOSP.AGE # New hospitalizations per week by age, a vector with 3 elements: HJ, HA, HO 'NEW.HOSP.AGE'.
#' @param PREVALENCE # Fraction of the population that already had the disease by age, a vector with 3 elements: PREVJ, PREVA, PREVO 'PREVALENCE'.
#' @param POP.DISTR # Population distributed by age, a vector with 3 elements: POPJ, POPA, POPO 'POP.DISTR'.
#' @param CONTACT.M # A contact matrix, must give as matrix 'CONTACT.M'.
#' @param EXPOSURE.PERIOD.DAYS # Average time between being infected and developing symptoms 'EXPOSURE.PERIOD.DAYS'.
#' @param SICKNESS.PERIOD.DAYS # Average time between being infectious and recovering for asymptomatic and mild 'SICKNESS.PERIOD.DAYS'.
#' @param SEVERE.PERIOD.DAYS # Average time between being infectious and recovering/dying for severe cases 'SEVERE.PERIOD.DAYS'.
#' @param CONT.REDUC.FRAC # Reduction on the expose of symptomatic (due to symptoms/quarantining) 'CONT.REDUC.FRAC'.
#' @param SEVERE.CONT.REDUC.FRAC # Reduction on the expose of severe cases (due to hospitalization) 'SEVERE.CONT.REDUC.FRAC'. 
#'                                 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
#' @param REL.INFEC.PRESYMP # relative infectiousness of pre-symptomatic individuals 'REL.INFEC.PRESYMP'.
#' @param ASYMPTOMATIC.FRAC # Fraction of asymptomatic cases in total cases'ASYMPTOMATIC.FRAC'.
#' @param SEVERITY.FRAC # Fraction of severe cases/hospitalizations in symptomatic cases (IHR) 'SEVERITY.FRAC'.
#' @param DEATH.FRAC # Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) 'DEATH.FRAC'.
#' @param V2.FRAC # Fraction of the infected people due to the second strain 'V2.FRAC'.

#'
#' @return A list where the first entry is BETA.RATE the average contact rate and success between unvaccinated susceptibles and infectious individuals
#'         and the second entry is the vector with the initial conditions.
#' @export
#' @import
#'
#' @examples

init_condits <- function(NEW.HOSP.AGE, PREVALENCE, POP.DISTR, CONTACT.M, EXPOSURE.PERIOD.DAYS,
                 SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                 SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                 SEVERITY.FRAC, DEATH.FRAC, V2.FRAC = 0.0001){
  
  POP.TOTAL <- sum(POP.DISTR)
  EXPOSURE.PERIOD.WEEKS <- EXPOSURE.PERIOD.DAYS/7.0
  SICKNESS.PERIOD.WEEKS <- SICKNESS.PERIOD.DAYS/7.0
  SEVERE.PERIOD.WEEKS <- SEVERE.PERIOD.DAYS/7.0

  # Wild strain populations 
  POP.E1 <- (1 - V2.FRAC) * NEW.HOSP.AGE * EXPOSURE.PERIOD.WEEKS / SEVERITY.FRAC # Exposed,
  POP.A1 <- (1 - V2.FRAC) * NEW.HOSP.AGE * ASYMPTOMATIC.FRAC * (1.0 - SEVERITY.FRAC) * SICKNESS.PERIOD.WEEKS / SEVERITY.FRAC # Asymptomatics
  POP.I1 <- (1 - V2.FRAC) * NEW.HOSP.AGE * (1.0 - SEVERITY.FRAC) * (1.0 - ASYMPTOMATIC.FRAC) * SICKNESS.PERIOD.WEEKS / SEVERITY.FRAC#Infectious with mild symptoms
  POP.H1 <- (1 - V2.FRAC) * NEW.HOSP.AGE * SEVERE.PERIOD.WEEKS # Hospitalized
  POP.C1 <- c(0.0, 0.0, 0.0) # Cases reported
  POP.R1 <- POP.DISTR * PREVALENCE # Recovered
  POP.D1 <- c(0.0, 0.0, 0.0) # Deaths

  # Variant strain population
  POP.E2 <- V2.FRAC * NEW.HOSP.AGE * EXPOSURE.PERIOD.WEEKS / SEVERITY.FRAC # Exposed,
  POP.A2 <- V2.FRAC * NEW.HOSP.AGE * ASYMPTOMATIC.FRAC * (1.0 - SEVERITY.FRAC) * SICKNESS.PERIOD.WEEKS / SEVERITY.FRAC # Asymptomatics
  POP.I2 <- V2.FRAC * NEW.HOSP.AGE * (1.0 - SEVERITY.FRAC) * (1.0 - ASYMPTOMATIC.FRAC) * SICKNESS.PERIOD.WEEKS / SEVERITY.FRAC#Infectious with mild symptoms
  POP.H2 <- V2.FRAC * NEW.HOSP.AGE * SEVERE.PERIOD.WEEKS # Hospitalized
  POP.C2 <- c(0.0, 0.0, 0.0) # Cases reported
  POP.R2 <- c(0.0, 0.0, 0.0) # Recovered
  POP.D2 <- c(0.0, 0.0, 0.0) # Deaths

  # Susceptibles remaining in the population
  POP.S <- POP.DISTR - (POP.E1 + POP.A1 + POP.I1 + POP.H1 + POP.R1 + POP.D1 + POP.E2 + POP.A2 + POP.I2 + POP.H2 + POP.R2 + POP.D2)

  POP0 <- c(POP.S, POP.E1, POP.A1, POP.I1, POP.H1, POP.C1, POP.R1, POP.D1, 
                   POP.E2, POP.A2, POP.I2, POP.H2, POP.C2, POP.R2, POP.D2)

  INFECTIVITY.VECTOR <- REL.INFEC.PRESYMP * POP.E1 + POP.A1 + CONT.REDUC.FRAC * POP.I1 + SEVERE.CONT.REDUC.FRAC * POP.H1

  BETA.RATE.VECTOR <- POP.E1 * POP.TOTAL / (7.0 * EXPOSURE.PERIOD.WEEKS * POP.S * CONTACT.M %*% INFECTIVITY.VECTOR)

  return(list(BETA.RATE = BETA.RATE.VECTOR[3], POP0 = POP0))
}
