#' reff_ngm.
#' 
#' R effective by Next-Generation Matrix Method.
#' 
#' The function receives P: the total population (an integer) 
#' and S: the age distribution of susceptibles (a vector with 3 elements: SJ, SA, SO) 
#' and returns an R effective by Next-Generation Matrix method. The rest are the parameters to the problem.
#' The function retrieves the eigenvalues of the NGM, usually, the first is the R_eff.
#'
#' @param P # Total population 'P'.
#' @param S # age distribution of susceptibles, a vector with 3 elements: SJ, SA, SO 'S'.
#' @param CONTACT.M # A contact matrix, must give as matrix 'CONTACT.M'.
#' @param BETA.RATE # # Average contact rate and success between unvaccinated susceptibles and infectious individuals 'BETA.RATE'.
#' @param EXPOSURE.PERIOD.DAYS # Average time between being infected and developing symptoms 'EXPOSURE.PERIOD.DAYS'.
#' @param SICKNESS.PERIOD.DAYS # Average time between being infectious and recovering for asymptomatic and mild 'SICKNESS.PERIOD.DAYS'.
#' @param SEVERE.PERIOD.DAYS # Average time between being infectious and recovering/dying for severe cases 'SEVERE.PERIOD.DAYS'.
#' @param CONT.REDUC.FRAC # Reduction on the expose of symptomatic (due to symptoms/quarantining) 'CONT.REDUC.FRAC'.
#' @param SEVERE.CONT.REDUC.FRAC # Reduction on the expose of severe cases (due to hospitalization) 'SEVERE.CONT.REDUC.FRAC'. 
#' 0 means the same level of exposure of mild cases and 1 means no expose whatsoever. 
#' @param REL.INFEC.PRESYMP # relative infectiousness of pre-symptomatic individuals 'REL.INFEC.PRESYMP'.
#' @param ASYMPTOMATIC.FRAC # Fraction of asymptomatic cases in total cases'ASYMPTOMATIC.FRAC'.
#' @param SEVERITY.FRAC # Fraction of severe cases/hospitalizations in symptomatic cases (IHR) 'SEVERITY.FRAC'.
#' @param DEATH.FRAC # Fraction of deaths in severe cases/hospitalizations of unvaccinated population (IHFR) 'DEATH.FRAC'.
#'
#' @return # Reff effective reproduction number, a number, if higher than 1 its epidemic, if lower than 1 its endemic.
#' @export
#' @import rARPACK, matlib
#'
#' @examples

library(rARPACK)
library(matlib)
reff_ngm <- function(P, S, CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
                 SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                 SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                 SEVERITY.FRAC, DEATH.FRAC){
  
  
  #NGM method
  #following the E,A,I,H order
  sml=diag(x=S)%*%(CONTACT.M)
  f=cbind(REL.INFEC.PRESYMP*sml,
                 sml,
                 (1-CONT.REDUC.FRAC)*sml,
                 (1-SEVERE.CONT.REDUC.FRAC)*sml)
  f_full=BETA.RATE/P*rbind(f,matrix(0, nrow = 9, ncol =12))
  unit<-c(1,1,1)
  v1=cbind(diag(1/EXPOSURE.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=9))
  v2=cbind(diag(-ASYMPTOMATIC.FRAC*(1-SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), diag(1/SICKNESS.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=6))
  v3=cbind(diag(-(1 - ASYMPTOMATIC.FRAC)*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3),
            diag(1/SICKNESS.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3))
  v4=cbind(diag(-(SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit),
           matrix(0,nrow=3,ncol=6),
           diag((1/SEVERE.PERIOD.DAYS)*unit))
  V=rbind(v1,v2,v3,v4)
  values<-eigs(f_full%*%inv(V), 1, which='LM', sigma=NULL )
  return(values[1])
}
