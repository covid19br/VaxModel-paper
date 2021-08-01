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
#' @import matlib
#'
#' @examples
reff_ngm <- function(P, S, CONTACT.M, BETA.RATE = 1., EXPOSURE.PERIOD.DAYS,
                     SICKNESS.PERIOD.DAYS, SEVERE.PERIOD.DAYS, CONT.REDUC.FRAC,
                     SEVERE.CONT.REDUC.FRAC, REL.INFEC.PRESYMP, ASYMPTOMATIC.FRAC,
                     SEVERITY.FRAC, DEATH.FRAC){
  
  
  #NGM method
  #following the E,A,I,H order
  sml=t(t(CONTACT.M)%*%diag(x=S))
  f=cbind(REL.INFEC.PRESYMP*sml,
          sml,
          (1-CONT.REDUC.FRAC)*sml,
          (1-SEVERE.CONT.REDUC.FRAC)*sml)
  f_full=BETA.RATE/P*rbind(f,matrix(0, nrow = 9, ncol =12))
  unit<-c(1,1,1)
  v1=cbind(diag(1/EXPOSURE.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=9))
  v2=cbind(diag(-ASYMPTOMATIC.FRAC/EXPOSURE.PERIOD.DAYS*unit), diag(1/SICKNESS.PERIOD.DAYS*unit),matrix(0,nrow=3,ncol=6))
  v3=cbind(diag(-(1 - ASYMPTOMATIC.FRAC)*(1 - SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3),
           diag(1/SICKNESS.PERIOD.DAYS*unit), matrix(0,nrow=3,ncol=3))
  v4=cbind(diag(-(1 - ASYMPTOMATIC.FRAC)*(SEVERITY.FRAC)/EXPOSURE.PERIOD.DAYS*unit),
           matrix(0,nrow=3,ncol=6),
           diag((1/SEVERE.PERIOD.DAYS)*unit))
  V=rbind(v1,v2,v3,v4)
  values<-eigen(f_full%*%inv(V))$values
  return(values[1])
}

#' aggregate_age_matrix.
#' 
#' A function to operate the average contacts with weighting between the age strucuture of a population.
#'
#' @param matrix N by N matrix with contacts between age bins 'matrix'.
#' @param aggregate_indices A list of vector of indexes of age bins that should be aggregated 'aggregate_indices'.
#' @param age_structure A vector of size N with population distribution to act as weights in a weighted average 'age_structure'.
#'
#' @return A new matrix with size \code{length(aggregate_index)} by \code{length(aggregate_index)} with average contacts weighted by population structure
#' @export
#'
#' @examples
aggregate_age_matrix <- function(matrix, aggregate_indices, age_structure){
  if (missing(age_structure))
  {
    age_structure = rep(c(1),each=dim(matrix)[1])
  }
  if(dim(matrix)[1] != dim(matrix)[2]){
    stop("matrix is not a square matrix")
  }
  if(length(age_structure) < dim(matrix)[1]){
    stop("age_structure is smaller than the linear dimension of matrix")
  }
  if(length(age_structure) > dim(matrix)[1]){
    warning("age_structure is bigger than linear size of matrix, aggregating the last values in one")
    age_structure[length(matrix)[1]] =  sum(age_structure[length(matrix)[1]:length(age_structure)])
  }
  new_matrix = matrix(0,nrow=length(aggregate_indices),ncol=length(aggregate_indices))
  for (i in 1:length(aggregate_indices)){
    for(j in 1:length(aggregate_indices)){
      sum_age_struc = 0;
      for(k in aggregate_indices[[i]]){
        for(l in aggregate_indices[[j]]){
          new_matrix[i,j] = new_matrix[i,j] + age_structure[k]*matrix[k,l]
        }
        sum_age_struc = sum_age_struc + age_structure[k]
      }
      new_matrix[i,j] = new_matrix[i,j]/sum_age_struc
    }
  }
  return(new_matrix)
}

#' plot.deaths.
#'
#' @param results 
#'
#' @return
#' @export
#' @import dplyr, ggplot2
#'
#' @examples
plot.deaths <- function(results) {
  Deaths <- results %>%
    select("vax.window.days", "first.dose.rel.effic", "reff0", contains("POP.D")) %>%
    mutate(total.deaths = rowSums(.[3:11]),
           vax.window.days = as.factor(vax.window.days)) %>%
    select("vax.window.days", "first.dose.rel.effic", "reff0", "total.deaths")
  
  p <- ggplot(Deaths) +
    geom_point(aes(x=first.dose.rel.effic, y=total.deaths,
                   col=vax.window.days)) +
    scale_color_discrete(name="spacing between\nvaccine doses",
                         labels=c("3 weeks" ,"12 weeks")) +
    scale_x_continuous(breaks=seq(0., 1., 0.1)) +
    facet_wrap(~reff0) +
    labs(x = "relative efficacy of first dose",
         y = "total deaths")
  print(p)
}

#' plot.simulation
#'
#' @param solution 
#' @param cols 
#'
#' @return
#' @export
#' @import dplyr, ggplot2
#'
#' @examples
plot.simulation <- function(solution, cols=c("Su", "Iu", "Ru", "Du", "Sv",
                                             "Iv", "Rv", "Dv", "Sw", "Iw",
                                             "Rw", "Dw")) {
  sol.df <- data.frame(
    Su = rowSums(solution[,(2:4)]),
    Iu = rowSums(solution[,(5:16)]),
    Ru = rowSums(solution[,(17:19)]),
    Du = rowSums(solution[,(20:22)]),
    Sv = rowSums(solution[,(23:25)]),
    Iv = rowSums(solution[,(26:37)]),
    Rv = rowSums(solution[,(38:40)]),
    Dv = rowSums(solution[,(41:43)]),
    Sw = rowSums(solution[,(44:46)]),
    Iw = rowSums(solution[,(47:59)]),
    Rw = rowSums(solution[,(59:61)]),
    Dw = rowSums(solution[,(63:64)])
  )
  sol.df <- sol.df[,cols]
  matplot(solution[,1], sol.df, type='l', lwd = 5)
  nn <- ncol(sol.df)
  legend("topright", colnames(sol.df), col=seq_len(nn), lty=seq_len(nn), cex=0.8)
}

#' opt_vax_rate.
#' 
#' Optimal vaccination rate function, 
#' its optimaze the problem of rollout vaccine solving a by non-linear programming solver.
#'
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'. 
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @import lpsolve
#'
#' @return
#' @export 
#'
#' @examples
opt_vax_rate <- function(VAX.INITIAL.STORAGE.NUM, VAX.PRODUCTION.RATE,
                         MAX.VAC.RATE, VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                         MAX.TIME.DAYS) {
  times <- seq(0, MAX.TIME.DAYS)
  
  alpha = 1 - SECOND.VAX.LOSS.FRAC
  V.T = max(VAX.WINDOW.DAYS * (alpha * MAX.VAC.RATE - VAX.PRODUCTION.RATE), 0)
  
  # boundary-free case
  if (VAX.INITIAL.STORAGE.NUM + MAX.TIME.DAYS * (VAX.PRODUCTION.RATE - MAX.VAC.RATE) -
      (MAX.TIME.DAYS - VAX.WINDOW.DAYS) * alpha * MAX.VAC.RATE >= V.T) {
    return(rep(MAX.VAC.RATE, length(times)))
  }
  
  # feasibility is trivial: V(T) > 0 if VAX.RATE == 0
  
  # solution using LP
  N <- length(times)
  M <- matrix(0, N, N)
  M[row(M)-col(M) >= 0] = 1
  M.2 <- matrix(0, N, N)
  M.2[seq(1 + VAX.WINDOW.DAYS, N),] <- M[seq(1, N - VAX.WINDOW.DAYS),]
  
  M.total <- - M - alpha*M.2
  
  # objective
  C.T <- colSums(M.total)
  ## inequality constraints (Ax <= b)
  A <- M.total
  # diagonal matrix offset by window
  offdiag <- matrix(0, nrow = N, ncol = N)
  offdiag[row(offdiag) - col(offdiag) == VAX.WINDOW.DAYS] = 1
  A <- rbind(- M.total,
             diag(nrow = N, ncol = N) + alpha * offdiag,
             M[N,]
  )
  b <- c(VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * times,
         rep(MAX.VAC.RATE, N),
         (VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * (times[N] + VAX.WINDOW.DAYS))/(1+alpha)
  )
  xopt <-  lp(direction="min",
              objective.in = C.T,
              const.mat = A,
              const.dir = rep("<=", length(b)),
              const.rhs = b)
  #print(paste("Vacinas restantes no tempo final:",
  #            sum(M.total[nrow(M.total),] * xopt$solution) + b[length(times)]))
  return(xopt$solution)
}

#' plot_vac_schedule.
#' 
#' A function to plot the rollout vaccine after the optimal solution.
#' 
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#' @param VAX.INITIAL.STORAGE.NUM # Number of READY vaccines on day one 'VAX.INITIAL.STORAGE.NUM'.
#' @param VAX.PRODUCTION.RATE # Number of vaccines produced per day 'VAX.PRODUCTION.RATE'.
#' @param MAX.VAC.RATE # Max number of vaccine applications per day 'MAX.VAC.RATE'.
#' @param VAX.WINDOW.DAYS # Time window between first and second vaccines 'VAX.WINDOWS.DAYS'. 
#' Its varies with the vaccine type.
#' @param SECOND.VAX.LOSS.FRAC # Fraction of people that take the first but not the second dose 'SECODNE.VAX.LOSS.FRAC'.
#' @param MAX.TIME.DAYS # time until end of vaccination schedule, in days 'MAX.TIME.DAYS'.
#' @import lpsolve
#'
#' @return
#' @export
#'
#' @examples
plot_vac_schedule <- function(OPT.VAX.RATE, VAX.INITIAL.STORAGE.NUM,
                              VAX.PRODUCTION.RATE, MAX.VAC.RATE,
                              VAX.WINDOW.DAYS, SECOND.VAX.LOSS.FRAC,
                              MAX.TIME.DAYS) {
  times <- seq(0, MAX.TIME.DAYS)
  alpha = 1 - SECOND.VAX.LOSS.FRAC
  V <- VAX.INITIAL.STORAGE.NUM + VAX.PRODUCTION.RATE * times -
    cumsum(OPT.VAX.RATE) -
    alpha * c(rep(0, VAX.WINDOW.DAYS),
              cumsum(OPT.VAX.RATE[seq(1, length(times) - VAX.WINDOW.DAYS)]))
  
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(times, OPT.VAX.RATE,
       xlab="time (days)",
       ylab="vaccination rate (doses/day)",
       type='l')
  lines(times, alpha * c(rep(0, VAX.WINDOW.DAYS),
                         OPT.VAX.RATE[-seq(length(times)+1-VAX.WINDOW.DAYS, length(times))]),
        type='l', lty=2)
  par(new = TRUE)
  plot(times, V, type = "l", axes = FALSE, bty = "n", lty = 3, xlab = "", ylab = "", col="blue")
  axis(side=4, at = pretty(round(range(V), -5)))
  mtext("Vaccine doses stored", side=4, line=3, col="blue")
}

#' interpolate.VAX.RATE.
#' 
#' A function to interpolate vaccine rate rollout give the vaccine rate
#'
#' @param OPT.VAX.RATE # Optimal vaccine rate from the `opt_vax_rate` 'OPT.VAX.RATE'.
#'
#' @return
#' @export
#'
#' @examples
interpolate.VAX.RATE <- function(OPT.VAX.RATE) {
  VAX.RATE <- function(t){
    if (t <= 0)
      return(0)
    tu <- ceiling(t)
    x <- tu - t
    s = 0.1
    if (x > 1-s && tu > 1)
      return((x)*OPT.VAX.RATE[tu] +
               (1-x)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu-1])/2)
    if (x < s)
      return((1-x/s)*OPT.VAX.RATE[tu] +
               (x/s)*(OPT.VAX.RATE[tu] + OPT.VAX.RATE[tu+1])/2)
    return(OPT.VAX.RATE[tu])
  }
}