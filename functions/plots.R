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