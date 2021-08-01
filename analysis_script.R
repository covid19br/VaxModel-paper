library(ggplot2)
library(dplyr)

# library(reconPlots)
vax.days <- c(21,49,84)
reffs <- c(0.9,1.0,1.1,1.2,1.3,1.4)
rel.effics <- seq(0.0,1.0,0.1)
death.compartments <- c("POP.D1","POP.D2","POP.D3","POP.Dv1","POP.Dv2","POP.Dv3","POP.Dw1","POP.Dw2","POP.Dw3")
aggregated_results <- c()
for (Days in vax.days) {
  for(R.eff in reffs){
    for (releffic in rel.effics) {
      file.name <- paste0("results/coronavac/sol_reff_CoronaVac_",Days,"_",R.eff,"_",releffic,".csv")
      # file.name <- paste(file.name,".csv",sep="")
      val <- read.table(file.name)
      deaths <- as.numeric(sum(val[301,death.compartments]))
      result <- c(releffic,deaths,Days/7,R.eff,"CoronaVac")
      result[1:4] <- as.numeric(result[1:4])
      # result$Vaccine <- "CoronaVac"
      
      aggregated_results <- rbind(aggregated_results,result)
    }
  }
}
for (Days in vax.days) {
  for(R.eff in reffs){
    for (releffic in rel.effics) {
      file.name <- paste0("results/astrazeneca/sol_reff_AstraZeneca_",Days,"_",R.eff,"_",releffic,".csv")
      # file.name <- paste(file.name,".csv",sep="")
      val <- read.table(file.name)
      deaths <- as.numeric(sum(val[301,death.compartments]))
      result <- c(releffic,deaths,Days/7,R.eff,"AZD1222")
      result[1:4] <- as.numeric(result[1:4])
      # result$Vaccine <- "CoronaVac"
      
      aggregated_results <- rbind(aggregated_results,result)
    }
  }
}
for (Days in vax.days) {
  for(R.eff in reffs){
    for (releffic in rel.effics) {
      file.name <- paste0("results/pfizer/sol_reff_Pfizer_",Days,"_",R.eff,"_",releffic,".csv")
      # file.name <- paste(file.name,".csv",sep="")
      val <- read.table(file.name)
      deaths <- as.numeric(sum(val[301,death.compartments]))
      result <- c(releffic,deaths,Days/7,R.eff,"BNT162b2")
      result[1:4] <- as.numeric(result[1:4])
      # result$Vaccine <- "CoronaVac"
      
      aggregated_results <- rbind(aggregated_results,result)
    }
  }
}
aggregated_results <- data.frame(aggregated_results)
colnames(aggregated_results) <- c("Relative.Efficacy","Deaths","Weeks","R.eff","Vaccine")
# aggregated_results$Weeks<- as.character(aggregated_results$Weeks)
rownames(aggregated_results) <- c()

####reading baseline (no vaccine)#####
baseline <- read.csv("results/sol_baseline2.csv")
baseline.df <- c()
for (reff in baseline$reff0){
  val <- filter(baseline, reff0 == reff)  
  deaths <- as.numeric(sum(val[,death.compartments]))
  baseline.df <- rbind(baseline.df, c(reff,deaths))
}
baseline.df <- as.data.frame(baseline.df)
rownames(baseline.df) <- c()
colnames(baseline.df) <- c("R.eff", "Deaths.baseline")
#create percentage column
aggregated_results$Relative.Efficacy <- unlist(aggregated_results$Relative.Efficacy)
aggregated_results$Weeks <- factor(unlist(aggregated_results$Weeks),levels = c(3,7,12))
aggregated_results$R.eff <- as.numeric(aggregated_results$R.eff)
aggregated_results$Deaths <- as.numeric(aggregated_results$Deaths)

aggregated_results <- left_join(aggregated_results,baseline.df,by="R.eff")
aggregated_results <-  mutate(aggregated_results,Reduction = 1-Deaths/Deaths.baseline)
aggregated_results$Deaths.baseline <- as.numeric(aggregated_results$Deaths.baseline)
aggregated_results<- as.data.frame(aggregated_results)
######################################
aggregated_results$Vaccine <- factor(aggregated_results$Vaccine,levels = c("CoronaVac","AZD1222","BNT162b2"))

# labels_reff = list("R.eff: 0.9" = "Rt = 0.9","Rt = 1.0","Rt = 1.1","Rt = 1.2","Rt = 1.3","Rt = 1.4")
colnames(aggregated_results)[4] <- "Rt"
g <- ggplot(aggregated_results,aes(x = Relative.Efficacy,y=Reduction,col=factor(Weeks),group = Weeks))+labs(x="Relative Efficacy of First Dose",y="Reduction in Total Number of Deaths",col="Weeks")+
  facet_grid( Rt~ Vaccine, labeller = labeller(Rt = label_both))+geom_point(size = 1.5)+geom_line(size = 1)+ scale_y_continuous(labels = scales::percent)+
    theme_bw() + theme(text = element_text(size =18))
# #####compute alpha critico
# alpha_crit <- c()
# weeks <- c(7,12) 
# for (reff in reffs) {
#   baseline <- select(filter(aggregated_results,R.eff == reff,Weeks == 3),c("Relative.Efficacy","Deaths"))
#   colnames(baseline) <- c("x","y")
#   for (week in weeks) {
#     comparison <- select(filter(aggregated_results,R.eff == reff,Weeks == week),c("Relative.Efficacy","Deaths"))
#     colnames(comparison) <- c("x","y")
#     value <- curve_intersect(baseline, comparison)
#     alpha_crit <- rbind(alpha_crit,c(week,reff,value$x))
#   }
# }
# alpha_crit <- as.data.frame(alpha_crit)
# colnames(alpha_crit) <- c("Weeks","R.eff","Alpha")
# g_alpha <- ggplot(alpha_crit,aes(x = Weeks,y=Alpha,col=R.eff)) + geom_point(size=3)

# alpha_grid <- c()
# for (alpha in seq(0.1,0.9,0.1)) {
#   for (reff in reffs){
#     dados <- select(filter(aggregated_results,Relative.Efficacy == alpha,R.eff == reff),c("Deaths","Weeks"))
#     vals <- c(alpha, reff, dados$Weeks[which.min(dados$Deaths)])
#     alpha_grid <- rbind(alpha_grid,vals)
#   }
# }
# alpha_grid <- as.data.frame(alpha_grid)
# colnames(alpha_grid) <- c("Alpha","R.eff","Weeks")
# rownames(alpha_grid) <- c()
# grid_plot <- ggplot(alpha_grid,aes(x=Alpha,y=R.eff,col=factor(Weeks))) + geom_tile(aes(fill=factor(Weeks)))
