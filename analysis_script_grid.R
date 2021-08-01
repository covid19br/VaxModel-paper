library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(reshape2)
library(metR)
library(pracma)
source("functions/parms_epi.R")
LIST.PRODUCTION.RATE <- POP.ESTADO.REL.FRAC*seq(5e5,2e6,length.out = 15)
LIST.PRODUCTION.RATE<- LIST.PRODUCTION.RATE[1:10]
list2 <- LIST.PRODUCTION.RATE[1]-(1:4)*(LIST.PRODUCTION.RATE[2] - LIST.PRODUCTION.RATE[1])
#list2 <- POP.ESTADO.REL.FRAC*seq(1e4,5e5,length.out = 5)
LIST.PRODUCTION.RATE <- c(list2,LIST.PRODUCTION.RATE)
LIST.PRODUCTION.RATE <- sort(LIST.PRODUCTION.RATE)
LIST.VAX.WINDOW.DAYS <- 7 * seq(3, 12)
LIST.FIRST.DOSE.REL.EFFIC = seq(0, 14)/14

death.compartments <- c("POP.D1","POP.D2","POP.D3","POP.Dv1","POP.Dv2","POP.Dv3","POP.Dw1","POP.Dw2","POP.Dw3")

values <- c()
for (PROD.RATE in LIST.PRODUCTION.RATE) {
  for(VAX.WINDOW.DAYS in LIST.VAX.WINDOW.DAYS){
    for(REL.EFFIC in LIST.FIRST.DOSE.REL.EFFIC){
      filename <- paste("./results/pfizer/sol_Pfizer_",PROD.RATE,"_",VAX.WINDOW.DAYS,"_1.1_",REL.EFFIC,".csv",sep = "")
      val <- try(read.table(filename))
      if(class(val) != "try-error"){
        deaths <- as.numeric(sum(val[301,death.compartments]))
        resultado <- as.numeric(c(PROD.RATE,VAX.WINDOW.DAYS,REL.EFFIC,deaths))
        resultado$Vaccine <- "BNT162b2"
        values <- rbind(values,resultado)
      }
    }
  }
}
for (PROD.RATE in LIST.PRODUCTION.RATE) {
  for(VAX.WINDOW.DAYS in LIST.VAX.WINDOW.DAYS){
    for(REL.EFFIC in LIST.FIRST.DOSE.REL.EFFIC){
      filename <- paste("./results/coronavac/sol_CoronaVac_",PROD.RATE,"_",VAX.WINDOW.DAYS,"_1.1_",REL.EFFIC,".csv",sep = "")
      val <- try(read.table(filename))
      if(class(val) != "try-error"){
        deaths <- as.numeric(sum(val[301,death.compartments]))
        resultado <- as.numeric(c(PROD.RATE,VAX.WINDOW.DAYS,REL.EFFIC,deaths))
        resultado$Vaccine <- "CoronaVac"
        values <- rbind(values,resultado)
      }
    }
  }
}
for (PROD.RATE in LIST.PRODUCTION.RATE) {
  for(VAX.WINDOW.DAYS in LIST.VAX.WINDOW.DAYS){
    for(REL.EFFIC in LIST.FIRST.DOSE.REL.EFFIC){
      filename <- paste("./results/astrazeneca/sol_AstraZeneca_",PROD.RATE,"_",VAX.WINDOW.DAYS,"_1.1_",REL.EFFIC,".csv",sep = "")
      val <- try(read.table(filename))
      if(class(val) != "try-error"){
        deaths <- as.numeric(sum(val[301,death.compartments]))
        resultado <- as.numeric(c(PROD.RATE,VAX.WINDOW.DAYS,REL.EFFIC,deaths))
        resultado$Vaccine <- "AZD1222"
        values <- rbind(values,resultado)
      }
    }
  }
}
values <- as.data.frame(values)
rownames(values) <- NULL
colnames(values) <- c("PRODUCTION.RATE","VAX.WINDOW.DAYS","REL.EFFIC","DEATHS","VACCINE")

filtered_values <-c()

vaccines <- c("CoronaVac","AZD1222","BNT162b2")
for (PROD.RATE in LIST.PRODUCTION.RATE) {
  for(REL_EFFIC in LIST.FIRST.DOSE.REL.EFFIC){
    for(VAC in vaccines){
      val <- filter(values,PRODUCTION.RATE == PROD.RATE,REL.EFFIC==REL_EFFIC, VACCINE == VAC)
      min_val <- val[which.min(val$DEATHS),]
      filtered_values <- rbind(filtered_values,min_val)
    }
  }
}
rownames(filtered_values) <- NULL
# grid_plot <- ggplot(filtered_values,aes(x=PRODUCTION.RATE,y=REL.EFFIC,col=factor(VAX.WINDOW.DAYS),fill=factor(VAX.WINDOW.DAYS)))+geom_tile()
# 
# data.loess <- loess(VAX.WINDOW.DAYS ~ REL.EFFIC*PRODUCTION.RATE,data=filtered_values,
#                     span=0.18,control=loess.control(surface = c( "interpolate"),
#                     statistics = c("approximate"),
#                     trace.hat = c("approximate"),
#                     cell = 0.2, iterations = 4, iterTrace = FALSE))
# ygrid <- seq(min(filtered_values$REL.EFFIC),max(filtered_values$REL.EFFIC),length.out = 200)
# xgrid <- seq(min(filtered_values$PRODUCTION.RATE),max(filtered_values$PRODUCTION.RATE),length.out = 200)
# data.fit <- expand.grid( PRODUCTION.RATE = xgrid,REL.EFFIC = ygrid)
# mtrx3d <- predict(data.loess,newdata = data.fit)
# contour(x = xgrid, y = ygrid, z = mtrx3d)

# mtrx.melt <- melt(mtrx3d,id.vars=c("PRODUCTION.RATE","REL.EFFIC"),measure.vars="VAX.WINDOW.DAYS")
# mtrx.melt$REL.EFFIC <- as.numeric(str_sub(mtrx.melt$REL.EFFIC, str_locate(mtrx.melt$REL.EFFIC, "=")[1,1] + 1))
# mtrx.melt$PRODUCTION.RATE <- as.numeric(str_sub(mtrx.melt$PRODUCTION.RATE, str_locate(mtrx.melt$PRODUCTION.RATE, "=")[1,1] + 1))
# names(mtrx.melt)<- c("PRODUCTION.RATE","REL.EFFIC","VAX.WINDOW.DAYS")
# plot1 <- ggplot(mtrx.melt, aes(y = REL.EFFIC, x = PRODUCTION.RATE, z = VAX.WINDOW.DAYS)) +
#   stat_contour()
# plot2 <- ggplot(mtrx.melt, aes(x = PRODUCTION.RATE, y = REL.EFFIC, z = VAX.WINDOW.DAYS)) +
#   stat_contour(geom = "polygon", aes(fill = ..level..,color = VAX.WINDOW.DAYS)) +
#   geom_tile(aes(fill = VAX.WINDOW.DAYS)) +
#   stat_contour(breaks = LIST.VAX.WINDOW.DAYS,bins = 22)+geom_label_contour(aes(z=VAX.WINDOW.DAYS),breaks=LIST.VAX.WINDOW.DAYS,skip=0,size=2.5)
# 
# plot4 <- ggplot(filter(filtered_values,PRODUCTION.RATE < 205000), aes(x = PRODUCTION.RATE, y = REL.EFFIC, z = factor(VAX.WINDOW.DAYS/7))) +
#   geom_tile(aes(fill = factor(VAX.WINDOW.DAYS/7)))+labs(y = "Relative Efficacy", x = "Production Rate", fill = "Weeks")
# mtrx.days <- mtrx.melt
# mtrx.days$VAX.WINDOW.DAYS <- round(mtrx.days$VAX.WINDOW.DAYS/7)
# plot5 <- ggplot(mtrx.days, aes(x = PRODUCTION.RATE, y = REL.EFFIC, z = factor(VAX.WINDOW.DAYS)))+geom_tile(aes(fill = factor(VAX.WINDOW.DAYS)))+labs(fill="Weeks")
# plot5
res.list <- meshgrid(seq(LIST.PRODUCTION.RATE[1],LIST.PRODUCTION.RATE[14],length.out = 300),seq(0,1,length.out = 300))
res.list$X <- as.vector(res.list$X)
res.list$Y <- as.vector(res.list$Y)
df <- filtered_values %>% filter(VACCINE == "CoronaVac") %>% select(names(filtered_values)[1:3])
df <- df[,c(1,3,2)]
df <- sapply(df,as.numeric)
df <- as.data.frame(df)
df.2 <- df %>% pivot_wider(names_from = PRODUCTION.RATE,values_from = VAX.WINDOW.DAYS) %>% as.data.frame()
df.2 <- df.2[,2:ncol(df.2)]
# df$VAX.WINDOW.DAYS <- as.numeric(df$VAX.WINDOW.DAYS)
# df.2 <- reshape(data = df,timevar = "PRODUCTION.RATE",idvar = "REL.EFFIC",direction = "wide",sep="_")
# df <- df[,!names(df) %in% c("REL.EFFIC")]
df.2 <- as.matrix(df.2)
# df.2 <- sapply(df.2,as.numeric)

# class(df.2) <- "numeric"
df.2 <- interp2(LIST.PRODUCTION.RATE,LIST.FIRST.DOSE.REL.EFFIC,df.2,res.list$X,res.list$Y)
df.2 <- data.frame(PRODUCTION.RATE = res.list$X, REL.EFFIC = res.list$Y, VAX.WINDOW.DAYS = df.2)
df.2$Weeks <- round(df.2$VAX.WINDOW.DAYS/7)
df.2$Vaccine <- "CoronaVac"
# mtrx.melt <- melt(df,id.vars=c("PRODUCTION.RATE","REL.EFFIC"),measure.vars="VAX.WINDOW.DAYS")
# mtrx.melt$REL.EFFIC <- as.numeric(str_sub(mtrx.melt$REL.EFFIC, str_locate(mtrx.melt$REL.EFFIC, "=")[1,1] + 1))
# mtrx.melt$PRODUCTION.RATE <- as.numeric(str_sub(mtrx.melt$PRODUCTION.RATE, str_locate(mtrx.melt$PRODUCTION.RATE, "=")[1,1] + 1))
# names(mtrx.melt)<- c("PRODUCTION.RATE","REL.EFFIC","VAX.WINDOW.DAYS")




ag.df <- c()
for(VAC in vaccines){
  df <- filtered_values %>% filter(VACCINE == VAC) %>% select(names(filtered_values)[1:3])
  df <- df[,c(1,3,2)]
  df <- sapply(df,as.numeric)
  df <- as.data.frame(df)
  df.2 <- df %>% pivot_wider(names_from = PRODUCTION.RATE,values_from = VAX.WINDOW.DAYS) %>% as.data.frame()
  df.2 <- df.2[,2:ncol(df.2)]
  # df$VAX.WINDOW.DAYS <- as.numeric(df$VAX.WINDOW.DAYS)
  # df.2 <- reshape(data = df,timevar = "PRODUCTION.RATE",idvar = "REL.EFFIC",direction = "wide",sep="_")
  # df <- df[,!names(df) %in% c("REL.EFFIC")]
  df.2 <- as.matrix(df.2)
  # df.2 <- sapply(df.2,as.numeric)
  
  # class(df.2) <- "numeric"
  df.2 <- interp2(LIST.PRODUCTION.RATE,LIST.FIRST.DOSE.REL.EFFIC,df.2,res.list$X,res.list$Y)
  df.2 <- data.frame(PRODUCTION.RATE = res.list$X, REL.EFFIC = res.list$Y, VAX.WINDOW.DAYS = df.2)
  df.2$Weeks <- round(df.2$VAX.WINDOW.DAYS/7)
  df.2$Vaccine <-VAC
  ag.df <- rbind(ag.df,df.2)
  
}
# df.1 <- filtered_values[,c(3,1,2)]
# df.1 <- reshape(data = df.1,timevar = "PRODUCTION.RATE",idvar = "REL.EFFIC",direction = "wide")
# df.1 <- df.1[,!names(df.1) %in% c("REL.EFFIC")]
# df.1 <- as.matrix(df.1)

# df.2 <- interp2(LIST.PRODUCTION.RATE,LIST.FIRST.DOSE.REL.EFFIC,df.1,res.list$X,res.list$Y)
# df.3 <- data.frame(PRODUCTION.RATE = res.list$X, REL.EFFIC = res.list$Y, VAX.WINDOW.DAYS = df.2)
# df.3$Weeks <- round(df.3$VAX.WINDOW.DAYS/7)
# df.3 <- filter(df.3,PRODUCTION.RATE > 100000)
df <- filter(ag.df,PRODUCTION.RATE < 200000)
df$Vaccine <- factor(df$Vaccine,levels = c("CoronaVac","AZD1222","BNT162b2"))
plot6 <- ggplot(df, aes(x = PRODUCTION.RATE, y = REL.EFFIC, z = factor(Weeks),group = Vaccine))+facet_wrap(~Vaccine,nrow = 1)+geom_raster(aes(fill = factor(Weeks)))+
        labs(x = "Production Rate", y = "Relative Efficacy",fill="Weeks")


