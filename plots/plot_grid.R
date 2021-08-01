library(ggplot2)


plot6 <- ggplot(vaccine_grid, aes(x = PRODUCTION.RATE/POP.TOTAL.NUM, y = REL.EFFIC, z = factor(Weeks)))+
  facet_wrap(~Vaccine,nrow = 1)+geom_raster(aes(fill = factor(Weeks)),interpolate = TRUE)+
  labs(x = "Production Rate (doses/population-day)", y = "Relative Efficacy of First Dose",fill="Weeks")+
  theme_bw()+theme(text = element_text(size=16))+scale_fill_viridis_d()+scale_x_continuous(labels = scales::percent_format(accuracy = 0.1))

