source('./functions/dist_age_states.R')
source('./functions/ihfr_estados.R')

## this file updates the data extracted from SIVEP database and saves to the
## DATA folder

# recent cases
recent <- dist_etaria_estados(data_ini = "2021-01-17", data_fim = "2021-02-13",
                              more_age_bins=F)
write.csv(recent, "./DATA/recent_cases_dist_age_states.csv", row.names=F)

# all cases
cases <- dist_etaria_estados(data_fim = "2021-02-06", more_age_bins=T)
write.csv(cases, "./DATA/cases_dist_age_states.csv", row.names=F)

# all deaths
deaths <- dist_etaria_estados(data_fim = "2021-02-06", more_age_bins=T, dist_obs=T)
write.csv(cases, "./DATA/deaths_dist_age_states.csv", row.names=F)

# IHFR
ihfr <- ihfr_estados()
write.csv(ihfr, "./DATA/ihfr_states.csv", row.names=F)
