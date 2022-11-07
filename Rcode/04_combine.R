# merge sr, pd, fd together

sp_div_df = left_join(mutate(coord_continents, x = as.character(x), y = as.character(y)),
                      select(pd_rapt, -rowid), 
                      by = c("x", "y")) 

sp_div_df = full_join(mutate(rap_fd, x = as.character(x), y = as.character(y)), 
                      sp_div_df,
                      by = c("x", "y", "country", "continent"))

# environmental data: includes elevation SD, ave_temp_current, ave_precip_current, ave_temp_paleo, ave_precip_paleo, ave_NPP

envi = readRDS("Data/output/envi.rds")

data_raptor = left_join(sp_div_df, mutate(envi, x = as.character(x), y = as.character(y)),
                        by = c("x", "y"))
data_raptor = mutate(data_raptor, x = as.numeric(x), y = as.numeric(y)) %>% 
  select(x, y, lon_wgs84, lat_wgs84, lon_eck4, lat_eck4, sr, everything())


data_raptor = filter(data_raptor, !is.na(sr), !is.na(continent))



