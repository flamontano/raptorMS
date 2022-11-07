## Plotting global patterns of observed and SES Functional and phylogenetic diversity of raptors ##


### OBSERVED VALUES ###

data_raptor.plots<-data_raptor[data_raptor$continent !="Antarctica",]

# grouping some values to improve contrast in visualization #

data_raptor.plots$mpd1<-data_raptor.plots$mpd
data_raptor.plots$mpd1[data_raptor.plots$mpd1 <75] <- 74.9

data_raptor.plots$mntd1<-data_raptor.plots$mntd
data_raptor.plots$mntd1[data_raptor.plots$mntd1 >95] <-96
data_raptor.plots$mntd1[data_raptor.plots$mntd1 < 12] <-11.9

data_raptor.plots$fdis_obs_niche2<-data_raptor.plots$fdis_obs_niche
data_raptor.plots$fdis_obs_niche2[data_raptor.plots$fdis_obs_niche2 < 0.2] <-0.15


library(ggplot2)

p.0 = plot_world_eqaul_area(color_polygon = "grey90") +
  viridis::scale_fill_viridis(direction = -1) +
  theme(legend.position = c(0.55,0.10),
        legend.direction = "horizontal", 
        legend.key.height = unit(0.5, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        plot.margin = margin(-0.5, -0.1, -0.5, -0.1, "cm"))

p_sprich = p.0 + geom_tile(data = filter(data_raptor.plots, !is.na(continent)),
                         aes(x = x, y = y, fill = sr), inherit.aes = F) +
   labs(fill = 'SR')
p_sprich


p = plot_world_eqaul_area(color_polygon = "grey90", fill_polygon = "white") +
  viridis::scale_fill_viridis(direction = -1) +
  theme(legend.position = c(0.55,0),
        legend.direction = "horizontal", 
        legend.key.height = unit(0.4, 'lines'),
        legend.key.width = unit(1.5, 'lines'),
        plot.margin = margin(-0.5, -0.1, -0.5, -0.1, "cm"))

# plot pd uroot
p_pduroot = p + geom_tile(data = data_raptor.plots,
                          aes(x = x, y = y, fill = pd.uroot),
                          inherit.aes = F) +
   labs(fill = 'PD')

p_mpd = p + geom_tile(data = data_raptor.plots,
                          aes(x = x, y = y, fill = mpd1),
                          inherit.aes = F) +
  labs(fill = 'MPD  ')+ # Change legend labels of continuous legend
  scale_fill_continuous(type = "viridis", direction= -1, breaks = c(80, 100, 120, 140, 160), labels = c("< 80", "100", "120", "140", "160"))

p_mntd = p + geom_tile(data = data_raptor.plots,
                          aes(x = x, y = y, fill = mntd1),
                          inherit.aes = F) +
  labs(fill = 'MNTD  ')+ # Change legend labels of continuous legend
  scale_fill_continuous(type = "viridis", direction= -1, breaks = c(10, 30, 50, 70, 90), labels = c("<10", "30", "50", "70", ">90"))


#plot fd

p_morph = p + geom_tile(data = data_raptor.plots,
                          aes(x = x, y = y, fill = fdis_obs_morph),
                          inherit.aes = F) +
   labs(fill = 'FD Morphology')

p_fdis_niche = p + geom_tile(data = data_raptor.plots,
                               aes(x = x, y = y, fill = fdis_obs_niche2),
                               inherit.aes = F) +
   labs(fill = 'FD Niche ')+
  scale_fill_continuous(type = "viridis", direction= -1, breaks = c(0, 0.15, 0.3, 0.45), labels = c("0.0", "0.15", "0.30", "0.45"))


p_fdis_diet = p + geom_tile(data = data_raptor.plots,
                             aes(x = x, y = y, fill = fdis_obs_diet),
                             inherit.aes = F) +
  labs(fill = 'FD Diet ')+ # Change legend labels of continuous legend
  scale_fill_continuous(type = "viridis", direction= -1, breaks = c(0, 0.03, 0.06, 0.09), labels = c("0.0", "0.03", "0.06", "0.09"))

p_fdis_forag = p + geom_tile(data = data_raptor.plots,
                             aes(x = x, y = y, fill = fdis_obs_forag),
                             inherit.aes = F) +
  labs(fill = 'FD Foraging  ')

p_fdis_dispersal = p + geom_tile(data = data_raptor.plots,
                               aes(x = x, y = y, fill = fdis_obs_dispe),
                               inherit.aes = F) +
   labs(fill = 'FD Vagility ')


## SES VALUES ###

pz = plot_world_eqaul_area(color_polygon = "white") +
   scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC") +
   theme(legend.position = c(0.55, 0),
         legend.direction = "horizontal", 
         legend.key.height = unit(0.5, 'lines'),
         legend.key.width = unit(1.5, 'lines'),
         plot.margin = margin(-0.5, -0.1, -0.5, -0.1, "cm"))

#phylogenetic

p_pduroot_z = pz + geom_tile(data = data_raptor.plots,
                             aes(x = x, y = y, fill = pd.uroot.z),
                             inherit.aes = F) +
   labs(fill = 'SES PD')
p_mpd_z = pz + geom_tile(data = data_raptor.plots,
                             aes(x = x, y = y, fill = mpd.z),
                             inherit.aes = F) +
  labs(fill = 'SES MPD')
p_mntd_z = pz + geom_tile(data = data_raptor.plots,
                             aes(x = x, y = y, fill = mntd.z),
                             inherit.aes = F) +
  labs(fill = 'SES MNTD')

p_morph_z = pz + geom_tile(data = data_raptor.plots,
                                  aes(x = x, y = y, fill = fdis_ses_morph),
                                  inherit.aes = F) +
   labs(fill = 'SES FD Morphology')

p_fdis_niche_z = pz + geom_tile(data = data_raptor.plots,
                               aes(x = x, y = y, fill = fdis_ses_niche),
                               inherit.aes = F) +
   labs(fill = 'SES FD Niche')

p_fdis_diet_z = pz + geom_tile(data = data_raptor.plots,
                                aes(x = x, y = y, fill = fdis_ses_diet),
                                inherit.aes = F) +
  labs(fill = 'SES FD Diet')
p_fdis_forag_z = pz + geom_tile(data = data_raptor.plots,
                                aes(x = x, y = y, fill = fdis_ses_forag),
                                inherit.aes = F) +
  labs(fill = 'SES FD Foraging ')

p_fdis_dispersal_z = pz + geom_tile(data = data_raptor.plots,
                                 aes(x = x, y = y, fill = fdis_ses_dispe),
                                 inherit.aes = F) +
   labs(fill = 'SES FD Vagility')


## FINAL FIGURES ##

library(cowplot)

p_obs = plot_grid(p_pduroot, 
                  p_mpd, 
                  p_mntd, 
                  p_morph, 
                  p_fdis_niche, 
                  p_fdis_diet, 
                  p_fdis_forag, 
                  p_fdis_dispersal, 
                  ncol = 2, labels = letters[2:9])

p_obs2 = plot_grid(p_sprich, p_obs, labels = c('a', ''), ncol = 1, rel_heights = c(0.3, 1))

#ggsave(filename = "Figures/newfig1_obs.pdf", plot = p_obs2, height = 15, width = 10)

p_ses = plot_grid(p_pduroot_z,
                  p_mpd_z, 
                  p_mntd_z,
                  p_morph_z,
                  p_fdis_niche_z,
                  p_fdis_diet_z,
                  p_fdis_forag_z,
                  p_fdis_dispersal_z,
                  ncol = 2, labels = letters[1:8])


#ggsave(filename = "Figures/newfig2_ses_alt.pdf", plot = p_ses, height = 15, width = 11)


