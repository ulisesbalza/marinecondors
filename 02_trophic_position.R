##### Formatting data for Trophic position estimation #####
library(tRophicPosition)
library(ggplot2)
library(ggpubr)
library(cowplot)

# example of how to produce baseline data

# Riccialdelli et al 2017 https://link.springer.com/article/10.1007/s00300-016-2007-x (For Fuego)
# mytilus_Riccialdelli_N <- rnorm(n=18, mean=11.9, sd=0.7)
# mytilus_Riccialdelli_C <- rnorm(n=18, mean=-18.7, sd=0.8)
# marine_baseline_Riccialdelli <- data.frame(mytilus_Riccialdelli_C, 
#                            mytilus_Riccialdelli_N)

#write.csv(marine_baseline_Riccialdelli, "data/marine_baseline_lapataia.csv")

# loading databases
data_fuego_all <- read.csv("data/condor_fuego_all.csv", sep = ",", header = T)
data_fuego_clusters <- read.csv("data/condor_fuego_by_clusters.csv", sep = ",", header = T)

data_patagonia <- read.csv("data/condor_patagonia.csv", sep = ",", header = T)
data_puna <- read.csv("data/condor_san_guillermo.csv", sep = ",", header = T)
data_auca <- read.csv("data/condor_auca.csv", sep = ",", header = T)
data_cordoba <- read.csv("data/condor_cordoba.csv", sep = ",", header = T)
data_payunia <- read.csv("data/condor_payunia.csv", sep = ",", header = T)
data_quijadas <- read.csv("data/condor_quijadas.csv", sep = ",", header = T)

# Defining consumers and baselines

# Generating raw numbers for chosen Trophic discrimination factor 
# Data from Kurle et al (2013) https://academic.oup.com/condor/article-abstract/115/3/492/5152863 
kurle_N = rnorm(n=30, mean=3.1, sd=0.2)
kurle_C = rnorm(n=30, mean=0.4, sd=0.4)

# FUEGIAN ALL
consumer_fuego <- extractIsotopeData(data_fuego_all, 
                                     consumersColumn = "FG",
                                     b1 = "terrestrial",
                                     b2 = "marine", 
                                     baselineColumn = "FG",
                                     deltaC= kurle_C,
                                     deltaN= kurle_N,
                                     groupsColumn = NULL)

# FUEGIAN BY CLUSTERS
consumer_fuego_clusters <- extractIsotopeData(data_fuego_clusters, 
                                              consumersColumn = "FG",
                                              b1 = "terrestrial",
                                              b2 = "marine", 
                                              baselineColumn = "FG",
                                              deltaC= kurle_C,
                                              deltaN= kurle_N,
                                              groupsColumn = NULL)


# PATAGONIA
consumer_patagonia <- extractIsotopeData(data_patagonia, 
                                         consumersColumn = "FG",
                                         b1 = "terrestrial",
                                         b2 = "marine", 
                                         baselineColumn = "FG",
                                         deltaC= kurle_C,
                                         deltaN= kurle_N,
                                         groupsColumn = NULL)



# SAN GUILLERMO
consumer_puna <- extractIsotopeData(data_puna, 
                                    consumersColumn = "FG",
                                    b1 = "terrestrial",
                                    baselineColumn = "FG",
                                    deltaC= kurle_C,
                                    deltaN= kurle_N,
                                    groupsColumn = NULL)

# AUCA MAHUIDA
consumer_auca <- extractIsotopeData(data_auca, 
                                    consumersColumn = "FG",
                                    b1 = "terrestrial",
                                    baselineColumn = "FG",
                                    deltaC= kurle_C,
                                    deltaN= kurle_N,
                                    groupsColumn = NULL)

# CORDOBA
consumer_cordoba <- extractIsotopeData(data_cordoba, 
                                       consumersColumn = "FG",
                                       b1 = "terrestrial",
                                       baselineColumn = "FG",
                                       deltaC= kurle_C,
                                       deltaN= kurle_N,
                                       groupsColumn = NULL)

# PAYUNIA
consumer_payunia <- extractIsotopeData(data_payunia, 
                                       consumersColumn = "FG",
                                       b1 = "terrestrial",
                                       baselineColumn = "FG",
                                       deltaC= kurle_C,
                                       deltaN= kurle_N,
                                       groupsColumn = NULL)

# QUIJADAS
consumer_quijadas <- extractIsotopeData(data_quijadas, 
                                        consumersColumn = "FG",
                                        b1 = "terrestrial",
                                        baselineColumn = "FG",
                                        deltaC= kurle_C,
                                        deltaN= kurle_N,
                                        groupsColumn = NULL)

#Example of the data that is going to be used for the TP estimations
# For instance, the ones depicted in Figure 2
plot(consumer_fuego$condor,
     b1="Terrestrial baseline",
     b2="Marine baseline",
     xylim= c(-30, -13, 0, 22))
plot(consumer_patagonia$condor,
     b1="Terrestrial baseline",
     b2="Marine baseline",
     xylim= c(-30, -13, 0, 22))

#### Trophic position multispecies model ####

# alpha = relative contribution of baseline 1 (terrestrial)
# lambda: baseline TP

# FUEGO ALL
condor_fuego <- multiSpeciesTP(consumer_fuego, model = "twoBaselinesFull",
                               lambda = 2, n.adapt = 50000, n.iter = 50000,
                               burnin = 50000, n.chains = 5, print = TRUE)

# FUEGO BY CLUSTERS
condor_fuego_clusters <- multiSpeciesTP(consumer_fuego_clusters, model = "twoBaselinesFull",
                               lambda = 2, n.adapt = 50000, n.iter = 50000,
                               burnin = 50000, n.chains = 5, print = TRUE)

# compare distributions of Fuego by clusters


# NORTHERN PATAGONIA
condor_patagonia <- multiSpeciesTP(consumer_patagonia, model = "twoBaselinesFull",
                                      lambda = 2, n.adapt = 50000, n.iter = 50000,
                                      burnin = 50000, n.chains = 5, print = TRUE)


# PUNA
condor_puna <- multiSpeciesTP(consumer_puna, model = "oneBaseline",
                               lambda = 2, n.adapt = 50000, n.iter = 50000,
                               burnin = 50000, n.chains = 5, print = TRUE)

# AUCA MAHUIDA
condor_auca <- multiSpeciesTP(consumer_auca, model = "oneBaseline",
                              lambda = 2, n.adapt = 50000, n.iter = 50000,
                              burnin = 50000, n.chains = 5, print = TRUE)

# CORDOBA
condor_cordoba <- multiSpeciesTP(consumer_cordoba, model = "oneBaseline",
                              lambda = 2, n.adapt = 50000, n.iter = 50000,
                              burnin = 50000, n.chains = 5, print = TRUE)

# PAYUNIA
condor_payunia <- multiSpeciesTP(consumer_payunia, model = "oneBaseline",
                              lambda = 2, n.adapt = 50000, n.iter = 50000,
                              burnin = 50000, n.chains = 5, print = TRUE)

# QUIJADAS
condor_quijadas <- multiSpeciesTP(consumer_quijadas, model = "oneBaseline",
                              lambda = 2, n.adapt = 50000, n.iter = 50000,
                              burnin = 50000, n.chains = 5, print = TRUE)


# Generate Figure 5

df_plot <- data.frame(site= c("San Guillermo", "Córdoba", "Quijadas",
                              "Payunia", "Auca Mahuida", "Patagonia", 
                              "Tierra del Fuego  
(terrestrial)", "Tierra del Fuego 
(marine)"),
                      tp_median= c(3.18,2.96,3.04,
                                   3.66,3.47,3.53,
                                   2.24, 3.27),
                      tp_lci= c(2.89,2.59,2.68,
                                3.46,3.21,3.35,
                                2.01,2.72),
                      tp_hci= c(3.47,3.33,3.39,
                                3.86,3.73,3.72,
                                2.91,3.85))

# plot
trophic_position <-  ggplot(df_plot, aes(x=reorder(site,tp_median), y=tp_median)) +
  geom_point(stat="identity", size=8, shape=15,
             color=c("#765E55", "#765E55", "#765E55",  "#765E55",
                     "#765E55", "#547AA5", "#765E55", "#547AA5")) +
  geom_errorbar(aes(ymin=tp_lci, ymax=tp_hci), 
                col="black", width=.2) +
  labs(title="", 
       subtitle = "",
       x="", y = "Estimated trophic position") +
  #  draw_image("data/vultur_fly.png", x = 0.6, y = 3.3, scale = .9) +
  geom_hline(yintercept=3, linetype="dashed", 
             color = "black", linewidth=0.8) +
  theme_classic(base_size = 14) +
  labs(color = "Marine input") +
  theme(legend.title = element_text(color = "black", size = 12),
        legend.position = "right",
        legend.background = element_rect(fill="lightgray", colour="black",
                                         linewidth=0.7, linetype="solid"))

trophic_position


