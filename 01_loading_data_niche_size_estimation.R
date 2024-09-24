##### Maintenance of a scavenger ancient food web in insular remote areas ####
### Loading isotopic data and preparing it for the modeling 

# Loading libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# descriptive data condors from Tierra del Fuego
condor_tdf_all <- read.csv("data/condor_tdf_todos.csv", sep = ",", header = T)

# summary of data
condor_tdf_all %>% group_by(Location) %>% 
  summarise(count = length(Location),
            mC = mean(d13C), 
            sdC = sd(d13C), 
            mN = mean(d15N), 
            sdN = sd(d15N))

# ellipses
p.ell <- 0.70 

Figure_3A <-  ggplot(data = condor_tdf_all, aes(d13C, d15N)) +
  geom_point(aes(fill = Location, color= Location), 
              stroke = 1.5, shape = 21,
             size = 4, alpha = 0.6) +
  facet_wrap(. ~ Location,
             strip.position = 'top', 
             labeller = as_labeller(c(IDLE='Isla de los Estados', 
                                      West='Western', 
                                      East='Eastern'))) +
  scale_fill_manual(values = c("IDLE" = "#DDD78D",
                               "West" = "#B87D4B",
                               "East" = "#A4DE82" )) +
  scale_color_manual(values = c("IDLE" = "#DDD78D",
                                "West" = "#B87D4B",
                                "East" = "#A4DE82" )) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  scale_y_continuous(breaks = seq(5, 20, by = 5), limits = c(5, 22)) +
  theme(text = element_text(size=15)) + 
  stat_ellipse(aes(group = Location, 
                   fill = Location, 
                   color = Location), 
               alpha = 0.2, 
               level = p.ell,
               type = "norm",
               geom = "polygon") + 
  labs(fill="Site", color = "Site") +
  theme_light() + 
  theme(text = element_text(size=15), 
        legend.position="none")

Figure_3A

# defining clusters
library(cluster)
library(factoextra)

condor_matrix <- data.frame(condor_tdf_all$d13C, condor_tdf_all$d15N)

dist.eucl <- dist(condor_matrix, method = "euclidean")

# optimal number of clusters
gap.stat <- clusGap(condor_matrix, FUNcluster = kmeans, K.max = 5)
fviz_gap_stat(gap.stat)

res.agnes <- agnes(x = condor_matrix,      
                   stand = FALSE,   
                   metric = "euclidean", 
                   method = "average"
)

fviz_dend(res.agnes, cex = 0.6, k = 2) # assing individuals to two clusters

grp <- cutree(res.agnes, k = 2)

condor_tdf_all$cluster <- grp
condor_tdf_all$cluster <- as.factor(condor_tdf_all$cluster) # add cluster to each individual
head(condor_tdf_all)

# clusters plot
Figure_3B <- ggplot(data = condor_tdf_all, aes(d13C, d15N)) +
  geom_point(aes(fill=cluster, color = cluster), stroke = 1.5, shape = 21,
             size = 4, alpha = 0.6) +
  facet_wrap(.~ cluster,
             strip.position = 'top', 
             labeller = as_labeller(c('1'='Marine (n=12)', 
                                      '2'='Terrestrial (n=31)'))) +
  scale_fill_manual(values = c("1" = "#547AA5",
                               "2" = "#765E55")) +
  scale_color_manual(values = c("1" = "#547AA5",
                                "2" = "#765E55")) +
  ylab(expression(paste(delta^{15}, "N (\u2030)")))+
  xlab(expression(paste(delta^{13}, "C (\u2030)"))) + 
  scale_y_continuous(breaks = seq(5, 20, by = 5), limits = c(5, 22)) +
  stat_ellipse(aes(group = cluster, 
                   fill = cluster, 
                   color = cluster), 
               alpha = 0.2, 
               type = "norm",
               geom = "polygon",
               show.legend = FALSE) + 
  theme_light() + 
  theme(text = element_text(size=15),
        legend.position="none")

Figure_3B


# describing clusters
library(compareGroups)
group <- compareGroups(Location ~ cluster, data=condor_tdf_all)
clustab <- createTable(group,digits=3, 
                       show.p.overall=FALSE)
clustab


### NicheROVER to estimate isotopic niche width of each population ####
library(nicheROVER)

condor_all <- read.csv("data/condor_all.csv", sep = ",", header = T)

table(condor_all$site) #samples per site

# 2-d projections of 10 niche regions
clrs <- c("black", "red", "blue", "orange", "purple", "green", "lightblue") # colors for each species
nsamples <- 2000
par <- tapply(1:nrow(condor_all), condor_all$site,
                   function(ii) niw.post(nsamples = nsamples, X = condor_all[ii,3:4]))

# format data for plotting function
data <- tapply(1:nrow(condor_all), condor_all$site, function(ii) X = condor_all[ii,3:4])

niche.plot(niche.par = par, niche.data = data, pfrac = .05,
           iso.names = expression(delta^{15}*N, delta^{13}*C, delta^{34}*S),
           col = clrs, xlab = expression("Isotope Ratio (per mil)"))

# posterior distribution of (mu, Sigma) for each group
nsamples <- 2000

# posterior distribution of niche size by group
size <- sapply(par, function(spec) {
  apply(spec$Sigma, 3, niche.size, alpha = .95)
})

# point estimate and standard error
rbind(est = colMeans(size),
      se = apply(size, 2, sd))

# boxplots
boxplot(size, col = "grey", pch = 16, cex = .5,
        ylab = "Isotopic Niche Size", xlab = "Population")

# plot with ggplot
niche <- stack(as.data.frame(size))
head(niche)

# Figure 4 ####

ggplot(niche, aes(x=reorder(ind, values), y=values)) + 
  geom_boxplot(width = .30, outlier.color = NA, outlier.shape = NA,
               fill= c("#547AA5", "#765E55", "#765E55",
                               "#765E55", "#765E55", "#765E55", "#547AA5")) +
  ggdist::stat_halfeye(adjust = .5, width = 1, justification = -.2, 
                       .width = 0, 
                       point_colour = NA,
                       fill = "gray") + 
  labs(title="", 
       subtitle = "",
       x="", y = expression("Isotopic Niche Size (‰)"^2)) +
  coord_cartesian(ylim = c(0,60)) +
  theme_classic(base_size = 12) 

