# Purpose: Statistical analysis changes in truffle-like fungal community composition between seasons and with changes in climate 

# Load packages 


library(phyloseq)
library(ggplot2)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(tibble)
library(indicspecies)
library(ggrepel)
library(pscl)
library(emmeans)
library(MASS)

# 1. Ordinations  ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data=prune_taxa(ecm.list, ITSrel.count2)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")] # Filter for truffle-like genera
data=prune_taxa(hypo.list, data)

data = prune_taxa(taxa_sums(data) > 0, data) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data) > 0, data) # Remove samples with zero OTUs 


iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

data.scores <- as.data.frame(scores(iMDS))
sort(data.scores$NMDS1) 

truffle_ecm <- subset_samples(truffle_ecm, !agrf. %in% c("1", "18","286")) # Remove outliers that obscure patterns in ordination

iDist <- distance((truffle_ecm), method="bray")
iMDS  <- ordinate(truffle_ecm, "NMDS", distance=iDist)

# Ordinate samples by season
ordination_df <- as.data.frame(scores(iMDS, display = "sites"))
ord_truffle_ecm <- cbind(ordination_df, sample_data(truffle_ecm))

ord_truffle_ecm$Season <- factor(ord_truffle_ecm$Season,levels = c("Summer", "Autumn", "Winter", "Spring"))

season_colors <- c("Summer"= "#EAA48F","Autumn" = "#D1B08F","Winter" = "#92c9c4","Spring" = "#95A58E") 

ordplot <- ggplot(ord_truffle_ecm, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3) +
  scale_colour_manual(values = season_colors)+
  theme_light() +
  theme(aspect.ratio = 1) +
  ggtitle("Seasonal Differences in Community Composition of Truffle-like ECM Fungi") +
  theme(plot.title = element_text(size = 16))
ordplot

find_hull <- function(df) df[chull(df$NMDS1, df$NMDS2), ]

# Calculate the polygons for each season
hulls <- ord_truffle_ecm %>%
  group_by(Season) %>%
  do(find_hull(.))

# Add polygons to ordination plot
ordplot + 
  geom_polygon(data = hulls, aes(x = NMDS1, y = NMDS2, fill = Season, group = Season), alpha = 0.05)  +
  theme(plot.title = element_text(size = 13)) +scale_fill_manual(values = season_colors) + theme_light()



# 2. PERMANOVAs  ----

sampledf <- data.frame(sample_data(truffle_ecm))
otu_table <- as(t(otu_table(truffle_ecm)), "matrix")
bray_dist <- vegdist(otu_table, method = "bray")

# PERMANOVA to test for differences in community composition between seasons 
permanova_result <- adonis2(bray_dist ~ Season, data = sampledf)
print(permanova_result)  

# Pairwise differences in community composition between seasons 
pairwise.adonis(bray_dist, sampledf$Season)

# Test for difference in group dispersion between seasons
beta_season <- betadisper(bray_dist, sampledf$Season)
permutest(beta_season) 
plot(beta_season, hull=FALSE, ellipse=TRUE) 
boxplot(beta_season)




# 3. Indicspecies ----

# Test for genera that are indicative of particular seasons 

sampledf <- data.frame(sample_data(truffle_ecm))
data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus

otu <- as.data.frame(t(otu_table(data_genus)))
tax <- as.data.frame(t(tax_table(data_genus))) %>% rownames_to_column(var = "Taxa")
genus <- tax %>% filter(Taxa == "Genus")
genus <- genus[-1]
genus_otu <- as.data.frame(rbind(genus, otu))
colnames(genus_otu) <- genus_otu[1,]
genus_otu <- genus_otu[-1,]
genus_otu <- genus_otu %>%
  mutate(across(everything(), as.numeric))

genus_otu_bind <- cbind(sampledf, genus_otu)
Season <- (sample_data(genus_otu_bind)$Season)


indicator_taxa <- multipatt(genus_otu, Season, func = "r", 
                 control = how(nperm=999)) 
summary(indicator_taxa)


# 4. DB-RDA ----

data_genus <- tax_glom(truffle_ecm, taxrank = "Genus") # Group taxa at the level of genus
otu_table <- as(t(otu_table(data_genus)), "matrix") 
sample_data <- data.frame(sample_data(data_genus))

col_names_indices <- data.frame(Index = seq_along(colnames(sample_data)), Column = colnames(sample_data))
print(col_names_indices)

sample_data[9:28] <- lapply(sample_data[9:28], function(x) as.numeric(as.character(x)))
sample_data[9:28]  <- scale(sample_data[9:28]) # Scale climate variables


bray_distance <- vegdist(otu_table, method = "bray")

#use ordistep to chose model by permutation tests.

null<- capscale(bray_distance ~ 1, data = sample_data, distance = "bray")

# Prepare a model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)
full_model <- capscale(bray_distance ~ Season + week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + fdsi_annual, data = sample_data, distance = "bray")

# Select model
model_select <- ordiR2step(null, scope = formula(full_model), R2scope = TRUE, pstep = 0.05)

formula(model_select)
print(model_select)
summary(model_select)

# Check for VIF
vif_values <- vif.cca(model_select)
print(vif_values)
                            
marginal_effects <- anova(model_select, by = "margin")
marginal_effects

anova_overall <- anova.cca(model_select, permutations = 999)
print(anova_overall)
db_rda_scores <- scores(model_select, display = "sites")
db_rda_env <- scores(model_select, display = "bp")
scores <- as.data.frame(db_rda_env)

# Combine the scores with the sample data
sites_data <- data.frame(db_rda_scores, sample_data) 
env_data <- data.frame(db_rda_env)


season_colors <- c("Summer"= "#EAA48F","Autumn" = "#D1B08F","Winter" = "#92c9c4","Spring" = "#95A58E") 

rownames(env_data)
env_data$name <- c("Spring","Summer","Winter",  "Twelve Month MINT","Weekly AVT", "Weekly Precip")

p <- ggplot(sites_data, aes(x = CAP1, y = CAP2, color = Season)) +
  geom_point(size = 3) 

find_hull <- function(df) df[chull(df$CAP1, df$CAP2), ]

# Calculate polygon for each season
hulls <- sites_data %>%
  group_by(Season) %>%
  do(find_hull(.))

# Add polygons to the plot
p <- p + 
  geom_polygon(data = hulls, aes(x = CAP1, y = CAP2, fill = Season, group = Season), alpha = 0.02) 
p +
  scale_fill_manual(values = season_colors) + scale_color_manual(values = season_colors)+
  geom_segment(data = env_data, aes(x = 0, xend = CAP1, y = 0, yend = CAP2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#000000", size = 0.8) +
  geom_text_repel(data = env_data, aes(x = CAP1, y = CAP2, label = name), 
                  color = "#000000", hjust =0, vjust = 0.0 ,force_pull = -0.06,
                  segment.size = 0.5,  # Controls the thickness of the lines
                  segment.linetype = "dashed", size = 4) + 
  theme_minimal() +
  theme(legend.position = "top") +
  labs(title = "db-RDA of Truffle-like ECM at Bellbird", x = "CAP1", y = "CAP2") 


# 5. Occurrence and relative abundance of truffle-like genera ----

## 5.1 Hysterangium ----


tax_table <- tax_table(ITSrel.count)
hysterangium_taxa <- taxa_names(tax_table[tax_table[, "Genus"] == "Hysterangium", ])
data <- prune_taxa(hysterangium_taxa, ITSrel.count) #

data = prune_taxa(taxa_sums(data) > 0, data) ## remove OTUs with zero count

otu_table_df <- as.data.frame(t(otu_table(data)))
sample_data_df <- as.data.frame(sample_data(data))
tax_table_df <- as.data.frame(tax_table(data))

hyst_abundance <- as.data.frame(rowSums(otu_table_df))
colnames(hyst_abundance)[1] <- "hyst_RelativeAbundance"

hyst_abundance_df <- data.frame(SampleID = rownames(hyst_abundance),hyst_RelativeAbundance = hyst_abundance)

hyst_abundance_data <- merge(sample_data_df, hyst_abundance_df, by.x = "row.names", by.y = "SampleID")
colnames(hyst_abundance)[1] <- "SampleID"

col_names_indices <- data.frame(Index = seq_along(colnames(hyst_abundance_data)), Column = colnames(hyst_abundance_data))
print(col_names_indices)

hyst_abun <- hyst_abundance_data %>% mutate_at(vars(10:29), as.numeric)
hyst_abun[, 10:29] <- scale(hyst_abun[, 10:29])


# Prepare a hurdle model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)

hurdle_hyst_full <- hurdle(hyst_RelativeAbundance ~  week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + wk_maxt_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = hyst_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_hyst_full, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)

# Remove variables that are not significant 

hurdle_hyst_1 <- hurdle(hyst_RelativeAbundance ~   three_mon_mint_av + wk_av_temp  , dist = "negbin", zero.dist = "binomial", data = hyst_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_hyst_1, direction = "both", 
                      trace = FALSE)
summary(step.model)

# Compare to null model

hurdle_hyst_null <- hurdle(hyst_RelativeAbundance ~   1   , dist = "negbin", zero.dist = "binomial", data = hyst_abun, link = "logit", na.action = "na.fail")
AIC(hurdle_hyst_1,hurdle_hyst_null) 


# Plot the relationship with each variable for occurrence and relative abundance, holding the other variable at its mean

hyst_prob0_wk_av_temp <- emmip(hurdle_hyst_1, ~ wk_av_temp,
                           at = list(wk_av_temp = seq(min(hyst_abun$wk_av_temp), max(hyst_abun$wk_av_temp), length.out = 50),
                                     three_mon_mint_av = mean(hyst_abun$three_mon_mint_av)
                           ), 
                           lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 


hyst_prob0_quarterly_mint <- emmip(hurdle_hyst_1, ~ three_mon_mint_av,
                                   at = list(three_mon_mint_av = seq(min(hyst_abun$three_mon_mint_av), max(hyst_abun$three_mon_mint_av), length.out = 50),
                                             wk_av_temp = mean(hyst_abun$wk_av_temp)
                                   ), 
                                   lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 

hyst_abun_quarterly_mint <- emmip(hurdle_hyst_1, ~ three_mon_mint_av,
                                  at = list(three_mon_mint_av = seq(min(hyst_abun$three_mon_mint_av), max(hyst_abun$three_mon_mint_av), length.out = 50),
                                            wk_av_temp = mean(hyst_abun$wk_av_temp)
                                  ), 
                                  lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 

hyst_abun_wk_av_temp <-emmip(hurdle_hyst_1, ~ wk_av_temp,
                             at = list(wk_av_temp = seq(min(hyst_abun$wk_av_temp), max(hyst_abun$wk_av_temp), length.out = 50),
                                       three_mon_mint_av = mean(hyst_abun$three_mon_mint_av)
                             ), 
                             lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 


var1_long <- hyst_abun_quarterly_mint %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)

var2_long <- hyst_abun_wk_av_temp %>% 
  mutate(variable = "wk_av_temp") %>% 
  dplyr::select(variable, yvar, SE, xvar = wk_av_temp)

hyst_combined_count <- bind_rows(var1_long, var2_long)


# Plot relative abundance 
ggplot(hyst_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance", title = "Hysterangium") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() #+
  #theme(legend.position = "none")  


var1_long <- hyst_prob0_quarterly_mint %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)

var2_long <- hyst_prob0_wk_av_temp  %>% 
  mutate(variable = "wk_av_temp") %>% 
  dplyr::select(variable, yvar, SE, xvar = wk_av_temp)

hyst_combined_zero <- bind_rows(var1_long, var2_long)
str(hyst_combined_zero$occurenceprob)

# Plot occurrence
ggplot(hyst_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Hysterangium Presence", title = "Hysterangium") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() #+
  #theme(legend.position = "none")  




## 5.2 Mesophellia ----


tax_table <- tax_table(ITSrel.count)
Mesophellia_taxa <- taxa_names(tax_table[tax_table[, "Genus"] == "Mesophellia", ])
data <- prune_taxa(Mesophellia_taxa, ITSrel.count) 

data = prune_taxa(taxa_sums(data) > 0, data) ## remove OTUs with zero count

otu_table_df <- as.data.frame(t(otu_table(data)))
sample_data_df <- as.data.frame(sample_data(data))
tax_table_df <- as.data.frame(tax_table(data))

meso_abundance <- as.data.frame(rowSums(otu_table_df))
colnames(meso_abundance)[1] <- "meso_RelativeAbundance"

meso_abundance_df <- data.frame(SampleID = rownames(meso_abundance),meso_RelativeAbundance = meso_abundance)

meso_abundance_data <- merge(sample_data_df, meso_abundance_df, by.x = "row.names", by.y = "SampleID")
colnames(meso_abundance)[1] <- "SampleID"

col_names_indices <- data.frame(Index = seq_along(colnames(meso_abundance_data)), Column = colnames(meso_abundance_data))
print(col_names_indices)

meso_abun <- meso_abundance_data %>% mutate_at(vars(10:29), as.numeric)
meso_abun[, 10:29] <- scale(meso_abun[, 10:29])



# Prepare a hurdle model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)

hurdle_meso_full <- hurdle(meso_RelativeAbundance ~ week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + wk_maxt_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = meso_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_meso_full, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)

# Remove variables that are not significant 

hurdle_meso_1 <- hurdle(meso_RelativeAbundance ~ wk_maxt_av, dist = "negbin", zero.dist = "binomial", data = meso_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_meso_1, direction = "both", 
                      trace = FALSE)
summary(step.model)


# Compare to null model

hurdle_meso_null <- hurdle(meso_RelativeAbundance ~   1   , dist = "negbin", zero.dist = "binomial", data = meso_abun, link = "logit", na.action = "na.fail")
AIC(hurdle_meso_1,hurdle_meso_null) 


# Plot the relationship with each variable for occurrence and relative abundance

meso_prob0_wk_maxt_av <- emmip(hurdle_meso_1, ~ wk_maxt_av,
                               at = list(wk_maxt_av = seq(min(meso_abun$wk_maxt_av), max(meso_abun$wk_maxt_av), length.out = 50)
                               ), 
                               lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 


meso_abun_wk_maxt_av <- emmip(hurdle_meso_1, ~ wk_maxt_av,
                                  at = list(wk_maxt_av = seq(min(meso_abun$wk_maxt_av), max(meso_abun$wk_maxt_av), length.out = 50)                                  ), 
                                  lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 


var1_long <- meso_abun_wk_av_temp %>% 
  mutate(variable = "wk_av_temp") %>% 
  dplyr::select(variable, yvar, SE, xvar = wk_av_temp)

meso_var1_long <- var1_long


# Plot relative abundance 
ggplot(meso_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance", title = "Mesophellia") +
  theme_minimal() +
  scale_color_manual(values = c( "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c( "#E7B945")) + theme_minimal() +
  theme(legend.position = "none")  


var2_long <- meso_prob0_wk_av_temp  %>% 
  mutate(variable = "wk_av_temp") %>% 
  dplyr::select(variable, yvar, SE, xvar = wk_av_temp)

meso_combined_zero <- var2_long


# Plot occurrence
ggplot(meso_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Mesophellia Presence", title = "Mesophellia") +
  theme_minimal() +
  scale_color_manual(values = c("#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c( "#E7B945")) + theme_minimal() +
  theme(legend.position = "none")  










## 5.3 Thaxterogaster ----


tax_table <- tax_table(ITSrel.count)
Thaxterogaster_taxa <- taxa_names(tax_table[tax_table[, "Genus"] == "Thaxterogaster", ])
data <- prune_taxa(Thaxterogaster_taxa, ITSrel.count) #

data = prune_taxa(taxa_sums(data) > 0, data) ## remove OTUs with zero count

otu_table_df <- as.data.frame(t(otu_table(data)))
sample_data_df <- as.data.frame(sample_data(data))
tax_table_df <- as.data.frame(tax_table(data))

thax_abundance <- as.data.frame(rowSums(otu_table_df))
colnames(thax_abundance)[1] <- "thax_RelativeAbundance"

thax_abundance_df <- data.frame(SampleID = rownames(thax_abundance),thax_RelativeAbundance = thax_abundance)

thax_abundance_data <- merge(sample_data_df, thax_abundance_df, by.x = "row.names", by.y = "SampleID")
colnames(thax_abundance)[1] <- "SampleID"

col_names_indices <- data.frame(Index = seq_along(colnames(thax_abundance_data)), Column = colnames(thax_abundance_data))
print(col_names_indices)

thax_abun <- thax_abundance_data %>% mutate_at(vars(10:29), as.numeric)
thax_abun[, 10:29] <- scale(thax_abun[, 10:29])


# Prepare a hurdle model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)

hurdle_thax_full <- hurdle(thax_RelativeAbundance ~  week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + wk_maxt_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = thax_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_thax_full, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)

# Remove variables that are not significant, and then stepwise remove vars that contribute the least to the model

hurdle_thax_1 <- hurdle(thax_RelativeAbundance ~   three_mon_mint_av      | three_mon_mint_av+ twelve_mon_mint_av  , dist = "negbin", zero.dist = "binomial", data = thax_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_thax_1, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(hurdle_thax_1)

# Compare to null model

hurdle_thax_null <- hurdle(thax_RelativeAbundance ~   1   , dist = "negbin", zero.dist = "binomial", data = thax_abun, link = "logit", na.action = "na.fail")
AIC(hurdle_thax_1,hurdle_thax_null) 


# Plot the relationship with each variable for occurrence and relative abundance, holding the other variable at its mean


thax_prob0_three_mon_mint_av <- emmip(hurdle_thax_1, ~ three_mon_mint_av,
                                      at = list(three_mon_mint_av = seq(min(thax_abun$three_mon_mint_av), max(thax_abun$three_mon_mint_av), length.out = 50)),        twelve_mon_mint_av = mean(thax_abun$twelve_mon_mint_av),               lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 


thax_prob0_twelve_mon_mint_av <- emmip(hurdle_thax_1, ~ twelve_mon_mint_av,
                                       at = list(twelve_mon_mint_av = seq(min(thax_abun$twelve_mon_mint_av), max(thax_abun$twelve_mon_mint_av), length.out = 50 )), three_mon_mint_av = mean(thax_abun$three_mon_mint_av),
                                       lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 


thax_abun_three_mon_mint_av <- emmip(hurdle_thax_1, ~ three_mon_mint_av,
                                     at = list(three_mon_mint_av = seq(min(thax_abun$three_mon_mint_av), max(thax_abun$three_mon_mint_av), length.out = 50)) , 
                                     lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 


thax_var1_long <- thax_abun_three_mon_mint_av %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)


# Plot relative abundance 

ggplot(thax_var1_long, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted relative abundance", title = "Thaxterogaster") +
  theme_minimal() +
  scale_color_manual(values = c("#799998")) +  # Adjust colors
  scale_fill_manual(values = c("#799998")) + theme_minimal() +
  theme(legend.position = "none")  



# Plot occurrence

var1_long <- thax_prob0_three_mon_mint_av %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)

var2_long <- thax_prob0_twelve_mon_mint_av  %>% 
  mutate(variable = "twelve_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = twelve_mon_mint_av)

thax_combined_zero <- bind_rows(var1_long, var2_long)

ggplot(thax_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Thaxterogaster Presence", title = "Thaxterogaster") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#92c9c4")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#92c9c4")) + theme_minimal() +
  theme(legend.position = "none")



## 5.4 Arcangeliella ----


tax_table <- tax_table(ITSrel.count)
Arcangeliella_taxa <- taxa_names(tax_table[tax_table[, "Genus"] == "Arcangeliella", ])
data <- prune_taxa(Arcangeliella_taxa, ITSrel.count) 
data = prune_taxa(taxa_sums(data) > 0, data) ## remove OTUs with zero count

otu_table_df <- as.data.frame(t(otu_table(data)))
sample_data_df <- as.data.frame(sample_data(data))
tax_table_df <- as.data.frame(tax_table(data))

arc_abundance <- as.data.frame(rowSums(otu_table_df))
colnames(arc_abundance)[1] <- "arc_RelativeAbundance"

arc_abundance_df <- data.frame(SampleID = rownames(arc_abundance),arc_RelativeAbundance = arc_abundance)

arc_abundance_data <- merge(sample_data_df, arc_abundance_df, by.x = "row.names", by.y = "SampleID")
colnames(arc_abundance)[1] <- "SampleID"

col_names_indices <- data.frame(Index = seq_along(colnames(arc_abundance_data)), Column = colnames(arc_abundance_data))
print(col_names_indices)

arc_abun <- arc_abundance_data %>% mutate_at(vars(10:29), as.numeric)
arc_abun[, 10:29] <- scale(arc_abun[, 10:29])


# Prepare a hurdle model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)

hurdle_arc_full <- hurdle(arc_RelativeAbundance ~  week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + wk_maxt_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = arc_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_arc_full, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)

# Remove variables that are not significant 

hurdle_arc_1 <- hurdle(arc_RelativeAbundance ~   three_mon_mint_av  + fdsi_annual | three_mon_mint_av 
                         
                         ,   , dist = "negbin", zero.dist = "binomial", data = arc_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_arc_1, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)


# Compare to null model

hurdle_arc_null <- hurdle(arc_RelativeAbundance ~   1   , dist = "negbin", zero.dist = "binomial", data = arc_abun, link = "logit", na.action = "na.fail")
AIC(hurdle_arc_1,hurdle_arc_null) 


# Plot the relationship with each variable for occurrence and relative abundance, holding the other variable at its mean

# Create a model with just quarterly min temperature in order to plot the zero component
hurdle_arc_plotzero <- hurdle(arc_RelativeAbundance ~       three_mon_mint_av, dist = "negbin", zero.dist = "binomial", data = arc_abun, link = "logit", na.action = "na.fail")

arc_prob0_three_mon_mint_av <- emmip(hurdle_arc_plotzero, ~ three_mon_mint_av,
                                     at = list(three_mon_mint_av = seq(min(arc_abun$three_mon_mint_av), max(arc_abun$three_mon_mint_av), length.out = 50)
                                     ), 
                                     lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 



arc_abun_quarterly_mint <- emmip(hurdle_arc_1, ~ three_mon_mint_av,
                                  at = list(three_mon_mint_av = seq(min(arc_abun$three_mon_mint_av), max(arc_abun$three_mon_mint_av), length.out = 50),
                                            fdsi_annual = mean(arc_abun$fdsi_annual)
                                  ), 
                                  lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 

arc_abun_fdsi_annual <-emmip(hurdle_arc_1, ~ fdsi_annual,
                             at = list(fdsi_annual = seq(min(arc_abun$fdsi_annual), max(arc_abun$fdsi_annual), length.out = 50),
                                       three_mon_mint_av = mean(arc_abun$three_mon_mint_av)
                             ), 
                             lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 


var1_long <- arc_abun_quarterly_mint %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)

var2_long <- arc_abun_fdsi_annual %>% 
  mutate(variable = "fdsi_annual") %>% 
  dplyr::select(variable, yvar, SE, xvar = fdsi_annual)

arc_combined_count <- bind_rows(var1_long, var2_long)


# Plot relative abundance 
ggplot(arc_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance", title = "Arcangeliella") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() #+
  #theme(legend.position = "none")  


var1_long <- arc_prob0_three_mon_mint_av %>% 
  mutate(variable = "three_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = three_mon_mint_av)


arc_combined_zero <- var1_long

# Plot occurrence
ggplot(arc_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Arcangeliella Presence", title = "Arcangeliella") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal()# +
  #theme(legend.position = "none")  




## 5.5 Austrogautieria ----


tax_table <- tax_table(ITSrel.count)
Austrogautieria_taxa <- taxa_names(tax_table[tax_table[, "Genus"] == "Austrogautieria", ])
data <- prune_taxa(Austrogautieria_taxa, ITSrel.count) #

data = prune_taxa(taxa_sums(data) > 0, data) ## remove OTUs with zero count

otu_table_df <- as.data.frame(t(otu_table(data)))
sample_data_df <- as.data.frame(sample_data(data))
tax_table_df <- as.data.frame(tax_table(data))

austro_abundance <- as.data.frame(rowSums(otu_table_df))
colnames(austro_abundance)[1] <- "austro_RelativeAbundance"

austro_abundance_df <- data.frame(SampleID = rownames(austro_abundance),austro_RelativeAbundance = austro_abundance)

austro_abundance_data <- merge(sample_data_df, austro_abundance_df, by.x = "row.names", by.y = "SampleID")
colnames(austro_abundance)[1] <- "SampleID"

col_names_indices <- data.frame(Index = seq_along(colnames(austro_abundance_data)), Column = colnames(austro_abundance_data))
print(col_names_indices)

austro_abun <- austro_abundance_data %>% mutate_at(vars(10:29), as.numeric)
austro_abun[, 10:29] <- scale(austro_abun[, 10:29])


# Prepare a hurdle model with all non correlated climate variables (selected in script 2.alpha_diversity.R , section 3.2)

hurdle_austro_full <- hurdle(austro_RelativeAbundance ~  week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + wk_maxt_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = austro_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_austro_full, direction = "both", 
                      trace = FALSE)
summary(step.model)
AIC(step.model)

# Remove variables that are not significant 

hurdle_austro_1 <- hurdle(austro_RelativeAbundance ~   twelve_mon_precip |  
                            twelve_mon_mint_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = austro_abun, link = "logit", na.action = "na.fail")

step.model <- stepAIC(hurdle_austro_1, direction = "both", 
                      trace = FALSE)


summary(step.model)
AIC(step.model)

# Compare to null model

hurdle_austro_null <- hurdle(austro_RelativeAbundance ~   1   , dist = "negbin", zero.dist = "binomial", data = austro_abun, link = "logit", na.action = "na.fail")
AIC(hurdle_austro_1,hurdle_austro_null) 


# Plot the relationship with each variable for occurrence and relative abundance, holding the other variable at its mean



# Create a model with just zero model variables in order to plot the zero component
hurdle_austro_plotzero <- hurdle(austro_RelativeAbundance ~       twelve_mon_mint_av + fdsi_annual, dist = "negbin", zero.dist = "binomial", data = austro_abun, link = "logit", na.action = "na.fail")


austro_prob0_twelve_mon_mint_av <- emmip(hurdle_austro_plotzero, ~ twelve_mon_mint_av,
                               at = list(twelve_mon_mint_av = seq(min(austro_abun$twelve_mon_mint_av), max(austro_abun$twelve_mon_mint_av), length.out = 50),
                                         fdsi_annual = mean(austro_abun$fdsi_annual)
                               ), 
                               lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 


austro_prob0_fdsi_annual <- emmip(hurdle_austro_plotzero, ~ fdsi_annual,
                                   at = list(fdsi_annual = seq(min(austro_abun$fdsi_annual), max(austro_abun$fdsi_annual), length.out = 50),
                                             twelve_mon_mint_av = mean(austro_abun$twelve_mon_mint_av)
                                   ), 
                                   lin.pred = FALSE, mode ="zero", CIs = TRUE, plotit = FALSE) 

austro_abun_twelve_mon_precip <- emmip(hurdle_austro_1, ~ twelve_mon_precip,
                                  at = list(twelve_mon_precip = seq(min(austro_abun$twelve_mon_precip), max(austro_abun$twelve_mon_precip))
                                  ), 
                                  lin.pred = FALSE, mode = "response", CIs = TRUE, plotit = FALSE) 



var1_long <- austro_abun_twelve_mon_precip %>% 
  mutate(variable = "twelve_mon_precip") %>% 
  dplyr::select(variable, yvar, SE, xvar = twelve_mon_precip)


austro_count <- var1_long


# Plot relative abundance 
ggplot(austro_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance", title = "Austrogautieria") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() #+
#theme(legend.position = "none")  


var1_long <- austro_prob0_fdsi_annual %>% 
  mutate(variable = "fdsi_annual") %>% 
  dplyr::select(variable, yvar, SE, xvar = fdsi_annual)

var2_long <- austro_prob0_twelve_mon_mint_av  %>% 
  mutate(variable = "twelve_mon_mint_av") %>% 
  dplyr::select(variable, yvar, SE, xvar = twelve_mon_mint_av)

austro_combined_zero <- bind_rows(var1_long, var2_long)
str(austro_combined_zero$occurenceprob)

# Plot occurrence
ggplot(austro_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Austrogautieria Presence", title = "Austrogautieria") +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() #+
#theme(legend.position = "none")  

### 5.6.1 Plotting relative abundance component -----



meso_var1_long$Genus <- "Mesophellia"
arc_combined_count$Genus <- " Arcangeliella"
thax_var1_long$Genus <- "Thaxterogaster"
hyst_combined_count$Genus <- "Hysterangium"
austro_count$Genus <- "Austrogautieria"


theme <-   theme(strip.text = element_text(size = 11, colour = "white"),
                 strip.background = element_rect(fill = "#799998", colour = "#799998"),
                 panel.border = element_rect(color = "#799998", fill = NA),
                 plot.background = element_blank()) 



meso <- ggplot(meso_var1_long, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance"
       #, title = "Mesophellia"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#E7B945")) + theme_minimal() +
  theme(legend.position = "none")  + theme  +
  ylim(0, 6500)+   facet_wrap(~Genus, scales = "free_y")



arc <- ggplot(arc_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance"
       #, title = "Arcangeliella"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#CF866F","#799998")) +  # Adjust colors
  scale_fill_manual(values = c("#CF866F","#799998")) + theme_minimal() +
  ylim(0, 1600) +
theme(legend.position = "none") + theme +   facet_wrap(~Genus, scales = "free_y") #brown = fdsi, blue = three month mint av

thax <- ggplot(thax_var1_long, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted relative abundance", 
       #title = "Thaxterogaster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#799998")) +  # Adjust colors
  scale_fill_manual(values = c("#799998")) + theme_minimal() +
  theme(legend.position = "none")  + theme+
  ylim(0, 1600)+   facet_wrap(~Genus, scales = "free_y")



hyst <- ggplot(hyst_combined_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted abundance" 
       # , title = "hysterangium"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() +
  theme(legend.position = "none")  + theme +
  ylim(0, 6500)+   facet_wrap(~Genus, scales = "free_y")


austro <- ggplot(austro_count, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Predicted relative abundance"
       #, title = "Austrogautiera"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#E0D1C3")) +  # Adjust colors
  scale_fill_manual(values = c("#E0D1C3")) + theme_minimal() +
  theme(legend.position = "none")  +
  ylim(0, 1000)+   facet_wrap(~Genus, scales = "free_y")


unique(austro_count$variable)




combined_plot <- (meso +  theme_light()+theme(axis.title.y = element_blank()) |
                            hyst +  theme_light()+theme(axis.title.y = element_blank())|
                            arc +  theme_light()+theme(axis.title.y = element_blank()) |   
                            thax +  theme_light()+theme(axis.title.y = element_blank()) |
                            austro +  theme_light()+theme(axis.title.y = element_blank()) 
)
combined_plot+ 
  plot_layout(guides = "collect") & 
  theme(
    panel.grid.major.x = element_line(color = "#EEE7E2"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EEE7E2"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(size = 12, colour = "white"),
    strip.background = element_rect(fill = "#799998", colour = "#799998"),
    legend.position = "none")







### 5.6.2 Plotting zero component -----



meso_combined_zero$Genus <- "Mesophellia"
arc_combined_zero$Genus <- " Arcangeliella"
thax_combined_zero$Genus <- "Thaxterogaster"
hyst_combined_zero$Genus <- "Hysterangium"
austro_combined_zero$Genus <- "Austrogautieria"

theme <-   theme(strip.text = element_text(size = 11, colour = "white"),
                 strip.background = element_rect(fill = "#799998", colour = "#799998"),
                 panel.border = element_rect(color = "#799998", fill = NA),
                 plot.background = element_blank()) 

meso <-ggplot(meso_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Mesophellia Presence"
       #, title = "Mesophellia") 
  )+
  theme_minimal() +
  scale_color_manual(values = c( "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#E7B945")) + theme_minimal() +
  theme(legend.position = "none")   +  theme + facet_wrap(~Genus, scales = "free_y") + ylim(-1.2,1.5)



arc <-ggplot(arc_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Arcangeliella Absence"
       #, title = "Arcangeliella"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#799998")) +  # Adjust colors
  scale_fill_manual(values = c("#799998")) + theme_minimal() +
  theme(legend.position = "none")+  theme + facet_wrap(~Genus, scales = "free_y")+ ylim(-1.2,1.5)



thax <-ggplot(thax_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Thaxterogaster Presence"
       #, title = "Thaxterogaster"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#92c9c4")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#92c9c4")) + theme_minimal() +
  theme(legend.position = "none") +  theme + facet_wrap(~Genus, scales = "free_y") + ylim(-1.2,1.5)



hyst <-ggplot(hyst_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Hysterangium Presence"
       #, title = "hysterangium"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("#799998", "#E7B945")) +  # Adjust colors
  scale_fill_manual(values = c("#799998", "#E7B945")) + theme_minimal() +
  theme(legend.position = "none") +  theme + facet_wrap(~Genus, scales = "free_y") + ylim(-1.2,1.5)



austro <- 
  ggplot(austro_combined_zero, aes(x = xvar, y = yvar, color = variable)) +
  geom_line(size = 0.7) + 
  geom_ribbon(aes(ymin = yvar - 2 * SE, ymax = yvar + 2 * SE, fill = variable), 
              alpha = 0.3, show.legend = FALSE) +
  labs(x = "Scaled climate variable", y = "Probability of Austrogautiera Presence"
       #, title = "Austrogautiera"
  ) +
  theme_minimal() +
  scale_color_manual(values = c( "#CF866F","#92c9c4")) +  # Adjust colors
  scale_fill_manual(values = c( "#CF866F","#92c9c4"))+
  theme(legend.position = "none") +  theme + facet_wrap(~Genus, scales = "free_y") + ylim(-1.2,1.5)









# to remove labels on the x axis (abundance vals), axis.text.y = element_blank()
combined_plot_chapter <- (meso +  theme_light()+theme(axis.title.y = element_blank()) |
                            hyst +  theme_light()+theme(axis.title.y = element_blank(), axis.text.y = element_blank())|
                            arc +  theme_light()+theme(axis.title.y = element_blank(), axis.text.y = element_blank()) |   
                            thax +  theme_light()+theme(axis.title.y = element_blank(), axis.text.y = element_blank())|
                            austro +  theme_light()+theme(axis.title.y = element_blank(), axis.text.y = element_blank()
                            ) 
)
combined_plot_chapter + 
  plot_layout(guides = "collect") & 
  theme(
    panel.grid.major.x = element_line(color = "#EEE7E2"),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "#EEE7E2"),
    panel.grid.minor.y = element_blank(),
    strip.text = element_text(size = 12, colour = "white"),
    strip.background = element_rect(fill = "#799998", colour = "#799998")
    ,    legend.position = "none"
  )








