# Purpose: Statistical analysis of differences in truffle-like fungal richness between seasons and in correlation with changes in climate 

# Load packages 


library(phyloseq)
library(ggplot2)
library(vegan)
library(ggpubr)
library(patchwork)
library(pscl)
library(caret)
library(MASS)
library(MuMIn)
library(car)
library(emmeans)

# 1. Filter for truffle-like ECM taxa ----

ecm.list=row.names(guild_table)[which(guild_table$primary_lifestyle=="ectomycorrhizal")] # Filter for ECM taxa
data=prune_taxa(ecm.list, ITSrel.count)

hypo.list=row.names(guild_table)[which(guild_table$Fruitbody_type=="gasteroid-hypogeous")] # Filter for truffle-like genera
data=prune_taxa(hypo.list, data)

data = prune_taxa(taxa_sums(data) > 0, data) # Remove OTUs with zero count
truffle_ecm <- prune_samples(sample_sums(data) > 0, data) # Remove samples with zero OTUs 


otu_table_df <- as.data.frame(otu_table(truffle_ecm))
otu_table_df <- t(otu_table_df)

spec_accum <- specaccum(otu_table_df)
spec_accum_df <- data.frame(Samples = spec_accum$sites,Richness = spec_accum$richness,SD = spec_accum$sd)

species_accumulation <- ggplot(spec_accum_df, aes(x = Samples, y = Richness)) +
  geom_line() +
  geom_ribbon(aes(ymin = Richness - SD, ymax = Richness + SD), alpha = 0.2) +
  labs(title = "Truffle-like ECM OTU Accumulation Curve",
       x = "Number of Samples",
       y = "OTU Richness") +
  theme_light()

species_accumulation

# 2. Calculate richness ----


meta = data.frame(sample_data(truffle_ecm))
filtered_reads <- sample_sums((truffle_ecm))
meta <- cbind(meta, filtered_reads)
diversity <-estimate_richness(truffle_ecm, measures=c("Observed", "Chao1", "Shannon", "Simpson"))
alpha_diversity <- cbind(diversity, meta)

# 2. Richness between seasons ----

# Negative binomial GLM 
nb_model_season <- glm.nb(Observed ~ Season , data = alpha_diversity)
summary(nb_model_season)

# Test differences with emmeans
emm <- emmeans(nb_model_season, ~ Season)
pairwise_results <- pairs(emm)
plot(emm)


# Plot seasonal differences

## Plot differences in richness

alpha_diversity$Season <- factor(alpha_diversity$Season, 
                                levels = c("Summer", "Autumn", "Winter", "Spring"))


richness_plot <- 
  ggplot(alpha_diversity, aes(x = Season, y = Observed)) +
  geom_boxplot(fill  = "#7C845C") +
  theme_minimal() +
  labs(title = "Seasonal Differences in Richness of Truffle-like ECM at Bellbird",
       x = "Season",
       y = "Observed Richness",
       color = "Year") +
  stat_compare_means(comparisons = list(c("Summer", "Spring"), c("Winter", "Autumn"), c("Autumn", "Spring")), method = "t.test", label = "p.signif") + # Based on emmeans results
  theme(
    axis.title.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    text = element_text(size = 15), legend.position = "none") + theme_light()+ coord_cartesian(ylim = c(0, 16))

richness_plot

## Plot differences in relative abundance

truffle_ecm_genus_relabun <- truffle_ecm %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt() %>%
  group_by(Season, Genus) %>%
  summarize(total_abun = sum(Abundance), .groups = "drop") %>%
  group_by(Season) %>%
  mutate(rel_abun = total_abun / sum(total_abun)) %>%
  ungroup()

colour_theme <- c("#C2CFD0", "#92c9c4", "#90D8D2", "#809B97", "#95A58E","#4a8780", "#93715c",  "#006973","#A67048","#c1ad9e","#D1B08F", "#B27440","#778773", "#4C2B16","#CF866F", "#EAA48F", "#e3ad9a","#F8D8D2"
)

genus_relabun$Season <- factor(genus_relabun$Season, levels = c("Summer","Autumn", "Winter", "Spring"))

relabun_plot <- ggplot(truffle_ecm_genus_relabun, aes(x = Season, y = rel_abun, fill = Genus)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = colour_theme) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL,
       y = "Relative Abundance",
       fill = "Genus") +
  guides(fill = guide_legend(ncol = 1))  +
  theme_light()+
  theme(
    legend.title = element_text(size = 10),   
    legend.text = element_text(size = 8),          
    legend.key.size = unit(0.4, "cm"),               
    legend.spacing.y = unit(0.2, "cm"),             
    legend.box.spacing = unit(0.2, "cm")   
  )

## Combine richness and relative abundance plots


richness_plot <- richness_plot + theme(plot.title = element_blank(), axis.title.x = element_blank())
richness_plot <- richness_plot  +
  plot_annotation(tag_levels = 'a')

combined_plot <- richness_plot / relabun_plot +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(tag_levels = 'a')

combined_plot



# 3. GLMs ----

## 3.1 Identify highly correlated variables ----

clim_vars <- alpha_diversity[, 13:33] # Select climate variables

cor_mat <- cor(clim_vars, use = "pairwise.complete.obs") # create the correlation matrix
high_cor <- findCorrelation(cor_mat, cutoff = 0.70, verbose = TRUE, names = TRUE) # Find high correlated pairs (correlation > 0.7)
non_cor_vars <- setdiff(names(clim_vars), high_cor) # List the climate variables that are not highly correlated

## 3.2 Select best model ----

# Scale climate variables

alpha_diversity <- alpha_diversity %>% mutate_at(vars(13:33), as.numeric)
alpha_diversity[, c(13:33)] <- scale(alpha_diversity[, c(13:33)])

# Use dredge to identify non-correlated climate variables that best predict observed richness of truffle-like ECM, using negative binomial regression


model <-glm.nb(Observed ~ week_precip + three_month_precip + three_mon_mint_av + twelve_mon_precip + twelve_mon_mint_av + wk_av_temp + mon_AHMI + fdsi_annual, data = alpha_diversity, na.action = "na.fail")


glm_nb <-dredge(model, beta= "sd", evaluate = TRUE, rank = "AIC",
                   m.lim = c(0,2),fixed = NULL ,    # Maximum of two terms 
                   extra=c("adjR^2", "BIC", "Cp"),
                   trace=TRUE)
View(glm_nb)

# Best model:
nb_model <- glm.nb(Observed ~  three_mon_mint_av + twelve_mon_mint_av, data = alpha_diversity)
summary(nb_model)

# Compare to null model
nb_model_null<- glm.nb(Observed ~  1 , data = alpha_diversity)
summary(nb_model_null)
AIC(nb_model,nb_model_null)
anova(nb_model,nb_model_null) # Selected model is sigificantly better than the null model

# Check for interaction
nb_modelInt <- glm.nb(Observed ~  three_mon_mint_av * twelve_mon_mint_av, data = alpha_diversity)
summary(nb_modelInt)
AIC(nb_model,nb_modelInt)
anova(nb_model,nb_modelInt) # Model with interaction terms provides a significantly better fit

# Check model assumptions
par(mfrow = c(2, 2))
plot(nb_modelInt)

## 3.3 Plots -----

### 3.3.1 Predicted richness for each variable ---- 

# For Quarterly Minimum Temperature: 
newdat_three_mon_mint_av <- data.frame(
  three_mon_mint_av = rep(seq(from = min(alpha_diversity$three_mon_mint_av), to = max(  alpha_diversity$three_mon_mint_av), length.out = 100), 1),
  twelve_mon_mint_av = mean(alpha_diversity$twelve_mon_mint_av)
)

newdat_three_mon_mint_av <- cbind(newdat_three_mon_mint_av, predict(nb_modelInt, newdat_three_mon_mint_av, type = "response"))
colnames(newdat_three_mon_mint_av)[3] <- "Predicted"

melt_newdat_three_mon_mint_av <- melt(newdat_three_mon_mint_av, id.vars = c("three_mon_mint_av", "twelve_mon_mint_av")
                                      , value.name = "Predicted")

# Plot for Quarterly Minimum Temperature
quarterly_temp_plot <- ggplot(melt_newdat_three_mon_mint_av, aes(x = three_mon_mint_av, y = Predicted, colour = variable)) +
  geom_line(color = "#92c9c4", size = 1) +
  labs(title = "Effect of Quarterly Minimum Temperature on Observed Richness", x = "Mean quarterly minimum temperature (°C)", y = "Predicted Richness") + theme_light()


# For Annual Minimum Temperature: 
newdat_twelve_mon_mint_av <- data.frame(
  three_mon_mint_av = mean(  alpha_diversity$three_mon_mint_av),
  twelve_mon_mint_av = rep(seq(from = min(  alpha_diversity$twelve_mon_mint_av), to = max(  alpha_diversity$twelve_mon_mint_av), length.out = 100), 1)
)

newdat_twelve_mon_mint_av <- cbind(newdat_twelve_mon_mint_av, predict(nb_modelInt, newdat_twelve_mon_mint_av, type = "response"))
colnames(newdat_twelve_mon_mint_av)[3] <- "Predicted"
melt_newdat_twelve_mon_mint_av <- melt(newdat_twelve_mon_mint_av, id.vars = c("three_mon_mint_av", "twelve_mon_mint_av"),
                                       variable.name = "Variable", value.name = "Predicted")
# Plot for Annual Minimum Temperature
annual_temp_plot <- ggplot(melt_newdat_twelve_mon_mint_av, aes(x = twelve_mon_mint_av, y = Predicted, colour = Variable)) +
  geom_line(color = "#92c9c4", size = 1) +
  labs(title = "Effect of Annual Minimum Temperature on Observed Diversity", x = "Mean annual minimum temperature (°C)", y = "Predicted Richness") + theme_light()


# For the interaction term: 
three_mon_mint_av_seq <- seq(from = min(  alpha_diversity$three_mon_mint_av), to = max(  alpha_diversity$three_mon_mint_av), length.out = 100)
twelve_mon_mint_av_seq <- seq(from = min(  alpha_diversity$twelve_mon_mint_av), to = max(  alpha_diversity$twelve_mon_mint_av), length.out = 100)

grid <- expand.grid(three_mon_mint_av = three_mon_mint_av_seq, twelve_mon_mint_av = twelve_mon_mint_av_seq)

grid$Predicted <- predict(nb_modelInt, newdata = grid, type = "response")


# Plot the interaction effect
interaction_plot <- ggplot(grid, aes(x = three_mon_mint_av, y = twelve_mon_mint_av, fill = Predicted)) +
  geom_tile() +
  labs(title = "Interaction Effect",
       x = "Mean quarterly minimum temperature (°C)",
       y = "Mean annual minimum temperature (°C)",
       fill = "Predicted Richness") +
  scale_fill_gradient(low = "#EAA48F", high = "#92C9C4") + theme_light()

quarterly_temp_plot <- quarterly_temp_plot + labs(title = NULL) + 
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))

annual_temp_plot    <- annual_temp_plot + labs(title = NULL, y = NULL) + 
  theme(axis.title.x = element_text(size = 10))

interaction_plot    <- interaction_plot + labs(title = NULL) + 
  theme(axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10))


model_plot <- (quarterly_temp_plot+ annual_temp_plot) / interaction_plot+
  plot_annotation(tag_levels = "a") 

model_plot



### 3.3.2 Visualise richness across months and years ---- 


alpha_diversity$newdate <- as.Date(alpha_diversity$Field.collection.date,"%d/%m/%Y")
alpha_diversity$julian <- as.numeric(format(alpha_diversity$newdate, "%j"))

alpha_diversity_summary <- alpha_diversity %>%
  mutate(month = format(as.Date(julian, origin = "1970-01-01"), "%m")) %>%
  group_by(month) %>%
  summarize(three_mon_mint_av = mean(three_mon_mint_av, na.rm = TRUE))

alpha_diversity_summary <- alpha_diversity_summary %>%
  mutate(julian_day = as.numeric(month) * 30 - 15)

alpha_diversity$Genus_months <- "Across the year"
alpha_diversity$Genus_years <- "Between years"

# Set breaks to only label every second month, while plotting richness by d of scat collection

days_in_months <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
month_middle_points <- cumsum(c(0, days_in_months[-length(days_in_months)])) + days_in_months / 2
selected_breaks <- month_middle_points[seq(1, length(month_middle_points), by = 2)]
selected_labels <- month.abb[seq(1, 12, by = 2)]

monthly_richness <- ggplot(alpha_diversity, aes(x = as.numeric(julian))) +
  geom_point(aes(y = three_mon_mint_av), color = "#006973", alpha = 0.2, size =2.5 ) +
  geom_smooth(aes(y = three_mon_mint_av), method = "loess", se = TRUE, color = "#006973", fill = "#006973", alpha = 0.2,size =1) +
  geom_point(aes(y = Observed), color = "#C4694C", alpha = 0.2,size =2.5) +
  geom_smooth(aes(y = Observed), method = "loess", se = TRUE, color = "#C4694C", fill = "#C4694C", alpha = 0.2,size =1) +
  scale_y_continuous(
    name = "Truffle-like ECM Richness",
    sec.axis = sec_axis(~., name = expression(paste("Quarterly Minimum Temperature (°C)")))
  ) +
  scale_x_continuous(
    name = NULL,
    breaks = selected_breaks,
    labels = selected_labels
  ) +
  theme(
    axis.title.y = element_text(color = "#C4694C", size = 14),
    axis.title.y.right = element_text(color = "#006973", face="bold", size = 14))+
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +  theme_light()

# Update scale_x_continuous with breaks and labels
scale_x_continuous(
  name = NULL,
  breaks = selected_breaks,
  labels = selected_labels
)

sec_axis_trans <- function(x) {
  (x - 9) / (10.7 - 9) * (12 - 5) + 5
}

# Create the inverse scaling function for the secondary axis labels
inv_sec_axis_trans <- function(x) {
  (x - 5) / (12 - 5) * (10.7 - 9) + 9
}

# Create the plot with adjusted secondary y-axis limits
annual_richness <- ggplot(alpha_diversity, aes(x = as.numeric(Year))) +
  geom_point(aes(y = Observed), color = "#C4694C", alpha = 0.2,size =2.5)+
  geom_smooth(aes(y = Observed), method = "loess", se = TRUE, color = "#C4694C", fill = "#C4694C", alpha = 0.2,size =1)  +
  geom_point(aes(y = sec_axis_trans(twelve_mon_mint_av)), color = "#006973", alpha = 0.2, size =2.5 )+
  geom_smooth(aes(y = sec_axis_trans(twelve_mon_mint_av)), method = "loess", se = TRUE, color = "#006973", fill = "#006973", alpha = 0.2,size =1) +
  scale_x_continuous(breaks = seq(1993, 2016, by = 3),  # Label every 5th year
                     labels = seq(1993, 2016, by = 3)) +
  scale_y_continuous(
    name = "Truffle-like ECM Richness",
    sec.axis = sec_axis(~inv_sec_axis_trans(.), name = expression(paste("Annual Minimum Temperature (°C)")))
  ) +
  theme(
    axis.title.y = element_text(color = "#C4694C"),
    axis.title.y.right = element_text(color = "#006973")
  ) +
  labs(x = NULL)+theme(
    axis.title.y = element_text(color = "#C4694C", size = 14),
    axis.title.y.right = element_text(color = "#006973", face="bold", size = 14),
    text = element_text(colour = "#C4694C", size = 14))+
  theme(plot.title = element_text(hjust = 0.5, size = 16)) +  theme_light()
summary(alpha_diversity$three_mon_mint_av)


theme <- theme(
  axis.title.y = element_text(color = "#C4694C", size = 12),
  axis.title.y.right = element_text(color = "#006973", face="bold", size = 12),
  plot.title = element_text(hjust = 0.5, size = 12)
)

monthly_richness<- monthly_richness + theme
annual_richness<- annual_richness + theme

combined_plot <- (monthly_richness + annual_richness) + plot_annotation(tag_levels = "a")
combined_plot


