library(tidyverse)
library(patchwork)
library(magrittr)
library(car)
library(effectsize)
library(margins)
library(emmeans)
library(ggeffects)
library(interactions)
####################################
# Load Data
####################################
dat <- read_csv("poorbeliefs.csv")

####################################
# Fit Model
####################################
model <- lm(Zagency ~ ZWpoor * Zzerosum + Zedu + Zincome + Democrat + ZBpoor, data = dat)

####################################
# Create Data for PDF of Residuals Plot
####################################
ggdat <- data.frame(e = rstudent(model))
ggdat.gaussian <- data.frame(x = seq(min(ggdat$e), max(ggdat$e), length.out = 1000), f = dnorm(seq(min(ggdat$e), max(ggdat$e), length.out = 1000), mean = 0, sd = 1))
####################################
# Create PDF of Residuals Plot
####################################
default.bins <- round(log2(nrow(ggdat) + 1))
p1 <- ggplot(ggdat, aes(x = e)) +
  geom_histogram(aes(y = after_stat(density)), bins = default.bins, fill = "grey30", color = "lightgray") +
  geom_density(aes(color = "Empirical"), linewidth = 1, show.legend = FALSE) +
  stat_density(aes(x = e, color = "Empirical"), geom = "line", position = "identity") +
  geom_line(data = ggdat.gaussian, aes(x = x, y = f, color = "Gaussian-Assumed"), linetype = "dashed", size = 1) +
  theme_bw() +
  xlab("Studentized Residual") +
  ylab("Density") +
  labs(color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual("", breaks = c("Empirical", "Gaussian-Assumed"), values = c("black", "lightslategray"))
####################################
# Create Data for CDF of Residuals Plot
####################################
e.cdf.func <- ecdf(rstudent(model))
e.cdf <- e.cdf.func(sort(rstudent(model)))
ggdat <- data.frame(e = sort(rstudent(model)), e.cdf = e.cdf)
ggdat.gaussian <- data.frame(x = seq(min(ggdat$e), max(ggdat$e), length.out = 1000), CDF = pnorm(seq(min(ggdat$e), max(ggdat$e), length.out = 1000), mean = 0, sd = 1))
####################################
# Create CDF of Residuals Plot
####################################
p2 <- ggplot(data = ggdat, aes(x = e)) +
  geom_step(aes(y = e.cdf, color = "Empirical")) +
  geom_line(data = ggdat.gaussian, aes(x = x, y = CDF, color = "Gaussian-Assumed"), linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0) +
  theme_bw() +
  xlab("Studentized Residual") +
  ylab("Cumulative Density") +
  labs(color = "") +
  theme(legend.position = "bottom") +
  scale_color_manual("", breaks = c("Empirical", "Gaussian-Assumed"), values = c("black", "lightslategray"))
####################################
# Create QQ Plot
####################################
p3 <- ggplot(data = ggdat, aes(sample = e)) +
  geom_qq() +
  geom_qq_line() +
  theme_bw() +
  xlab("Gaussian Quantiles") +
  ylab("Sample Quantiles")
####################################
# Create Data for Fitted vs Residual Plot
####################################
ggdat <- data.frame(x = fitted(model), e = rstudent(model))
ggdat.out3 <- ggdat %>%
  filter(abs(e) > 3)
####################################
# Create Fitted vs Residual Plot
####################################
p4 <- ggplot(data = ggdat, aes(x = x, y = e)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab(bquote("Fitted Values" ~ (hat(Y)))) +
  ylab("Studentized Residual") +
  theme_bw() +
  geom_hline(yintercept = c(-3, 3), color = "black", linetype = "dotted", size = 0.75) +
  geom_point(data = ggdat.out3, aes(x = x, y = e), fill = "black", shape = 8)
####################################
# Print Plots
####################################
pdf("cooleyresid.pdf", width = 6, height = 4)
(p1 + p2 + p3 + p4) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
dev.off()


####################################
# Create Data for Leverage Plot
####################################
d <- model$model
n <- nrow(d)
p <- length(coef(model))
d <- d %>%
  mutate(obs = 1:n)
ggdat <- d %>%
  mutate(h.values = hatvalues(model))
####################################
# Create Leverage Plot
####################################
p1 <- ggplot(data = ggdat, aes(x = obs)) +
  geom_linerange(aes(ymin = 0, ymax = h.values)) +
  theme_bw() +
  xlab("Observation Number") +
  ylab("Leverage") +
  geom_hline(yintercept = 2 * p / n, linetype = "dotted", color = "gray65", size = 0.75) +
  geom_hline(yintercept = 3 * p / n, linetype = "dotted", color = "black", size = 0.75)
####################################
# Create Data for Cook's D Plot
####################################
ggdat <- d %>%
  mutate(cook.d = cooks.distance(model))
####################################
# Create Cook's D Plot
####################################
p2 <- ggplot(data = ggdat, aes(x = obs)) +
  geom_linerange(aes(ymin = 0, ymax = cook.d)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)) +
  xlab("Observation Number") +
  ylab("Cook's Distance") +
  geom_hline(yintercept = qf(p = 0.1, df1 = p, df2 = n - p), linetype = "dotted", color = "gray65", size = 0.75) +
  geom_hline(yintercept = qf(p = 0.5, df1 = p, df2 = n - p), linetype = "dotted", color = "black", size = 0.75)
####################################
# Create Data for DFFITS Plot
####################################
ggdat <- d %>%
  mutate(dffits = dffits(model))
####################################
# Create DFFITS Plot
####################################
p3 <- ggplot(data = ggdat, aes(x = obs)) +
  geom_linerange(aes(ymin = 0, ymax = dffits)) +
  theme_bw() +
  xlab("Observation Number") +
  ylab("DFFITs") +
  geom_hline(yintercept = c(-2 * sqrt(p / n), 2 * sqrt(p / n)), linetype = "dotted", size = 0.75, color = "gray65") +
  geom_hline(yintercept = c(-2, 2), linetype = "dotted", size = 0.75, color = "black")
####################################
# Create Data for Residual Plot
####################################
ggdat <- data.frame(obs = d$obs, y = rstudent(model))
ggdat.out2 <- ggdat %>%
  filter(abs(y) > 2)
ggdat.out3 <- ggdat %>%
  filter(abs(y) > 3)
####################################
# Create Residual Plot
####################################
p4 <- ggplot(data = ggdat, aes(x = obs, y = y)) +
  geom_point(shape = 1) +
  geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
  xlab("Observation Number") +
  ylab("Studentized Residual") +
  theme_bw() +
  geom_hline(yintercept = c(-3, 3), color = "black", linetype = "dotted", size = 0.75) +
  geom_hline(yintercept = c(-2, 2), color = "gray65", linetype = "dotted", size = 0.75) +
  geom_point(data = ggdat.out3, aes(x = obs, y = y), shape = 8)
####################################
# Print Plots
####################################
pdf("cooleyoutlier.pdf", width = 6, height = 4)
(p1 | p2) / (p3 | p4)
dev.off()


####################################
# ANOVA Table
####################################
anova.table <- data.frame(Anova(model, type = "III"))

####################################
# Clean Up Labels for Printing
####################################
anova.table <- anova.table %>%
  mutate(Term = rownames(.)) %>%
  relocate(Term) %>%
  mutate(EffectSize = c(NA, epsilon_squared(model)$Epsilon2_partial, NA)) %>%
  set_rownames(NULL) %>%
  set_colnames(c("Term", "SS (Type III)", "df", "F", "p-value", "Partial Epsilon-Squared"))

####################################
# Print Table
####################################
anova.table %>%
  mutate_if(is.numeric, round, 4) %>%
  mutate(`p-value` = ifelse(`p-value` < 0.0001, "<0.0001", `p-value`))

####################################
# Calculate Marginal Effects
####################################
m <- mean(unlist(model$model["Zzerosum"]), na.rm = T)
s <- sd(unlist(model$model["Zzerosum"]), na.rm = T)
modvarat <- list(c(round(m - s, 2), round(m + s, 2)))
names(modvarat) <- c("Zzerosum")
mod.emtrends <- data.frame(emtrends(object = model, spec = "Zzerosum", var = "ZWpoor", at = modvarat))
mod.emtrends[, "Zzerosum"] <- ifelse(mod.emtrends[, "Zzerosum"] == round(m - s, 2), "Low (Mean - 1SD)", "High (Mean + 1SD)")
####################################
# Clean Up Labels for Printing
####################################
colnames(mod.emtrends)[-(1)] <- c(paste("Slope of", "ZWpoor"), "SE", "df", "Lower CI", "Upper CI")
emtrend.test <- test(emtrends(object = model, spec = "Zzerosum", var = "ZWpoor", at = modvarat))
mod.emtrends <- mod.emtrends %>%
  mutate(`p-value` = emtrend.test$p.value, `t ratio` = emtrend.test$t.ratio) %>%
  relocate(c(`t ratio`, `p-value`), .after = df)
####################################
# Print Table
####################################
mod.emtrends %>%
  mutate_if(is.numeric, round, 4) %>%
  mutate(`p-value` = ifelse(`p-value` < 0.0001, "<0.0001", `p-value`))

####################################
# Calculate Marginal Effects Contrasts
####################################
m.mod <- mean(unlist(model$model["Zzerosum"]), na.rm = T)
s.mod <- sd(unlist(model$model["Zzerosum"]), na.rm = T)
m.var <- mean(unlist(model$model["ZWpoor"]), na.rm = T)
s.var <- sd(unlist(model$model["ZWpoor"]), na.rm = T)
modvarat <- list(c(round(m.mod - s.mod, 2), round(m.mod + s.mod, 2)), c(round(m.var - s.var, 2), round(m.var + s.var, 2)))
names(modvarat) <- c(("Zzerosum"), "ZWpoor")
mod.emtrends <- emtrends(object = model, spec = "Zzerosum", var = "ZWpoor", at = modvarat)
mod.emtrends <- add_grouping(mod.emtrends, "moderator", "Zzerosum", c(paste("(Low ", "Zzerosum", ")", sep = ""), paste("(High ", "Zzerosum", ")", sep = "")))
mod.emtrends <- emmeans(mod.emtrends, spec = "moderator", var = "ZWpoor", at = modvarat)
mod.emtrendcontrast <- data.frame(pairs(mod.emtrends))
mod.emtrendcontrastci <- data.frame(confint(pairs(mod.emtrends)))
mod.emtrendcontrast$lower.CL <- mod.emtrendcontrastci$lower.CL
mod.emtrendcontrast$upper.CL <- mod.emtrendcontrastci$upper.CL
####################################
# Clean Up Labels for Printing
####################################
mod.emtrendcontrast <- mod.emtrendcontrast %>%
  set_rownames(NULL) %>%
  set_colnames(c("Contrast", "Estimate", "SE", "df", "t ratio", "p-value", "Lower CI", "Upper CI"))

####################################
# Print Table
####################################
mod.emtrendcontrast %>%
  mutate_if(is.numeric, round, 4) %>%
  mutate(`p-value` = ifelse(`p-value` < 0.0001, "<0.0001", `p-value`))

####################################
# Effects of Interest
####################################
m.mod <- mean(unlist(model$model["Zzerosum"]), na.rm = T)
s.mod <- sd(unlist(model$model["Zzerosum"]), na.rm = T)
meffectsfor <- c("ZWpoor", paste("Zzerosum", "[", round(m.mod - s.mod, 2), ",", round(m.mod + s.mod, 2), "]", sep = ""))
mod.emmeans <- ggemmeans(model = model, terms = meffectsfor, interval = "confidence")

####################################
# Plot
####################################
pdf("cooleymeffs.pdf", width = 6, height = 3)
ggdat <- data.frame(mod.emmeans) %>%
  mutate(group = ifelse(group == -1, "Low (Mean - 1SD)", "High (Mean + 1SD)"))
ggplot(data = ggdat, aes(x = x, y = predicted, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.4) +
  xlab("White-Poor Beliefs (Z)") +
  ylab("Predicted Agency (Z)") +
  scale_color_grey("Zero-Sum Beliefs (Z)") +
  scale_fill_grey("Zero-Sum Beliefs (Z)") +
  theme_bw()
dev.off()

####################################
# Johnson Neyman
####################################
jn <- johnson_neyman(model = model, pred = !!sym("ZWpoor"), modx = !!sym("Zzerosum"), control.fdr = TRUE, plot = FALSE)
####################################
# Plot
####################################
p <- jn$plot + xlab("Zero-Sum Beliefs (Z)") + ylab(paste("Slope of ", "White-Poor Beliefs (Z)", sep = "")) + ggtitle("Johnson-Neyman Plot") + scale_color_grey("Slope of White-Poor Beliefs (Z)") + scale_fill_grey("Slope of White-Poor Beliefs (Z)") + theme_bw() +
  geom_vline(xintercept = 0.1246184, color="black")

pdf("cooleyjn.pdf", width = 6, height = 3)
p
dev.off()
