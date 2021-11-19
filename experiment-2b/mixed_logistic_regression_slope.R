library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(lme4)
library(car)
library(ggplot2)
library(sjPlot)
library(reghelper)


files <- dir("./results/alphabet/approved", ".*csv$", full.names = TRUE)  # Specify data folder (alphabet or imagenet)
dat <- c()

for (i in 1: length(files)){
  d <- fread(files[i], header = TRUE)
  d <- d[, c("TF", "size", "rate", "pos", "ori")]

  d$size <- as.character(d$size)
  d$rate <- as.character(d$rate)
  d$pos <- as.factor(d$pos)
  d$arrange <- paste(d$size, d$rate, sep = "")
  d$id <- i  # Glmer requires numeric variable for subject id

  d %>%
    mutate(TF = gsub(TF, pattern = "FALSE", replacement = "0", ignore.case = TRUE)) %>%
    mutate(TF = gsub(TF, pattern = "TRUE", replacement = "1", ignore.case = TRUE)) ->
    d
  d$TF <- as.numeric(d$TF)  # Glmer logistic regression requires 0/1 numeric variable

  d %>%
    mutate(arrange = gsub(arrange, pattern = "11", replacement = "[1, 1, 1]", ignore.case = TRUE)) %>%
    mutate(arrange = gsub(arrange, pattern = "12", replacement = "[1, 2, 4]", ignore.case = TRUE)) %>%
    mutate(arrange = gsub(arrange, pattern = "21", replacement = "[2, 2, 2]", ignore.case = TRUE)) %>%
    mutate(arrange = gsub(arrange, pattern = "22", replacement = "[2, 4, 8]", ignore.case = TRUE)) ->
    d

  d %>%
    mutate(pos = gsub(pos, pattern = "1", replacement = "3.83")) %>%
    mutate(pos = gsub(pos, pattern = "0", replacement = "1.41")) %>%
    mutate(pos = gsub(pos, pattern = "2", replacement = "8.06")) ->
    d
  d$pos <- as.numeric(d$pos)

  dat <- rbind(dat, d)
}
dat <- na.omit(dat)
colnames(dat)[1] <- "Accuracy"


# Random slope glmm
f1 <- glmer(Accuracy ~ pos + arrange + pos:arrange + (1 + pos + arrange|id), data = dat, family = "binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
f2 <- glmer(Accuracy ~ poly(pos, 2, raw = TRUE) + arrange + poly(pos, 2, raw = TRUE):arrange + (1 + poly(pos, 2, raw = TRUE) + arrange|id), data = dat, family = "binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  # Polynomial model


summary(f1)  # Wald tests for parameter estimates
Anova(f1)  # Type II Wald chisquare test by car package
plot(f1)  # Residual plot
ranef(f1)  # Distribution of individual estimates (as deviation from the global estimate)

summary(f2)  # Wald tests for parameter estimates
Anova(f2)  # Type II Wald chisquare test by car package
plot(f2)  # Residual plot
ranef(f2)  # Distribution of individual estimates (as deviation from the global estimate)


# Visualize predicted models
plot_model(f1, type="pred", terms=c("pos", "arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732"))
plot_model(f2, type="pred", terms=c("pos[all]", "arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732"))


# Predict slopes
predicted_slopes <- simple_slopes(f1)
# Visualize simple slopes
predicted_slopes[10:13,] %>%
  rename("PredictedSlope"=3, "Std.Error"=4) %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-Std.Error*1.96, ymax=PredictedSlope+Std.Error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18)) +
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)) +
  theme(legend.title=element_text(size=18),legend.text=element_text(size=18)) +
  scale_y_continuous(limits=c(-0.5, 0.5)) ->
  g
g


# Visualize quadratic slopes
# alphabet
PredictedSlope <- c(-0.018145, -0.011213, -0.010163, -0.020050)
row_std_error <- c(0.014135, 0.013714, 0.013744, 0.013791)

arrange <- c("[1,1,1]", "[1,2,4]", "[2,2,2]", "[2,4,8]")
slope_df <- data.frame(PredictedSlope, row_std_error, arrange)

slope_df %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-row_std_error*1.96, ymax=PredictedSlope+row_std_error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18)) +
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)) +
  theme(legend.title=element_text(size=18),legend.text=element_text(size=18)) +
  scale_y_continuous(limits=c(-0.1, 0.1)) ->
  g
g


# Visualize average Acc by arrangement (all positions)
dat_11 <- dat[arrange=="[1, 1, 1]",]
dat_12 <- dat[arrange=="[1, 2, 4]",]
dat_21 <- dat[arrange=="[2, 2, 2]",]
dat_22 <- dat[arrange=="[2, 4, 8]",]

rbind(
  mutate(describe(dat_11[, "Accuracy"])[, c("mean", "se")], arrange="[1, 1, 1]"),
  mutate(describe(dat_12[, "Accuracy"])[, c("mean", "se")], arrange="[1, 2, 4]"),
  mutate(describe(dat_21[, "Accuracy"])[, c("mean", "se")], arrange="[2, 2, 2]"),
  mutate(describe(dat_22[, "Accuracy"])[, c("mean", "se")], arrange="[2, 4, 8]")
) %>%
  rename("averageAcc"="mean") %>%
  ggplot(aes(x=arrange, y=averageAcc)) +
  geom_errorbar(aes(ymin=averageAcc-se*1.96, ymax=averageAcc+se*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18)) +
  theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)) +
  theme(legend.title=element_text(size=18),legend.text=element_text(size=18)) +
  scale_y_continuous(limits=c(0.3, 0.5)) ->
  g_ave
g_ave

