library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(lme4)
library(car)
library(ggplot2)
library(sjPlot)
library(reghelper)

files <- dir("../../results/imagenet/approved", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)
files <- dir("../../results/alphabet/approved", ".*csv$", full.names=TRUE)
dat <- c()

for (i in 1: length(files)){
  
  d <- fread(files[i], header=TRUE)
  d <- d[, c("TF", "size", "rate", "pos", "ori")]
  
  d$size <- as.character(d$size)
  d$rate <- as.character(d$rate)
  d$pos <- as.factor(d$pos)
  d$arrange <- paste(d$size, d$rate, sep="")
  d$id <- i  # Glmer requires numeric variable for subject id
  
  d %>%
    mutate(TF=gsub(TF, pattern="FALSE", replacement="0", ignore.case=TRUE)) %>%
    mutate(TF=gsub(TF, pattern="TRUE", replacement="1", ignore.case=TRUE)) ->
    d
  d$TF <- as.numeric(d$TF)  # Glmer logistic regression requires 0/1 numeric variable
  
  d %>%
    mutate(arrange=gsub(arrange, pattern="11", replacement="[1, 1, 1]", ignore.case=TRUE)) %>%
    mutate(arrange=gsub(arrange, pattern="12", replacement="[1, 2, 4]", ignore.case=TRUE)) %>%
    mutate(arrange=gsub(arrange, pattern="21", replacement="[2, 2, 2]", ignore.case=TRUE)) %>%
    mutate(arrange=gsub(arrange, pattern="22", replacement="[2, 4, 8]", ignore.case=TRUE)) ->
    d
 
  d %>%
    mutate(pos=gsub(pos, pattern="1", replacement="3.83")) %>%
    mutate(pos=gsub(pos, pattern="0", replacement="1.41")) %>%
    mutate(pos=gsub(pos, pattern="2", replacement="8.06")) ->
    d
  d$pos <- as.numeric(d$pos)
  
  dat <- rbind(dat, d)
  
}

dat <- na.omit(dat)
colnames(dat)[1] <- "Accuracy"

filter(count(dat, Accuracy, id), Accuracy==TRUE) %>%
  mutate(Acc=n/108) ->
  average_acc
average_acc <- mean(average_acc$Acc)

# Random slope glmm
f0 <- glmer(Accuracy ~ poly(pos, 2, raw=TRUE) + arrange + (1 + poly(pos, 2, raw=TRUE) + arrange|id), data=dat, family="binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun= e5)))
f1 <- glmer(Accuracy ~ pos + arrange + pos:arrange + (1 + pos + arrange|id), data=dat, family="binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
f2 <- glmer(Accuracy ~ I(pos^2) + pos + arrange + I(pos^2):arrange + pos:arrange + (1 + I(pos^2) + pos + arrange|id), data=dat, family="binomial")  # Polynomial model
f2 <- glmer(Accuracy ~ poly(pos, 2, raw=TRUE) + arrange + poly(pos, 2, raw=TRUE):arrange + (1 + poly(pos, 2, raw=TRUE) + arrange|id), data=dat, family="binomial", control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))  # Polynomial model

summary(f0)

summary(f1)  # Wald tests for parameter estimates
Anova(f1)  # Type II Wald chisquare test by car package
plot(f1)  # Residual plot
ranef(f1)  # Distribution of individual estimates (as deviation from the global estimate)

summary(f2)  # Wald tests for parameter estimates
Anova(f2)  # Type II Wald chisquare test by car package
plot(f2)  # Residual plot
ranef(f2)  # Distribution of individual estimates (as deviation from the global estimate)

# Visualize predicted models
plot_model(f0, type="pred", terms=c("pos[all]", "arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732"))
plot_model(f1, type="pred", terms=c("pos", "arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  scale_x_continuous(breaks=c(1.41,3.83,8.06), limits=c(1.2,8.2)) +
  theme(axis.title.x=element_text(size=36), axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentricity (deg)", y="Accuracy")
plot_model(f2, type="pred", terms=c("pos[all]", "arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732"))

# Predict slopes
predicted_slopes <- simple_slopes(f1)
# Visualize simple slopes
# predicted_slopes[10:13,] %>%
#   rename("PredictedSlope"=3, "Std.Error"=4) %>%
#   ggplot(aes(x=arrange, y=PredictedSlope)) +
#   geom_errorbar(aes(ymin=PredictedSlope-Std.Error*1.96, ymax=PredictedSlope+Std.Error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
#   geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
#   theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18)) +
#   theme(axis.text.x=element_text(size=18),axis.text.y=element_text(size=18)) +
#   theme(legend.title=element_text(size=18),legend.text=element_text(size=18)) +
#   scale_y_continuous(limits=c(-0.5, 0.5)) ->
#   g
# g
arrange <- c("[1,1,1]", "[1,2,4]", "[2,2,2]", "[2,4,8]")
predicted_slopes$arrange <- as.numeric(predicted_slopes$arrange)
predicted_slopes[10:13,] %>%
  rename("PredictedSlope"=2, "Std.Error"=3) %>%
  cbind(., arrange) %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-Std.Error*1.96, ymax=PredictedSlope+Std.Error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  scale_y_continuous(limits=c(-0.5, 0.5)) +
  labs(x="Arrange", y="First-Order Coefficient") ->
  g
g

# Visualize quadratic slopes
# imagenet
PredictedSlope <- c(0.013783, -0.02814, -0.0005626, 0.014130)
row_std_error <- c(0.018379, 0.01668, 0.0167775, 0.016100)
# alphabet
PredictedSlope <- c(-0.018145, -0.011213, -0.010163, -0.020050)
row_std_error <- c(0.014135, 0.013714, 0.013744, 0.013791)

arrange <- c("[1,1,1]", "[1,2,4]", "[2,2,2]", "[2,4,8]")
slope_df <- data.frame(PredictedSlope, row_std_error, arrange)

slope_df %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-row_std_error*1.96, ymax=PredictedSlope+row_std_error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  scale_y_continuous(limits=c(-0.1, 0.1)) +
  labs(x="Arrange", y="Second-Order Coefficient") ->
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
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  # scale_y_continuous(limits=c(0, 0.4)) +  # Imagenet
  scale_y_continuous(breaks=seq(0, 0.6, by=0.1), limits=c(0, 0.6)) +  # Alphabet
  labs(x="Arrange", y="Accuracy") ->
  g_ave
g_ave

# Simple trend
gs111 <- fixef(f1)[2]
gs124 <- fixef(f1)[2] + fixef(f1)[6]
gs222 <- fixef(f1)[2] + fixef(f1)[7]
gs248 <- fixef(f1)[2] + fixef(f1)[8]

ls111 <- ranef(f1)$id[2] + gs111
ls124 <- ranef(f1)$id[2] + gs124
ls222 <- ranef(f1)$id[2] + gs222
ls248 <- ranef(f1)$id[2] + gs248

ls111 <- as.data.frame(ls111)
ls124 <- as.data.frame(ls124)
ls222 <- as.data.frame(ls222)
ls248 <- as.data.frame(ls248)

colnames(ls111) <- "slope"
colnames(ls124) <- "slope"
colnames(ls222) <- "slope"
colnames(ls248) <- "slope"

ls111$arrange <- "[1, 1, 1]"
ls124$arrange <- "[1, 2, 4]"
ls222$arrange <- "[2, 2, 2]"
ls248$arrange <- "[2, 4, 8]"

slopes <- rbind(ls111, ls124, ls222, ls248)

g <- ggplot(slopes, aes(x=arrange, y=slope)) + geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="#EC4C4D", width=0.2)
g