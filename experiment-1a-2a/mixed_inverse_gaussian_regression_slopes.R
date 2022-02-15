library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(lmerTest)
library(lme4)
library(car)
library(ggplot2)
library(sjPlot)
library(ggeffects)
library(reghelper)

files <- dir("../../results/imagenet/approved", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)
files <- dir("../../results/alphabet/approved", ".*csv$", full.names=TRUE)
dat <- c()

for (i in 1: length(files)){
  
  d <- fread(files[i], header=TRUE)
  d <- d[, c("TF", "size", "rate", "pos", "ori", "reactionTime")]
  
  d$size <- as.character(d$size)
  d$rate <- as.character(d$rate)
  d$pos <- as.factor(d$pos)
  d$arrange <- paste(d$size, d$rate, sep = "")
  d$id <- i  # Glmer requires numeric variable for subject id

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

# When using imagenet results, execute below one line
dat <- subset(dat, id != 5 & id != 16 & id != 19 & id != 45)  # Extract participants who's median of the reactionTime < 10

dat <- na.omit(dat)  # Fill NA
filter(count(dat, TF, id), TF==TRUE) %>% 
  mutate(Acc=n/108) ->
  average_acc
average_acc <- mean(average_acc$Acc)
average_rt <- mean(tapply(dat$reactionTime, dat$id, mean))

dat <- dat[dat$TF=="TRUE",]  # Extract TF==TRUE
dat <- mutate(dat, pos=pos/10)
# dat <- mutate(dat, reactionTime=reactionTime/10)


# Random slope glmm
# No interaction
f0 <- glmer(reactionTime ~ poly(pos, 2, raw = TRUE) + arrange + (1 + poly(pos, 2, raw = TRUE) + arrange|id), data = dat, family = inverse.gaussian(link="identity"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  # Polynomial model
f1 <- glmer(reactionTime ~ pos + arrange + pos:arrange + (1 + pos + arrange|id), data=dat, family=inverse.gaussian(link="identity"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  # Imagenet
# f1 <- glmer(reactionTime ~ pos + arrange + pos:arrange + (1 + pos + arrange|id), data=dat, family=inverse.gaussian(link="identity"))  # Alphabet
f2 <- glmer(reactionTime ~ poly(pos, 2, raw = TRUE) + arrange + poly(pos, 2, raw = TRUE):arrange + (1 + poly(pos, 2, raw = TRUE) + arrange|id), data = dat, family = inverse.gaussian(link="identity"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))  # Polynomial model

summary(f1)  # Wald tests for parameter estimates
Anova(f1)  # Type II Wald chisquare test by car package
plot(f1)  # Residual plot
ranef(f1)  # Distribution of individual estimates (as deviation from the global estimate)

summary(f2)  # Wald tests for parameter estimates
Anova(f2)  # Type II Wald chisquare test by car package
plot(f2)  # Residual plot
ranef(f2)  # Distribution of individual estimates (as deviation from the global estimate)

# Visualize predicted models
plot_model(f0, type="pred", terms=c("pos [all]","arrange"), colors=c("#E93732", "#5EC93B", "#F07A3B", "#1244F4"))
plot_model(f1, type="pred", terms=c("pos","arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732"))
# Imagenet
plot_model(f2, type="pred", terms=c("pos [all]","arrange"), colors=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  scale_x_continuous(breaks=c(0.141,0.383,0.806), limits=c(0.12,0.82), labels=c(1.41,3.83,8.06)) +
  theme(axis.title.x=element_text(size=36), axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
  theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  labs(x="Eccentricity (deg)", y="Reaction Time (sec)")
# Alphabet
plot_model(f2, type="pred", terms=c("pos [all]","arrange"), colors=c("#E93732", "#5EC93B", "#F07A3B", "#1244F4")) +
  scale_x_continuous(breaks=c(0.141,0.383,0.806), limits=c(0.12,0.82), labels=c(1.41,3.83,8.06)) +
  theme(axis.title.x=element_text(size=36), axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
  theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  labs(x="Eccentricity (deg)", y="Reaction Time (sec)") +
  guides(colour=guide_legend(reverse=TRUE))

# Predict slopes
predicted_slopes <- simple_slopes(f1)
# Visualize simple slopes
# predicted_slopes[10:13,] %>%
#   rename("PredictedSlope"=3, "Std.Error"=4) %>%
#   mutate(PredictedSlope=PredictedSlope/10) %>%
#   mutate(Std.Error=Std.Error/10) %>%
#   ggplot(aes(x=arrange, y=PredictedSlope)) +
#   geom_errorbar(aes(ymin=PredictedSlope-Std.Error*1.96, ymax=PredictedSlope+Std.Error*1.96), size=1.2, width = 0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
#   geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
#   theme(axis.title.x = element_text(size=18),axis.title.y = element_text(size=18)) +
#   theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18)) +
#   theme(legend.title = element_text(size=18),legend.text = element_text(size=18)) +
#   scale_y_continuous(limits=c(-0.20, 0.20)) ->
#   g
# g
arrange <- c("[1,1,1]", "[1,2,4]", "[2,2,2]", "[2,4,8]")
predicted_slopes$arrange <- as.numeric(predicted_slopes$arrange)
predicted_slopes[10:13,] %>%
  rename("PredictedSlope"=2, "Std.Error"=3) %>%
  mutate(PredictedSlope=PredictedSlope/10) %>%
  mutate(Std.Error=Std.Error/10) %>%
  cbind(., arrange) %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-Std.Error*1.96, ymax=PredictedSlope+Std.Error*1.96), size=1.2, width=0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  # scale_y_continuous(limits=c(-0.2, 0.2)) +. # Imagenet
  scale_y_continuous(limits=c(-0.1, 0.1)) +  # Alphabet
  labs(x="Arrange", y="First-Order Coefficient") ->
  g
g

# Visualize quadratic slopes
# summary() -> Fixed effects -> poly(pos,2,raw=TRUE)2
# imagenet
PredictedSlope <- c(-2.0922/100, 7.3798/100, 1.5678/100, 2.9177/100)
row_std_error <- c(1.3593/100, 1.1090/100, 1.0284/100, 1.0231/100)
# alphabet
PredictedSlope <- c(0.50005/100, 1.6637/100, 1.1211/100, 2.8002/100)
row_std_error <- c(0.51506/100, 0.5333/100, 0.5098/100, 0.4667/100)

arrange <- c("[1,1,1]", "[1,2,4]", "[2,2,2]", "[2,4,8]")
slope_df <- data.frame(PredictedSlope, row_std_error, arrange)
slope_df %>%
  ggplot(aes(x=arrange, y=PredictedSlope)) +
  geom_errorbar(aes(ymin=PredictedSlope-row_std_error*1.96, ymax=PredictedSlope+row_std_error*1.96), size=1.2, width = 0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  # scale_y_continuous(limits=c(-0.10, 0.10)) +  # Imagenet
  scale_y_continuous(limits=c(-0.05, 0.05)) +  # Alphabet
  labs(x="Arrange", y="Second-Order Coefficient") ->
  g
g

# Visualize reactionTime by arrangement
dat[arrange=="[1, 1, 1]",] %>%
  ggplot(aes(x=pos, y=reactionTime)) +
  geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="#EC4C4D", width=0.2) ->
  g_11
g_11
dat[arrange=="[1, 2, 4]",] %>%
  ggplot(aes(x=pos, y=reactionTime)) +
  geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="#EC4C4D", width=0.2) ->
  g_12
g_12
dat[arrange=="[2, 2, 2]",] %>%
  ggplot(aes(x=pos, y=reactionTime)) +
  geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="#EC4C4D", width=0.2) ->
  g_21
g_21
dat[arrange=="[2, 4, 8]",] %>%
  ggplot(aes(x=pos, y=reactionTime)) +
  geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="#EC4C4D", width=0.2) ->
  g_22
g_22

# Visualize average reactionTime by arrangement
dat_11 <- dat[arrange=="[1, 1, 1]",]
dat_12 <- dat[arrange=="[1, 2, 4]",]
dat_21 <- dat[arrange=="[2, 2, 2]",]
dat_22 <- dat[arrange=="[2, 4, 8]",]
rbind(
  mutate(aggregate(reactionTime~pos, data=dat_11, FUN=mean), arrange="[1, 1, 1]"),
  mutate(aggregate(reactionTime~pos, data=dat_12, FUN=mean), arrange="[1, 2, 4]"),
  mutate(aggregate(reactionTime~pos, data=dat_21, FUN=mean), arrange="[2, 2, 2]"),
  mutate(aggregate(reactionTime~pos, data=dat_22, FUN=mean), arrange="[2, 4, 8]")
) %>%
  ggplot(aes(x=pos, y=reactionTime, group=arrange, colour=arrange)) +
  geom_line() + 
  scale_color_manual(values=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) ->
  g
g

# Visualize average reactionTime by arrangement (all positions)
rbind(
  mutate(describe(dat_11[, "reactionTime"])[, c("mean", "se")], arrange="[1, 1, 1]"),
  mutate(describe(dat_12[, "reactionTime"])[, c("mean", "se")], arrange="[1, 2, 4]"),
  mutate(describe(dat_21[, "reactionTime"])[, c("mean", "se")], arrange="[2, 2, 2]"),
  mutate(describe(dat_22[, "reactionTime"])[, c("mean", "se")], arrange="[2, 4, 8]")
) %>%
  rename("averageRT"="mean") %>%
  ggplot(aes(x=arrange, y=averageRT)) +
  geom_point() +
  geom_errorbar(aes(ymin=averageRT-se*1.96, ymax=averageRT+se*1.96), size=1.2, width = 0.4, colour=c("#A1C8F5", "#F4B281", "#8DE4A0", "#F29F9A")) +
  geom_point(size=3.5, colour=c("#1244F4", "#F07A3B", "#5EC93B", "#E93732")) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  # scale_y_continuous(limits=c(2, 4.5)) +  # Imagenet  
  scale_y_continuous(limits=c(1, 2.5)) +  # Alphabet
  labs(x="Arrange", y="Average Reaction Time (sec)") ->
  g_ave
g_ave

# Quadratic trend
gs111 <- fixef(f2)[3] 
gs124 <- gs111 + fixef(f2)[8]
gs222 <- gs111 + fixef(f2)[10]
gs248 <- gs111 + fixef(f2)[12]

ls111 <- ranef(f2)$id[3] + gs111
ls124 <- ranef(f2)$id[3] + gs124
ls222 <- ranef(f2)$id[3] + gs222
ls248 <- ranef(f2)$id[3] + gs248

ls111 <- as.data.frame(ls111)
ls124 <- as.data.frame(ls124)
ls222 <- as.data.frame(ls222)
ls248 <- as.data.frame(ls248)

colnames(ls111) <- "quadratic_slope"
colnames(ls124) <- "quadratic_slope"
colnames(ls222) <- "quadratic_slope"
colnames(ls248) <- "quadratic_slope"

ls111$arrange <- "[1, 1, 1]"
ls124$arrange <- "[1, 2, 4]"
ls222$arrange <- "[2, 2, 2]"
ls248$arrange <- "[2, 4, 8]"

qslopes <- rbind(ls111, ls124, ls222, ls248)

g <- ggplot(qslopes, aes(x=arrange, y=quadratic_slope)) + geom_point() +
  stat_summary(fun.y="mean", fun.ymin=function(x)mean(x), fun.ymax=function(x)mean(x), geom="crossbar", colour="black", width=0.2) 
g
