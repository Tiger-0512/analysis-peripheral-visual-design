library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(lme4)
library(car)
library(ggplot2)
library(sjPlot)
library(ggeffects)
library(reghelper)
library(gridExtra)


files <- dir("../approved", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)
dat <- c()

for (i in 1: length(files)) {
  d <- fread(files[i], header=TRUE)
  d <- filter(d, ActualTrials.thisTrialN>=0)
  d <- d[, c("TF", "size_0", "size_1", "size_2", "pos", "ori", "reactionTime")]
  
  d$size_0 <- as.character(d$size_0)
  d$size_1 <- as.character(d$size_1)
  d$size_2 <- as.character(d$size_2)
  d$pos <- as.factor(d$pos)
  d$arrange <- paste(d$size_0, d$size_1, d$size_2, c("", ""))
  d$id <- i  # Glmer requires numeric variable for subject id

  d %>%
    mutate(pos=gsub(pos, pattern="1", replacement="3.83")) %>%
    mutate(pos=gsub(pos, pattern="0", replacement="1.41")) %>%
    mutate(pos=gsub(pos, pattern="2", replacement="8.06")) ->
    d
  d$pos <- as.numeric(d$pos)
  
  dat <- rbind(dat, d)
}
dat <- subset(dat, id!=21 & id!=23 & id!=24)  # Extract participants who's median of the reactionTime < 10

filter(count(dat, TF, id), TF==TRUE) %>% 
  mutate(Acc=n/324) ->
  average_acc
average_acc <- mean(average_acc$Acc)
average_rt <- mean(tapply(dat$reactionTime, dat$id, mean))

dat <- dat[dat$TF=="TRUE",]  # Extract TF==TRUE

# f2 <- glmer(reactionTime ~ poly(pos, 2, raw = TRUE) + arrange + poly(pos, 2, raw = TRUE):arrange + (1 + poly(pos, 2, raw = TRUE) + arrange|id), data = dat, family = Gamma(), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))  # Polynomial model
# f2 <- glmer(reactionTime ~ poly(pos, 2, raw = TRUE) + arrange + poly(pos, 2, raw = TRUE):arrange + (1 + poly(pos, 2, raw = TRUE) + arrange|id), data = dat, family = inverse.gaussian(), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))  # Polynomial model

# Show the average reaction time of each participant in each arrangement
basic_sizes <- c(1.28, 1.44, 5.34)
to_rate <- c(2/3, 1, 3/2)
Arrange <- c()
for (i in 1: 3) {
  for (j in 1: 3) {
    for (k in 1: 3) {
      str <- paste(
        round(basic_sizes[1]*to_rate[i], digits=2), 
        round(basic_sizes[2]*to_rate[j], digits=2),
        round(basic_sizes[3]*to_rate[k], digits=2)
      )
      Arrange <- c(Arrange, str)
    }
  }
}

aggregate(reactionTime ~ pos+arrange, dat, mean) %>%
  cbind(., sort(rep(Arrange, times=3))) ->
  class_mean
colnames(class_mean) <- c("pos", "arrange", "reactionTime", "Arrange")
line_0 <- ggplot(class_mean[1:27,], aes(x=pos, y=reactionTime, color=Arrange)) +
  geom_line() +
  geom_point(size=4) +
  scale_y_continuous(limits=c(2.0, 4.0)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=28)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentircity (deg)", y="Reaction Time (sec)", color="Arrange")
line_1 <- ggplot(class_mean[28:54,], aes(x=pos, y=reactionTime, color=Arrange)) +
  geom_line() +
  geom_point(size=4) +
  scale_y_continuous(limits=c(2.0, 4.0)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=28)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentircity (deg)", y="Reaction Time (sec)", color="Arrange")
line_2 <- ggplot(class_mean[55:81,], aes(x=pos, y=reactionTime, color=Arrange)) +
  geom_line() +
  geom_point(size=4) +
  scale_y_continuous(limits=c(2.0, 4.0)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=28)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentircity (deg)", y="Reaction Time (sec)", color="Arrange")
gridExtra::grid.arrange(line_0,line_1,line_2, ncol=3)

dat <- mutate(dat, pos=pos/10)
dat <- mutate(dat, reactionTime=reactionTime/10)
dat <- na.omit(dat)  # Fill NA

# Regression of each arrange
df_iter <- expand.grid(0:2, 0:2, 0:2)
df_iter$Var1 <- as.character(df_iter$Var1)
df_iter$Var2 <- as.character(df_iter$Var2)
df_iter$Var3 <- as.character(df_iter$Var3)
df_iter$arrange <- paste(df_iter$Var1, df_iter$Var2, df_iter$Var3, c("", ""))
df_iter <- df_iter[order(df_iter$arrange),]

# First-order regression
dat_coefficients <- data.frame()
for(i in 1: nrow(df_iter)) {
  cur_df <- subset(dat, arrange==df_iter[i,]$arrange)
  nam <- paste("plot_",i,sep="")
  f <- glmer(reactionTime ~ pos + (1 + pos|id), data=cur_df, family=Gamma(), control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10e5)))
  assign(nam, f)
  
  coefficient <- summary(f)$coefficients[2,]
  coefficient["arrange"] = df_iter[i,]$arrange
  dat_coefficients <- rbind(dat_coefficients, coefficient)
}
# Second-order regression
dat_coefficients_q <- data.frame()
for(i in 1: nrow(df_iter)) {
  cur_df <- subset(dat, arrange==df_iter[i,4])
  nam_quad <- paste("plot_quad_",i,sep="")
  f_q <- glmer(reactionTime ~ poly(pos, 2, raw=TRUE) + (1 + poly(pos, 2, raw=TRUE)|id), data=cur_df, family=Gamma(), control=glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 2e7)))
  assign(nam_quad, f_q)
  
  coefficient <- summary(f_q)$coefficients[3,]
  coefficient["arrange"] = df_iter[i,]$arrange
  dat_coefficients_q <- rbind(dat_coefficients_q, coefficient)
}

# Plot first-order regression
grid.arrange(
  plot_model(plot_1, type="pred", terms="pos") + ggtitle("0.851, 0.958, 3.556") + ylim(0, 0.4),
  plot_model(plot_2, type="pred", terms="pos") + ggtitle("100") + ylim(0, 0.4),
  plot_model(plot_3, type="pred", terms="pos") + ggtitle("200") + ylim(0, 0.4),
  plot_model(plot_4, type="pred", terms="pos") + ggtitle("010") + ylim(0, 0.4),
  plot_model(plot_5, type="pred", terms="pos") + ggtitle("110") + ylim(0, 0.4),
  plot_model(plot_6, type="pred", terms="pos") + ggtitle("210") + ylim(0, 0.4),
  plot_model(plot_7, type="pred", terms="pos") + ggtitle("020") + ylim(0, 0.4),
  plot_model(plot_8, type="pred", terms="pos") + ggtitle("120") + ylim(0, 0.4),
  plot_model(plot_9, type="pred", terms="pos") + ggtitle("220") + ylim(0, 0.4)
)
grid.arrange(
  plot_model(plot_10, type="pred", terms="pos") + ggtitle("001") + ylim(0, 0.4),
  plot_model(plot_11, type="pred", terms="pos") + ggtitle("101") + ylim(0, 0.4),
  plot_model(plot_12, type="pred", terms="pos") + ggtitle("201") + ylim(0, 0.4),
  plot_model(plot_13, type="pred", terms="pos") + ggtitle("011") + ylim(0, 0.4),
  plot_model(plot_14, type="pred", terms="pos") + ggtitle("111") + ylim(0, 0.4),
  plot_model(plot_15, type="pred", terms="pos") + ggtitle("211") + ylim(0, 0.4),
  plot_model(plot_16, type="pred", terms="pos") + ggtitle("021") + ylim(0, 0.4),
  plot_model(plot_17, type="pred", terms="pos") + ggtitle("121") + ylim(0, 0.4),
  plot_model(plot_18, type="pred", terms="pos") + ggtitle("221") + ylim(0, 0.4)
)
grid.arrange(
  plot_model(plot_19, type="pred", terms="pos") + ggtitle("002") + ylim(0, 0.4),
  plot_model(plot_20, type="pred", terms="pos") + ggtitle("102") + ylim(0, 0.4),
  plot_model(plot_21, type="pred", terms="pos") + ggtitle("202") + ylim(0, 0.4),
  plot_model(plot_22, type="pred", terms="pos") + ggtitle("012") + ylim(0, 0.4),
  plot_model(plot_23, type="pred", terms="pos") + ggtitle("112") + ylim(0, 0.4),
  plot_model(plot_24, type="pred", terms="pos") + ggtitle("212") + ylim(0, 0.4),
  plot_model(plot_25, type="pred", terms="pos") + ggtitle("022") + ylim(0, 0.4),
  plot_model(plot_26, type="pred", terms="pos") + ggtitle("122") + ylim(0, 0.4),
  plot_model(plot_27, type="pred", terms="pos") + ggtitle("222") + ylim(0, 0.4)
)
# Plot second-order regression
grid.arrange(
  plot_model(plot_quad_1, type="pred", terms="pos [all]") + ggtitle("000") + ylim(0, 0.4),
  plot_model(plot_quad_2, type="pred", terms="pos [all]") + ggtitle("100") + ylim(0, 0.4),
  plot_model(plot_quad_3, type="pred", terms="pos [all]") + ggtitle("200") + ylim(0, 0.4),
  plot_model(plot_quad_4, type="pred", terms="pos [all]") + ggtitle("010") + ylim(0, 0.4),
  plot_model(plot_quad_5, type="pred", terms="pos [all]") + ggtitle("110") + ylim(0, 0.4),
  plot_model(plot_quad_6, type="pred", terms="pos [all]") + ggtitle("210") + ylim(0, 0.4),
  plot_model(plot_quad_7, type="pred", terms="pos [all]") + ggtitle("020") + ylim(0, 0.4),
  plot_model(plot_quad_8, type="pred", terms="pos [all]") + ggtitle("120") + ylim(0, 0.4),
  plot_model(plot_quad_9, type="pred", terms="pos [all]") + ggtitle("220") + ylim(0, 0.4)
)
grid.arrange(
  plot_model(plot_quad_10, type="pred", terms="pos [all]") + ggtitle("001") + ylim(0, 0.4),
  plot_model(plot_quad_11, type="pred", terms="pos [all]") + ggtitle("101") + ylim(0, 0.4),
  plot_model(plot_quad_12, type="pred", terms="pos [all]") + ggtitle("201") + ylim(0, 0.4),
  plot_model(plot_quad_13, type="pred", terms="pos [all]") + ggtitle("011") + ylim(0, 0.4),
  plot_model(plot_quad_14, type="pred", terms="pos [all]") + ggtitle("111") + ylim(0, 0.4),
  plot_model(plot_quad_15, type="pred", terms="pos [all]") + ggtitle("211") + ylim(0, 0.4),
  plot_model(plot_quad_16, type="pred", terms="pos [all]") + ggtitle("021") + ylim(0, 0.4),
  plot_model(plot_quad_17, type="pred", terms="pos [all]") + ggtitle("121") + ylim(0, 0.4),
  plot_model(plot_quad_18, type="pred", terms="pos [all]") + ggtitle("221") + ylim(0, 0.4)
)
grid.arrange(
  plot_model(plot_quad_19, type="pred", terms="pos [all]") + ggtitle("002") + ylim(0, 0.4),
  plot_model(plot_quad_20, type="pred", terms="pos [all]") + ggtitle("102") + ylim(0, 0.4),
  plot_model(plot_quad_21, type="pred", terms="pos [all]") + ggtitle("202") + ylim(0, 0.4),
  plot_model(plot_quad_22, type="pred", terms="pos [all]") + ggtitle("012") + ylim(0, 0.4),
  plot_model(plot_quad_23, type="pred", terms="pos [all]") + ggtitle("112") + ylim(0, 0.4),
  plot_model(plot_quad_24, type="pred", terms="pos [all]") + ggtitle("212") + ylim(0, 0.4),
  plot_model(plot_quad_25, type="pred", terms="pos [all]") + ggtitle("022") + ylim(0, 0.4),
  plot_model(plot_quad_26, type="pred", terms="pos [all]") + ggtitle("122") + ylim(0, 0.4),
  plot_model(plot_quad_27, type="pred", terms="pos [all]") + ggtitle("222") + ylim(0, 0.4)
)

# f <- glmer(
#   reactionTime ~ size_0 + size_1 + size_2 + size_0:size_1 + size_1:size_2 + size_2:size_0 + (1 + size_0 + size_1 + size_2|id),
#   data=dat,
#   family=inverse.gaussian(link="identity"),
#   control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10e7))
# )

# summary(f1)  # Wald tests for parameter estimates
# Anova(f1)  # Type II Wald chisquare test by car package
# plot(f1)  # Residual plot
# ranef(f1)  # Distribution of individual estimates (as deviation from the global estimate)


# Visualize simple slopes
colnames(dat_coefficients) <- c("First-Order Coefficient", "Std.Error", "t value", "Pr(>|z|)", "arrange")
dat_coefficients$`First-Order Coefficient`<- as.numeric(dat_coefficients$`First-Order Coefficient`)
dat_coefficients$Std.Error <- as.numeric(dat_coefficients$Std.Error)

dat_coefficients %>%
  cbind(., Arrange) %>%
  mutate(`First-Order Coefficient`=`First-Order Coefficient`/100) %>%
  mutate(Std.Error=Std.Error/100) %>%
  ggplot(aes(x=Arrange, y=`First-Order Coefficient`)) +
  geom_errorbar(aes(ymin=`First-Order Coefficient`-Std.Error*1.96, ymax=`First-Order Coefficient`+Std.Error*1.96)) +
  geom_point(size=4) +
  scale_y_continuous(limits=c(-0.05, 0.05)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=20, angle=70, hjust=1),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=24),legend.text=element_text(size=24)) ->
  g_simple_slope
g_simple_slope


# Visualize quadratic slopes
colnames(dat_coefficients_q) <- c("Second-Order Coefficient", "Std.Error", "t value", "Pr(>|z|)", "arrange")
dat_coefficients_q$`Second-Order Coefficient` <- as.numeric(dat_coefficients_q$`Second-Order Coefficient`)
dat_coefficients_q$Std.Error <- as.numeric(dat_coefficients_q$Std.Error)
dat_coefficients_q <- cbind(dat_coefficients_q, Arrange)

dat_coefficients_q %>%
  mutate(`Second-Order Coefficient`=`Second-Order Coefficient`/1000) %>%
  mutate(Std.Error=Std.Error/1000) %>%
  ggplot(aes(x=Arrange, y=`Second-Order Coefficient`)) +
  geom_errorbar(aes(ymin=`Second-Order Coefficient`-Std.Error*1.96, ymax=`Second-Order Coefficient`+Std.Error*1.96)) +
  geom_point(size=4) +
  scale_y_continuous(limits=c(-0.03, 0.03)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=20, angle=70, hjust=1),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=24),legend.text=element_text(size=24)) ->
  g_quadratic_slope
g_quadratic_slope
