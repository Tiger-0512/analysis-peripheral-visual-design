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


files <- dir("./results/imagenet/approved", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)
dat <- c()
for (i in 1: length(files)){

  d <- fread(files[i], header=TRUE)
  d <- d[, c("TF", "size", "rate", "pos", "ori", "reactionTime")]

  d$stimSize <- d$size * d$rate ^ d$pos
  d$invS <- 1 / (pi * (d$stimSize / 2) ^ 2)
  d$size_1 <- d$size
  d$size_2 <- d$size * d$rate
  d$size_3 <- d$size * d$rate * d$rate
  d$size <- as.character(d$size)
  d$rate <- as.character(d$rate)
  d$pos <- as.character(d$pos)
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

  dat <- rbind(dat, d)
}
dat <- dat[dat$TF=="TRUE",]  # Extract TF==TRUE
# When using imagenet results, execute below one line
dat <- subset(dat, id != 5 & id != 16 & id != 19 & id != 45)  # Extract participants who's median of the reactionTime < 10
dat <- na.omit(dat)  # Fill NA

dat_2 <- aggregate(reactionTime ~ arrange + pos + stimSize, dat, mean)
dat_2$pos <- as.character(dat_2$pos)

dat_2 %>%
  ggplot(aes(x=stimSize, y=reactionTime, color=pos, shape=arrange)) +
  geom_point()

dat_3 <- aggregate(reactionTime ~ pos + stimSize, dat, mean)
dat_3$pos <- as.character(dat_3$pos)

dat_3 %>%
  ggplot(aes(x=stimSize, y=reactionTime, color=pos)) +
  geom_line() +
  scale_color_manual(values=c("#E4211D", "#377EB7", "#4DB04B")) +
  scale_y_continuous(limits=c(1, 5.5))


# Regression
f1 <- glmer(reactionTime ~ stimSize + pos + (1 + stimSize + pos|id), data=dat, family=inverse.gaussian(link="identity"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
f2 <- glmer(reactionTime ~ invS + pos + invS:pos + (1 + invS + pos|id), data=dat, family=inverse.gaussian(link="log"), control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))


plot_model(f1, type="pred", terms=c("stimSize", "pos"))
plot_model(f2, type="pred", terms=c("invS", "pos")) +
  scale_x_continuous(breaks=c(0.125,0.25,0.5,1.00,2.00))
plot_model(f2, type="pred", terms=c("invS", "pos"), colors=c("#1245F5", "#F17B3B", "#5FC93A")) +  # From seaborn "bright" color palette
  scale_x_continuous(name="S", trans="reverse", breaks=c(0.125,0.25,0.5,1.00,2.00), labels=c(8,4,2,1,0.5))


# log(μ) = β0 + β1x1(S^-1) + β2x2(pos=3.83, 0or1) + β3x3(pos=8.06, 0or1) + β4x1x2 + β5x1x3
plot(function(x1) exp(1.018450 + 0.231975 * x1), xlim=c(0,1.5), ylim=c(2.5,4.5), col=2, xlab="invS", ylab="reactionTime")
plot(function(x2) exp(1.010612 + 0.306478 * x2), xlim=c(0,1.5), col=4, add=TRUE)
plot(function(x3) exp(1.192104 + 0.167547 * x3), xlim=c(0,1.5), col=3, add=TRUE)
plot(function(x) 3.319 + 0 * x, xlim=c(0,1.5), col=1, add=TRUE)

