library(data.table)
library(dplyr)
library(psych)
library(stringr)
library(lme4)
library(car)
library(ggplot2)
library(ggeffects)
library(readr)
library(reghelper)
library(optimx)


files <- dir("./approved", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)
files <- dir("./approved-3sec", ".*csv$", full.names=TRUE)  # Specify data folder (alphabet or imagenet)

dat_first <- c()
dat_second <- c()

# Initialize DataFrame
dat_extracted <- data.frame(matrix(rep(NA, 5), nrow=1))[numeric(0), ]
for (i in 1: length(files)) {
  d <- fread(files[i], header=TRUE)

  if ("surveyCode" %in% colnames(d)) {
    d <- select(d, -surveyCode)
  }

  d_first <- filter(d, ActualTrials1.thisRepN>=0)  # Task1
  d_first$id <- c(rep(i, nrow(d_first)))
  d_second <- filter(d, rankImages_a.thisRepN>=0 | ActualTrials2.thisRepN>=0)  # Task2
  
  # Extract answer row
  d_first <- mutate(d_first, ans=gsub(ans_mouse.clicked_name, pattern="[^0-9\\.]", replacement="", ignore.case=TRUE))
  d_second <- mutate(d_second, ans=gsub(ans_mouse2.clicked_name, pattern="[^0-9\\.]", replacement="", ignore.case=TRUE))
  
  for (j in 1: nrow(d_first)) {
    current_col <- d_first[j,]
    current_path <- current_col$stimuli  # Current Stimulus Path
    splited_path <- str_split(current_path, pattern="/", simplify=TRUE)
    arrange <- splited_path[3]
    presented_class <- splited_path[4]
    presented_No <- parse_number(splited_path[5])

    target_path <- paste("imagenet/one_level/", presented_class, "/", presented_class, "_", presented_No, ".png", sep="")
    df_second_extracted <- filter(d_second, stimuli==target_path)
    order <- which(df_second_extracted$ans==current_col$ans)

    dat_extracted <- rbind(dat_extracted, c(presented_class, current_col$ans, order, arrange, i))
  }
  dat_first <- rbind(dat_first, d_first)
  dat_second <- rbind(dat_second, d_second)
}
colnames(dat_extracted) <- c("presented_class", "ans", "order", "arrange", "id")  # Change Row Names
dat_extracted <- mutate(dat_extracted, level=str_sub(ans, start=1, end=1))
dat_extracted <- mutate(dat_extracted, arrange_level=str_c(arrange, "_", level))
dat_extracted$order <- as.numeric(dat_extracted$order)
dat_extracted$id <- as.integer(dat_extracted$id)


# Create Graph of Cumulative Rate
dat_order <- count(dat_extracted, order, arrange)
dat_order["rate"] <- ""
dat_order <- mutate(dat_order, rate=n/count(dat_extracted)$n*2)

# Calculate Cumulative Rate
dat_order$order <- as.integer(dat_order$order)
dat_order <- dat_order[order(dat_order$order),]
dat_order <- mutate(group_by(dat_order, arrange), cumulative_rate=cumsum(rate))

dat_order %>%
  ggplot(aes(x=order, y=cumulative_rate, color=arrange)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=c(1:12))  +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits=c(-0.01, 1.01))->
  p
p


# Plot Each Order Rate
dat_order$rate <- as.double(dat_order$rate)
dat_order %>%
  ggplot(aes(x=order, y=rate, color=arrange, shape=arrange)) +
  geom_point() +
  scale_x_continuous(breaks=c(1:12)) +
  scale_y_continuous(breaks=c(0, 0.2, 0.4), limits=c(0, 0.4)) ->
  p_each_order_rate
p_each_order_rate

# Plot Each Order Rate With Errorbar
dat_errorbar <- count(dat_extracted, id, order, arrange)
dat_errorbar$order <- as.integer(dat_errorbar$order)
dat_errorbar <- dat_errorbar[order(dat_errorbar$order),]
dat_errorbar["rate"] <- ""
dat_errorbar <- mutate(dat_errorbar, rate=n/count(dat_extracted, id)[1,]$n*2)
dat_errorbar$rate <- as.double(dat_errorbar$rate)
dat_errorbar %>%
  group_by(order, arrange) %>%
  summarise(mean=mean(rate), sd=sd(rate)) ->
  dat_errorbar
dat_errorbar %>%
  ggplot(aes(x=order, y=mean, color=arrange, shape=arrange)) +
  geom_errorbar(aes(ymax=mean+sd/sqrt(56-1)*1.96, ymin=mean-sd/sqrt(56-1)*1.96)) +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1:12)) +
  scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5), limits=c(-0.05, 0.5)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Ranking", y="Selection Rate", color="Arrange", shape="Arrange") +
  scale_color_hue(labels=c(magnified="Magnified", standard="Standard")) +
  scale_shape(labels=c(magnified="Magnified", standard="Standard")) +
  guides(colour=guide_legend(reverse=TRUE), shape=guide_legend(reverse=TRUE)) ->
  p_errorbar
p_errorbar


# Create Graph of Each Level's Cumulative Rate
dat_order_each_level <- count(dat_extracted, arrange_level, arrange, level, order)
dat_order_each_level["rate"] <- 0

index_list <- count(dat_extracted, arrange_level)
for (i in 1: nrow(index_list)) {
  current_index <- index_list[i,]$arrange_level
  dat_order_each_level[dat_order_each_level$arrange_level==current_index,] <- 
    mutate(dat_order_each_level[dat_order_each_level$arrange_level==current_index,], rate=n/(index_list[i,]$n))
}
dat_order_each_level$order <- as.integer(dat_order_each_level$order)
dat_order_each_level <- dat_order_each_level[order(dat_order_each_level$order),]
dat_order_each_level <- mutate(group_by(dat_order_each_level, arrange_level), cumulative_rate=cumsum(rate))
dat_order_each_level %>%
  mutate(level=gsub(level, pattern="1", replacement="3.83")) %>%
  mutate(level=gsub(level, pattern="0", replacement="1.41")) %>%
  mutate(level=gsub(level, pattern="2", replacement="8.06")) %>%
  ggplot(aes(x=order, y=cumulative_rate, color=arrange, shape=level)) +
  geom_line() +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1:12))  +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits=c(-0.01, 1.01)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=30)) +
  theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  labs(x="Ranking", y="Cumulative Rate", color="Arrange", shape="Level") +
  scale_color_hue(labels=c(magnified="Magnified", standard="Standard")) +
  guides(colour=guide_legend(reverse=TRUE)) ->
  p_order_each_level
p_order_each_level

# Plot Each Order Rate (Each Level)
dat_order_each_level$rate <- as.double(dat_order_each_level$rate)
dat_order_each_level %>%
  mutate(level=gsub(level, pattern="1", replacement="3.83")) %>%
  mutate(level=gsub(level, pattern="0", replacement="1.41")) %>%
  mutate(level=gsub(level, pattern="2", replacement="8.06")) %>%
  ggplot(aes(x=order, y=rate, color=arrange, shape=level)) +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1:12)) +
  scale_y_continuous(breaks=c(0, 0.1, 0.2, 0.3, 0.4), limits=c(0, 0.4)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Ranking", y="Selection Rate", color="Arrange", shape="Level") +
  scale_color_hue(labels=c(magnified="Magnified", standard="Standard")) +
  guides(colour=guide_legend(reverse=TRUE)) ->
  p_each_order_rate_level
p_each_order_rate_level


# Create Graph of Each Participant's Reaction Time in Task1
dat_first[, c("reactionTime", "id", "ans", "magnified")] %>%
  mutate(level=str_sub(ans, start=1, end=1)) ->
  dat_first_extracted
dat_first_extracted$magnified <- as.factor(dat_first_extracted$magnified)
# Scatter & Line Plot
dat_first_rt <- aggregate(reactionTime~id+level+magnified, data=dat_first_extracted, FUN=mean)
dat_first_rt$level <- as.integer(dat_first_rt$level)
dat_first_rt %>%
  ggplot(aes(x=level, y=reactionTime)) +
  geom_jitter(height=0, width =0.1) +
  stat_summary(fun=mean, geom="line", size=0.5) +
  facet_wrap(~magnified) ->
  p_reactiontime
p_reactiontime

# Participants Whose Reaction Time < 6 Sec
dat_first_rt_6 <- filter(dat_first_rt, reactionTime<6)
dat_first_rt_6 %>%
  ggplot(aes(x=level, y=reactionTime)) +
  geom_jitter(height=0, width =0.1) +
  stat_summary(fun=mean, geom="line", size=0.5) +
  facet_wrap(~magnified)


# Create Graph of Each Level's Cumulative Rate (Reaction Time < 6)
dat_first_rt_6 %>%
  mutate(magnified=gsub(magnified, pattern="0", replacement="standard")) %>%
  mutate(magnified=gsub(magnified, pattern="1", replacement="magnified")) ->
  dat_first_rt_6
dat_extracted_rt_6 <- c()
for (i in 1: nrow(dat_first_rt_6)) {
  current_col <- dat_first_rt_6[i,]
  d_current <- filter(dat_extracted, id==current_col$id, level==current_col$level, arrange==current_col$magnified)
  dat_extracted_rt_6 <- rbind(dat_extracted_rt_6, d_current)
}

dat_order_each_level_rt_6 <- count(dat_extracted_rt_6, arrange_level, order)
dat_order_each_level_rt_6["rate"] <- ""

index_list <- count(dat_extracted_rt_6, arrange_level)
for (i in 1: nrow(index_list)) {
  current_index <- index_list[i,]$arrange_level
  dat_order_each_level_rt_6[dat_order_each_level_rt_6$arrange_level==current_index,] <- 
    mutate(dat_order_each_level_rt_6[dat_order_each_level_rt_6$arrange_level==current_index,], rate=n/(index_list[i,]$n))
}
dat_order_each_level_rt_6$order <- as.integer(dat_order_each_level_rt_6$order)
dat_order_each_level_rt_6 <- dat_order_each_level_rt_6[order(dat_order_each_level_rt_6$order),]
dat_order_each_level_rt_6 <- mutate(group_by(dat_order_each_level_rt_6, arrange_level), cumulative_rate=cumsum(rate))

dat_order_each_level_rt_6 %>%
  ggplot(aes(x=order, y=cumulative_rate, color=arrange_level)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=c(1:12))  +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits=c(-0.01, 1.01))->
  p_order_each_level_rt_6
p_order_each_level_rt_6

# Create Graph of Each Level's Cumulative Rate (Reaction Time >= 6)
filter(dat_first_rt, reactionTime>=6) %>%
  mutate(magnified=gsub(magnified, pattern="0", replacement="standard")) %>%
  mutate(magnified=gsub(magnified, pattern="1", replacement="magnified")) ->
  dat_first_rt_greater_6
dat_extracted_rt_greater_6 <- c()
for (i in 1: nrow(dat_first_rt_greater_6)) {
  current_col <- dat_first_rt_greater_6[i,]
  d_current <- filter(dat_extracted, id==current_col$id, level==current_col$level, arrange==current_col$magnified)
  dat_extracted_rt_greater_6 <- rbind(dat_extracted_rt_greater_6, d_current)
}

dat_order_each_level_rt_greater_6 <- count(dat_extracted_rt_greater_6, arrange_level, order)
dat_order_each_level_rt_greater_6["rate"] <- ""

index_list <- count(dat_extracted_rt_greater_6, arrange_level)
for (i in 1: nrow(index_list)) {
  current_index <- index_list[i,]$arrange_level
  dat_order_each_level_rt_greater_6[dat_order_each_level_rt_greater_6$arrange_level==current_index,] <- 
    mutate(dat_order_each_level_rt_greater_6[dat_order_each_level_rt_greater_6$arrange_level==current_index,], rate=n/(index_list[i,]$n))
}
dat_order_each_level_rt_greater_6$order <- as.integer(dat_order_each_level_rt_greater_6$order)
dat_order_each_level_rt_greater_6 <- dat_order_each_level_rt_greater_6[order(dat_order_each_level_rt_greater_6$order),]
dat_order_each_level_rt_greater_6 <- mutate(group_by(dat_order_each_level_rt_greater_6, arrange_level), cumulative_rate=cumsum(rate))

dat_order_each_level_rt_greater_6 %>%
  ggplot(aes(x=order, y=cumulative_rate, color=arrange_level)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks=c(1:12))  +
  scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1), limits=c(-0.01, 1.01))->
  p_order_each_level_rt_greater_6
p_order_each_level_rt_greater_6


# Create Graph of Rate of Choosing Each Level's Stimulus in Task1
dat_extracted %>%
  count(arrange, level, id) %>%
  mutate(rate=n/count(dat_extracted, arrange, id)[1,]$n) %>%
  group_by(level, arrange) %>%
  summarise(mean=mean(rate), sd=sd(rate)) %>%
  mutate(level=gsub(level, pattern="1", replacement="3.83")) %>%
  mutate(level=gsub(level, pattern="0", replacement="1.41")) %>%
  mutate(level=gsub(level, pattern="2", replacement="8.06")) ->
  dat_each_level
dat_each_level$level <- as.numeric(dat_each_level$level)
# Point + ErrorBar
dat_each_level %>%
  ggplot(aes(x=level, y=mean, color=arrange, shape=arrange)) +
  geom_errorbar(aes(ymax=mean+sd/sqrt(56-1)*1.96, ymin=mean-sd/sqrt(56-1)*1.96)) +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1.41, 3.83, 8.06)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=24),legend.text=element_text(size=24)) +
  labs(x="Eccentricity (deg)", y="Selection Rate", color="Arrange", shape="Arrange") ->
  p_each_level
p_each_level
# Point + Line
dat_each_level %>%
  ggplot(aes(x=level, y=mean, color=arrange, shape=arrange)) +
  geom_line() +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1.41, 3.83, 8.06)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentricity (deg)", y="Selection Rate", color="Arrange", shape="Arrange") +
  scale_color_hue(labels=c(standard="Standard", magnified="Magnified")) +
  scale_shape(labels=c(standard="Standard", magnified="Magnified")) +
  guides(shape=guide_legend(reverse=TRUE), colour=guide_legend(reverse=TRUE)) ->
  p_each_level
p_each_level

# Create Graph of Rate of Choosing Each Level's Stimulus in Task1 (reactionTime < 6)
filter(dat_first_extracted, reactionTime<6) %>%
  count(magnified, level) %>%
  mutate(rate=n/count(filter(dat_first_extracted, reactionTime<6))$n*2) ->
  dat_first_extracted_rt_6
dat_first_extracted_rt_6$magnified <- as.factor(dat_first_extracted_rt_6$magnified)
dat_first_extracted_rt_6$level <- as.numeric(dat_first_extracted_rt_6$level)
dat_first_extracted_rt_6 %>%
    ggplot(aes(x=level, y=rate, color=magnified)) +
    geom_point() +
    geom_line() ->
    p_each_level_rt_6
p_each_level_rt_6

# Create Graph of Rate of Choosing Each Level's Stimulus in Task1 (reactionTime < 6)
filter(dat_first_extracted, reactionTime>=6) %>%
  count(magnified, level) %>%
  mutate(rate=n/count(filter(dat_first_extracted, reactionTime>=6))$n*2) ->
  dat_first_extracted_rt_greater_6
dat_first_extracted_rt_greater_6$magnified <- as.factor(dat_first_extracted_rt_greater_6$magnified)
dat_first_extracted_rt_greater_6$level <- as.numeric(dat_first_extracted_rt_greater_6$level)
dat_first_extracted_rt_greater_6 %>%
  ggplot(aes(x=level, y=rate, color=magnified)) +
  geom_point() +
  geom_line() ->
  p_each_level_rt_greater_6
p_each_level_rt_greater_6


# Regression of Cumulative Curve
dat_order_each_id <- count(dat_extracted, order, arrange, level, arrange_level, id)
dat_order_each_id["rate"] <- 0

index_list <- count(dat_extracted, arrange_level, arrange, level, id)
for (i in 1: nrow(index_list)) {
  current_index <- index_list[i,]
  dat_order_each_id[dat_order_each_id$arrange_level==current_index$arrange_level & dat_order_each_id$id==current_index$id,] <- 
    mutate(dat_order_each_id[dat_order_each_id$arrange_level==current_index$arrange_level & dat_order_each_id$id==current_index$id,], rate=n/(index_list[i,]$n))

  # If order row does not exist, add it to dataframe
  for (j in 0: 12) {
    if (nrow(filter(dat_order_each_id, arrange_level==current_index$arrange_level, id==current_index$id, order==j)) == 0) {
      if (j == 0) {
        added_row <- data.frame(j, current_index$arrange, current_index$level, current_index$arrange_level, current_index$id, 0, 0)
      } else {
        prev_row <- filter(dat_order_each_id, arrange_level==current_index$arrange_level, id==current_index$id, order==j-1)
        added_row <- data.frame(j, current_index$arrange, current_index$level, current_index$arrange_level, current_index$id, 0, 0)
      }
      names(added_row) <- colnames(dat_order_each_id)
      dat_order_each_id <- rbind(dat_order_each_id, added_row)
    }
  }
}
# dat_order_each_id$order <- as.integer(dat_order_each_id$order)
dat_order_each_id <- dat_order_each_id[order(dat_order_each_id$order),]
dat_order_each_id <- mutate(group_by(dat_order_each_id, arrange_level, id), cumulative_rate=cumsum(rate))

# Plot Example of The Area of The Cumulative Rate
dat_t <- data.frame(order=c(0, 12, 12), mean=c(0, 0, 12))
dat_order_each_id %>%
  group_by(order) %>%
  filter(arrange=="magnified", level==0) %>%
  summarise(mean=mean(cumulative_rate)) %>%
    ggplot(aes(x=order, y=mean)) +
    geom_polygon(fill="lightgray") +
    geom_point(size=4) +
    geom_line() +
    scale_x_continuous(breaks=0:12, limits=c(0, 12)) +
    theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
    theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
    labs(x="Ranking", y="Cumulative Rate") ->
    p_example
p_example

# Calculate Cumulative Curve Area
dat_cc_area <- data.frame()
for (i in 1: nrow(index_list)) {
  current_index <- index_list[i,]
  filter(dat_order_each_id, arrange_level==current_index$arrange_level, id==current_index$id) %>%
    mutate(prev_cr=lag(cumulative_rate, n=1)) %>%
    mutate(area=(cumulative_rate+prev_cr)/2) ->
    current_df
  dat_cc_area <- rbind(dat_cc_area, c(current_index$arrange, current_index$level, current_index$arrange_level, current_index$id, sum(current_df$area, na.rm=TRUE)))
}
colnames(dat_cc_area) <- c("arrange", "level", "arrange_level", "id", "area")  # Change Row Names
dat_cc_area$area <- as.numeric(dat_cc_area$area)
dat_cc_area$id <- as.numeric(dat_cc_area$id)
dat_cc_area %>%
  mutate(area_normalized=(area-6)/6) %>%
  mutate(level=gsub(level, pattern=1, replacement=3.83)) %>%
  mutate(level=gsub(level, pattern=0, replacement=1.41)) %>%
  mutate(level=gsub(level, pattern=2, replacement=8.06)) ->
  dat_cc_area
dat_cc_area$level <- as.double(dat_cc_area$level)

#Raw Data Plot
dat_cc_area %>%
  group_by(level, arrange) %>%
  summarise(mean=mean(area_normalized), sd=sd(area_normalized)) ->
  dat_cc_errorbar
# Point + ErrorBar
dat_cc_errorbar %>%
  ggplot(aes(x=level, y=mean, color=arrange, shape=arrange)) +
  geom_errorbar(aes(ymax=mean+sd/sqrt(56-1)*1.96, ymin=mean-sd/sqrt(56-1)*1.96)) +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1.41, 3.83, 8.06)) +
  scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(-0.1, 1.1)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=24),legend.text=element_text(size=24)) +
  labs(x="Eccentricity (deg)", y="Area Under Cumulative Curve (Normalized)", color="Arrange", shape="Arrange")
# Point + Line
dat_cc_errorbar %>%
  ggplot(aes(x=level, y=mean, color=arrange, shape=arrange)) +
  geom_line() +
  geom_point(size=4) +
  scale_x_continuous(breaks=c(1.41, 3.83, 8.06)) +
  scale_y_continuous(breaks=seq(0, 1, by=0.2), limits=c(-0.1, 1.1)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30),legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentricity (deg)", y="Area Under Cumulative Curve (Normalized)", color="Arrange", shape="Arrange") +
  scale_color_hue(labels=c(standard="Standard", magnified="Magnified")) +
  scale_shape(labels=c(standard="Standard", magnified="Magnified")) +
  guides(shape=guide_legend(reverse=TRUE), colour=guide_legend(reverse=TRUE))

# Swapped DataFrame
dat_cc_area %>%
  mutate(arrange=gsub(arrange, pattern="magnified", replacement="tmp")) %>%
  mutate(arrange=gsub(arrange, pattern="standard", replacement="magnified")) %>%
  mutate(arrange=gsub(arrange, pattern="tmp", replacement="standard")) ->
  dat_cc_area_swap

# Simple Regression
f_cc <- lmer(area_normalized ~ arrange + level + arrange:level + (1 + arrange + level|id), data=dat_cc_area)
f_cc <- lmer(area_normalized ~ arrange + level + arrange:level + (1 + arrange + level|id), data=dat_cc_area, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))  # 3 sec
f_cc_swap <- lmer(area_normalized ~ arrange + level + arrange:level + (1 + arrange + level|id), data=dat_cc_area_swap)
plot_model(f_cc, type="pred", terms=c("level", "arrange")) +
  scale_x_continuous(breaks=c(1.41,3.83,8.06), limits=c(1.2,8.2)) +
  theme(axis.title.x=element_text(size=36), axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
  theme(legend.position=c(1,1), legend.justification=c(1,1)) +
  labs(x="Eccentricity (deg)", y="Area Under Cumulative Curve (Normalized)", color="Arrange") +
  guides(colour=guide_legend(reverse=TRUE))

# Polynomial Regression
f_cc_q <- lmer(area_normalized ~ poly(level, 2, raw=TRUE) + arrange + poly(level, 2, raw=TRUE):arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_cc_area)
f_cc_q <- lmer(area_normalized ~ poly(level, 2, raw=TRUE) + arrange + poly(level, 2, raw=TRUE):arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_cc_area, control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
f_cc_q_swap <- lmer(area_normalized ~ poly(level, 2, raw=TRUE) + arrange + poly(level, 2, raw=TRUE):arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_cc_area_swap)
plot_model(f_cc_q, type="pred", terms=c("level [all]", "arrange"))

# Simple Slopes
arrange <- c("Magnified", "Standard")
data.frame(rbind(summary(f_cc)$coefficients[3,], summary(f_cc_swap)$coefficients[3,])) %>%
  cbind(., arrange) ->
  simple_coefficients_cc
simple_coefficients_cc %>%
  ggplot(aes(x=arrange, y=Estimate, color=arrange)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymax=Estimate+`Std..Error`*1.96, ymin=Estimate-`Std..Error`*1.96), width=0.2) +
  scale_y_continuous(breaks=seq(-0.05, 0.05, length=5), limits=c(-0.05, 0.05)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.position="none") +
  labs(x="Arrange", y="First-Order Coefficient") +
  scale_x_discrete(limit=c('Standard', 'Magnified')) ->
  p_simple_slopes_cc
p_simple_slopes_cc

# Quadratic Slopes
data.frame(rbind(summary(f_cc_q)$coefficients[3,], summary(f_cc_q_swap)$coefficients[3,])) %>%
  cbind(., arrange) ->
  quadratic_coefficients_cc
quadratic_coefficients_cc %>%
  ggplot(aes(x=arrange, y=Estimate, color=arrange)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymax=Estimate+`Std..Error`*1.96, ymin=Estimate-`Std..Error`*1.96), width=0.2) +
  scale_y_continuous(breaks=seq(-0.05, 0.05, length=5), limits=c(-0.05, 0.05)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.position="none") +
  labs(x="Arrange", y="Second-Order Coefficient") +
  scale_x_discrete(limit=c('Standard', 'Magnified')) ->
  p_quadratic_slopes_cc
p_quadratic_slopes_cc


# Regression of Rate of Choosing Each Level's Stimulus in Task1
dat_logistic <- expand.grid(arrange=c("magnified", "standard"), level=c(0, 1, 2), id=c(1:57))
dat_logistic["n"] <- 0
index_list <- count(dat_extracted, arrange, level, id)
for (i in 1:nrow(index_list)) {
  cur_index <- index_list[i,]
  dat_logistic[dat_logistic$arrange==cur_index$arrange & dat_logistic$level==cur_index$level & dat_logistic$id==cur_index$id,]$n <-
    cur_index$n
}
dat_logistic %>%
  mutate(level=gsub(level, pattern=1, replacement=3.83)) %>%
  mutate(level=gsub(level, pattern=0, replacement=1.41)) %>%
  mutate(level=gsub(level, pattern=2, replacement=8.06)) ->
  dat_logistic
dat_logistic$level <- as.numeric(dat_logistic$level)

dat_logistic %>%
  mutate(arrange=gsub(arrange, pattern="magnified", replacement="tmp")) %>%
  mutate(arrange=gsub(arrange, pattern="standard", replacement="magnified")) %>%
  mutate(arrange=gsub(arrange, pattern="tmp", replacement="standard")) ->
  dat_logistic_swap
  
# Simple Regression
f_logistic <- glmer(cbind(n, 80-n) ~ arrange + level + arrange:level + (1 + arrange + level|id), data=dat_logistic, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
f_logistic_swap <- glmer(cbind(n, 80-n) ~ arrange + level + arrange:level + (1 + arrange + level|id), data=dat_logistic_swap, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
# No Interaction
f_logistic_ni <- glmer(cbind(n, 80-n) ~ arrange + level + (1 + arrange + level|id), data=dat_logistic, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
plot_model(f_logistic, type="pred", terms=c("level", "arrange"))

# Polynomial Regression
f_logistic_q <- glmer(cbind(n, 80-n) ~ poly(level, 2, raw=TRUE) + arrange + poly(level, 2, raw=TRUE):arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_logistic, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
f_logistic_q_swap <- glmer(cbind(n, 80-n) ~ poly(level, 2, raw=TRUE) + arrange + poly(level, 2, raw=TRUE):arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_logistic_swap, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
# No Interaction
f_logistic_q_ni <- glmer(cbind(n, 80-n) ~ poly(level, 2, raw=TRUE) + arrange + (1 + poly(level, 2, raw=TRUE) + arrange|id), data=dat_logistic, family="binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
plot_model(f_logistic_q, type="pred", terms=c("level [all]", "arrange")) +
  scale_x_continuous(breaks=c(1.41,3.83,8.06), limits=c(1.2,8.2)) +
  theme(axis.title.x=element_text(size=36), axis.title.y=element_text(size=30)) +
  theme(axis.text.x=element_text(size=30), axis.text.y=element_text(size=30)) +
  theme(legend.title=element_text(size=30), legend.text=element_text(size=30)) +
  theme(legend.position=c(1,0), legend.justification=c(1,0)) +
  labs(x="Eccentricity (deg)", y="Selection Rate", color="Arrange") +
  guides(colour=guide_legend(reverse=TRUE))

# Simple Slopes
arrange <- c("Magnified", "Standard")
data.frame(rbind(summary(f_logistic)$coefficients[3,], summary(f_logistic_swap)$coefficients[3,])) %>%
  cbind(., arrange) ->
  simple_coefficients
simple_coefficients %>%
  ggplot(aes(x=arrange, y=Estimate, color=arrange, shape=arrange)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymax=Estimate+`Std..Error`*1.96, ymin=Estimate-`Std..Error`*1.96), width=0.2) +
  scale_y_continuous(breaks=seq(-0.6, 0.6, length=7), limits=c(-0.6, 0.6)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30)) +
  theme(legend.position="none") +
  labs(x="Arrange", y="First-Order Coefficient") +
  scale_x_discrete(limit=c('Standard', 'Magnified')) ->
  p_simple_slopes
p_simple_slopes

# Quadratic Slopes
data.frame(rbind(summary(f_logistic_q)$coefficients[3,], summary(f_logistic_q_swap)$coefficients[3,])) %>%
  cbind(., arrange) ->
  quadratic_coefficients
quadratic_coefficients %>%
  ggplot(aes(x=arrange, y=Estimate, color=arrange, shape=arrange)) +
  geom_point(size=4) +
  geom_errorbar(aes(ymax=Estimate+`Std..Error`*1.96, ymin=Estimate-`Std..Error`*1.96), width=0.2) +
  scale_y_continuous(breaks=seq(-0.1, 0.1, length=11), limits=c(-0.08, 0.08)) +
  theme(axis.title.x=element_text(size=36),axis.title.y=element_text(size=36)) +
  theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30))  +
  theme(legend.position="none") +
  labs(x="Arrange", y="Second-Order Coefficient") +
  scale_x_discrete(limit=c('Standard', 'Magnified')) ->
  p_quadratic_slopes
p_quadratic_slopes
