

# Install packages
if(!require('reshape2')) install.packages('reshape2'); library(reshape2)
if(!require('devtools')) install.packages('devtools'); library(devtools)
if(!require('ggpubr')) install.packages('ggpubr'); library(ggpubr)
if(!require('cowplot')) install.packages('cowplot'); library(cowplot)
if(!require('boot')) install.packages('boot'); library(boot)
if(!require('caret')) install.packages('caret'); library(caret)
library(brms)
library(gridExtra)
library(lme4)
library(tidyverse)
library(rstan)
library(sjstats)

# Get rid of scientific notation numbers
options(scipen=999)



#......................................................................................
# CFR by age --------------------------------------------------------------
#......................................................................................

CFR.Age <- read_csv("data/cfr_age.csv")

# Need to split out data to each year of age so we can merge it back together

CFR.Age.mult.imput <- CFR.Age %>% mutate(row = 1:nrow(CFR.Age)) %>% group_by(Ref, Location, Region, Decade, row, Age.min, Age.max, Cases, Deaths) %>% do({
  Age = .$Age.min:.$Age.max
  data.frame(Age=Age, CFR= .$CFR, Cases.yr=.$Cases/.$Age.years, Deaths.yr=.$Deaths/.$Age.years)
}) %>% ungroup

CFR.Age.mult.imput$CFR <- CFR.Age.mult.imput$CFR*100

# Check that this worked correctly
sum(CFR.Age.mult.imput$Cases.yr)
sum(CFR.Age.mult.imput$Deaths.yr)

sum(CFR.Age$Cases)
sum(CFR.Age$Deaths)



# Set up age groups
CFR.Age.mult.imput$Group <- NA
CFR.Age.mult.imput$Group[which(CFR.Age.mult.imput$Age<5)] <- "0-4y"
CFR.Age.mult.imput$Group[which(CFR.Age.mult.imput$Age>=5 & CFR.Age.mult.imput$Age<20)] <- "5-19y"
CFR.Age.mult.imput$Group[which(CFR.Age.mult.imput$Age>=20)] <- ">=20y"
CFR.Age.mult.imput$Group <- factor(CFR.Age.mult.imput$Group, levels=c("0-4y", "5-19y", ">=20y"))

# Sum cases and deaths for the split ages back to the new groups & Calc CFR
CFR.Age.grouped <- CFR.Age.mult.imput %>% group_by(Ref, Location, Region, Decade, Group) %>%
  summarise(Cases.group=round(sum(Cases.yr, na.rm=TRUE),1), Deaths.group=round(sum(Deaths.yr, na.rm=TRUE),1))
CFR.Age.grouped <- as.data.frame(CFR.Age.grouped)
CFR.Age.grouped$CFR <- round(CFR.Age.grouped$Deaths.group / CFR.Age.grouped$Cases.group,2)
CFR.Age.grouped$CFR[CFR.Age.grouped$Deaths.group==0] <- 0
CFR.Age.grouped$CFR <- CFR.Age.grouped$CFR * 100
CFR.Age.grouped <- rename(CFR.Age.grouped, Cases=Cases.group, Deaths=Deaths.group)


(p2 <- ggplot(CFR.Age.mult.imput, aes(x=Group, y=CFR)) +geom_boxplot()+xlab("Age")+ylab("Case Fatality Ratio")+scale_x_discrete(labels=c("1"= "0-4", "2"="5-19", "3"="≥20")))
(p2b <- ggplot(CFR.Age.grouped, aes(x=Group, y=CFR)) +geom_boxplot()+xlab("Age")+ylab("Case Fatality Ratio"))#+scale_x_discrete(labels=c("1"= "0-4", "2"="5-19", "3"="≥20")))




# CFR by Age OR --------------------------------------------------------------

CFR.Age.grouped %>% filter(Group=="0-4y") %>% summarize(mean=mean(CFR),  ql = quantile(CFR,.025), qu = quantile(CFR,.975))
CFR.Age.grouped %>% filter(Group=="5-19y") %>% summarize(mean=mean(CFR), ql = quantile(CFR,.025), qu = quantile(CFR,.975))
CFR.Age.grouped %>% filter(Group==">=20y") %>% summarize(mean=mean(CFR),  ql = quantile(CFR,.025), qu = quantile(CFR,.975))

# Median & IQR
CFR.Age.grouped %>% filter(Group=="0-4y") %>% summarize(median=median(CFR),  ql = quantile(CFR,.25), qu = quantile(CFR,.75))
CFR.Age.grouped %>% filter(Group=="5-19y") %>% summarize(median=median(CFR), ql = quantile(CFR,.25), qu = quantile(CFR,.75))
CFR.Age.grouped %>% filter(Group==">=20y") %>% summarize(median=median(CFR),  ql = quantile(CFR,.25), qu = quantile(CFR,.75))

# add study and group weights
tot.cases <- sum(CFR.Age.grouped$Cases)
CFR.Age.grouped <- CFR.Age.grouped %>% group_by(Ref) %>% mutate(Study.Wt = sum(Cases)/tot.cases) %>% ungroup()
CFR.Age.grouped <- CFR.Age.grouped %>% group_by(Group) %>% mutate(Group.Tot = sum(Cases)) %>% ungroup()
CFR.Age.grouped <- CFR.Age.grouped %>% mutate(Group.Wt = Cases/Group.Tot,
                                              CFR.contr.grp = CFR*Group.Wt,
                                              CFR.contr.study = CFR*Study.Wt)
# Calc weighted group means
CFR.Age.grouped %>% group_by(Group) %>% summarize(sum(CFR.contr.grp, na.rm = T))
CFR.Age.grouped %>% group_by(Group) %>% summarize(sum(CFR.contr.study, na.rm = T))

# Convert to long
CFR.Age.grouped <- CFR.Age.grouped %>% mutate(row = 1:nrow(CFR.Age.grouped))
CFR.Age.long <- CFR.Age.grouped %>% group_by(row) %>%
  do({ left_join(data_frame(row=.$row, died = c(rep(0,.$Cases-.$Deaths), rep(1,.$Deaths))),.,by='row') })



mod <- glm(died ~ Group, data=CFR.Age.long, family = binomial(link = "logit"))
exp(coef(mod))[2]
exp(cbind("Odds ratio" = coef(mod), confint.default(mod, level = 0.95)))

# Calculate the probabilities
probs1 <- predict(mod, data.frame(Group=unique(CFR.Age.long$Group)), type = "response")
names(probs1) <- unique(CFR.Age.long$Group)
probs1

prob1.se <- predict(mod, data.frame(Group=unique(CFR.Age.long$Group)), se.fit=TRUE)

probs1 <- cbind(exp(prob1.se$fit) / (1+exp(prob1.se$fit)),
                exp(prob1.se$fit - 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit - 1.96*prob1.se$se.fit)),
                exp(prob1.se$fit + 1.96*prob1.se$se.fit) / (1+exp(prob1.se$fit + 1.96*prob1.se$se.fit)))
row.names(probs1) <- unique(CFR.Age.long$Group)
probs1





# Simple logistic regression to look at ORs
mod <- glm(died~Group, data=CFR.Age.long, family = binomial(link = "logit"))
summary(mod)
coef(mod)
exp(coef(mod))
inv.logit(coef(mod)[1])
inv.logit(coef(mod)[1] + coef(mod)[2])
inv.logit(coef(mod)[1] + coef(mod)[3])
exp(cbind("Odds ratio" = coef(mod), confint.default(mod, level = 0.95)))





# Model with mixed effect for Location/Decade (because some studies had multiple populations)  ---------------------------------------
# - Used log-binomial model to calculate the RR

CFR.Age.long$LocDecade <- factor(paste0(CFR.Age.long$Location, "-", CFR.Age.long$Decade))
CFR.Age.long$Group <- relevel(CFR.Age.long$Group, ref=">=20y") #set ref level

mod <- lme4::glmer(died ~ Group + (1|LocDecade), data=CFR.Age.long, family = poisson(link = "log"))
summary(mod)
exp(fixef(mod))
exp(cbind("Risk ratio" = fixef(mod), confint(mod, method="Wald", level = 0.95)[-1,]))[-1,]

# Convert to Probabilities
(ef.1 <- data.frame(effect(c("Group"), mod)))

# Check variances
get_re_var(mod) # risidual variance
re_var(mod)     # all variances
re_var(mod, adjusted = TRUE)     # all variances
r2(mod)

#prop variance between studies:
0.551 / (2.378 + 0.551)

#........................................................................................
# Vaccination & CFR - brms Estimation -------------------------------------------------------
#........................................................................................

mod_brms <- brms::brm(died ~ Group + (1|LocDecade),
                 data=CFR.Age.long,
                 family = 'bernoulli',
                 prior = set_prior('normal(0, 100)'),
                 iter = 100,
                 chains = 4,
                 cores = 4
                 )
summary(mod_brms)
exp(fixef(mod_brms))
tidy_stan(mod_brms)
icc(mod_brms, ppd=TRUE)
r2(mod_brms)

make_stancode(died ~ Group + (1|LocDecade),
              data=CFR.Age.long,
              family = 'binomial',
              prior = set_prior('normal(0, 100)'))


#........................................................................................
# Vaccination & CFR - Stan Estimation -------------------------------------------------------
#........................................................................................

dmy <- dummyVars(died ~ Group, data=CFR.Age.long)
fixed.dummies <- data.frame(predict(dmy, newdata = CFR.Age.long))

CFR.Age.stan <- list(N=as.integer(nrow(CFR.Age.long)),
                     J=as.integer(length(unique(CFR.Age.long$LocDecade))),
                     id=as.integer(CFR.Age.long$LocDecade),
                     K=3,
                     X=as.matrix(fixed.dummies),
                     Y=as.integer(CFR.Age.long$died))

table(CFR.Age.long$Ref)



# MODEL ---------------------------------------------------------------

fileName <- "source/CFR_mixed_logistic.stan"
ret <- stanc(fileName) # Check Stan file
ret_sm <- stan_model(stanc_ret = ret) # Compile Stan code
save(fileName, ret, ret_sm, file="source/CFR_mixed_logistic_compiled.Rdata")
load(file="source/CFR_mixed_logistic_compiled.Rdata")




# Run Model ---------------------------------------------------------------
rstan_options(auto_write=TRUE)
options (mc.cores=3) # Run on multiple cores

CFR.Age.stan.fit <- sampling(ret_sm, warmup=500, iter=1000, seed=123, data=CFR.Age.stan,
                             chains=3, control=list(adapt_delta=.95))

save(CFR.Age.stan.fit, file='results/CFR.Age.stan.Study.Rdata')
load(file='results/CFR.Age.stan.Study.Rdata') # CFR.Age.stan.fit

CFR.Age.stan.fit

summary(CFR.Age.stan.fit)
r2(mod_brms)




stan.fits <- rstan::extract(CFR.Age.stan.fit)

CFR0 <- stan.fits$CFR0 # 0-4y
CFR1 <- stan.fits$CFR1 # 5-19y
CFR2 <- stan.fits$CFR2 # 20+y

cfr_age <- data.frame("age1"=CFR0, "age2"=CFR1, "age3"=CFR2)
write.csv(cfr_age, file='results/cfr_age_results.csv', row.names = FALSE)

# 0-44
summary(CFR0); quantile(CFR0, probs = c(0.025, 0.975))
# 5-19y
summary(CFR1); quantile(CFR1, probs = c(0.025, 0.975))
# 20+y
summary(CFR2); quantile(CFR2, probs = c(0.025, 0.975))


# Relative Risks
RR0to1 <-  CFR1 / CFR0
summary(RR0to1); quantile(RR0to1, probs = c(0.025, 0.975))

RR0to2 <-  CFR2 / CFR0
summary(RR0to2); quantile(RR0to2, probs = c(0.025, 0.975))








