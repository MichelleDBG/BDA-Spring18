require(lme4)
require(lattice)

spider <- read.csv("P:/My Documents/BDA_Spring2018/spider_wolff_gorb_2013.csv")
spider$group <- factor(paste0(spider$leg,spider$type))
g.mus <- tapply(spider$friction,spider$group,mean)
boxplot(spider$friction ~ spider$type + spider$leg, col = c("grey90","grey40"), las=2, outline=F,
        main = "Friction distributions for push vs. pull across different leg sets")
stripchart(friction ~ type + leg, data = spider, method = "jitter", pch = 16, vertical = T, ylab = "",
           ylim = c(0, 2), xlim = c(0.5, 4.5), las = 2, cex=0.55, add=T)
points(1:8,g.mus,pch=16,cex=1.25,col="royalblue")

# fit models
f <- lm(friction ~ type + leg, data = spider) #5 paramters 1: intercept, 1 type, 3 for leg
f.i <- lm(friction ~ type + leg + type:leg, data=spider) # 8 parameters
fme <- lmer(friction ~ type + leg + (type | leg),spider) # only the variation around leg for pull/push
summary(lmer(friction ~ (leg | type), spider))
summary(lmer(friction ~ (type| leg ), spider))

# explore summaries
summary(f)
summary(f.i)
summary(fme) # get random effects, and fixed effects
# in random effects get standard deviation around the variance around the variance of the random effects.
# fixed effects: then might change coefficients, but still trying to minimize residuals, in y directions


# check assumptions - normality of residuals
qqnorm(resid(f), main="lm without interaction")
qqline(resid(f))

#Do we want distributions on these paramters or no? Do we get mroe with more parameters or cost?
#some skew
qqnorm(resid(f.i), main="lm with interaction")
qqline(resid(f.i))

qqnorm(resid(fme), main="mixed effects")
qqline(resid(fme))

hist(resid(f))
hist(resid(f.i))
hist(resid(fme))

# check for homogeneity of variance
plot(fitted.values(f),resid(f),main="lm",xlab="fitted",ylab="residuals")
abline(h=0)

plot(fitted.values(f.i),resid(f.i),main="lmi",xlab="fitted",ylab="residuals")
abline(h=0)

plot(fitted.values(fme),resid(fme),main="lme",xlab="fitted",ylab="residuals")
abline(h=0)

###########################################################################################
data(sleepstudy)

# linear model
lm1 <- lm(Reaction ~ Days, sleepstudy)
plot(sleepstudy$Days, sleepstudy$Reaction, xlab="Days no sleep", 
     ylab="Reaction Time (ms)")
abline(lm1)

summary(lm1)

qqnorm(resid(lm1), main = "linear model")
qqline(resid(lm1))
hist(resid(lm1))
plot(fitted.values(lm1),resid(lm1),main="",xlab="fitted",ylab="residuals")
abline(h=0)

xyplot(Reaction ~ Days | Subject, sleepstudy, type = c("g","p","r"),
       index = function(x,y) coef(lm(y ~ x))[1],
       xlab = "Days of sleep deprivation",
       ylab = "Average reaction time (ms)", aspect = "xy")

# simple variance components appraoch
fmday <- lmer(Reaction ~ (1|Days), sleepstudy) # days across subject
summary(fmday)
dotplot(ranef(fmday,condVar=T))
 #show distribution of variance, from grand mean, offset to day, then allow varaition around dot mean of that day

#There is variation, how much of the variance is accounted for by day, the rest is by person, stuff we didnt' measure
# Get first componant divided by sum of all variance. 
# the first part of random effects variance over all variation in reaction time (over all subjects)

# see variance components
(vday <- as.data.frame(VarCorr(fmday)))
vday$vcov[1]/sum(vday$vcov)

fmsub <- lmer(Reaction ~ (1|Subject), sleepstudy) # subject across days
summary(fmsub)
dotplot(ranef(fmsub,condVar=T))

# When look by subject, all change in reaction time per person, the rest is left over
# multiple repeat measuers per person is somewhat meaningless. 
# see variance components
(vsub <- as.data.frame(VarCorr(fmsub)))
vsub$vcov[1]/sum(vsub$vcov)

# mixed effect model approach, random intercepts for subjects
# fit slope by days but allow residuals to vary by subject, randomly move intercept 
# based on variation subject o subject. 
fm1 <- lmer(Reaction ~ Days + (1|Subject), data=sleepstudy) # random intercept
summary(fm1)
summary(lm1)
dotplot(ranef(fm1,condVar=T)) # 95% CI for intercepts

# check assumptions
qqnorm(resid(lm1),main="lm")
qqline(resid(lm1))
qqnorm(resid(fm1),main="lme-rand int")
qqline(resid(fm1))
hist(resid(lm1),main="lm",xlim=c(-150,150))
hist(resid(fm1),main="lme-rand int",xlim=c(-150,150))

lm.i <- lm(Reaction~Days*Subject,data=sleepstudy)
qqnorm(resid(lm.i),main="lme-interaction")
qqline(resid(lm.i))
hist(resid(lm.i),main="lm-interaction",xlim=c(-150,150))
plot(fitted.values(lm1),resid(lm1),main="lm",xlab="fitted",
     ylab="residuals")
abline(h=0)
plot(fitted.values(fm1),resid(fm1),main="lm-random int",
     xlab="fitted",ylab="residuals")
abline(h=0)
plot(fitted.values(lm.i),resid(fm1),main="lm-interaction",
     xlab="fitted",ylab="residuals")
abline(h=0)

# now not have 2 additional paramters for each subject, let be random,
#random effect of slope by subject days conditioned on subject
# mixed effect model, random intercepts and slopes for subjects
fm2 <- lmer(Reaction ~ Days + (Days|Subject), data=sleepstudy)
# random intercept and slope
summary(fm2)
summary(lm1)
dotplot(ranef(fm2,condVar=T)) #95% CI for slopes and intercepts
#how slope and intercept change, vary per subject. variation, residuals, now far fewer
#parameters - if interested in those individuals then interaction, but if want to generalize
# to another set of plots, then let random effects, then don't overfit. 

qqnorm(resid(lm1),main="lm")
qqline(resid(lm1))
qqnorm(resid(fm1),main="lme-rand int")
qqline(resid(fm1))
qqnorm(resid(fm2),main="lme-rand slope+int")
qqline(resid(fm2))
qqnorm(resid(lm.i),main="lm-interactions")
qqline(resid(lm.i))
hist(resid(lm1),main="lm",xlim=c(-150,150))
hist(resid(fm1),main="lme-rand int",xlim=c(-150,150))
hist(resid(fm2),main="lme-rand slope+int",xlim=c(-150,150))
hist(resid(lm.i),main="lm-interaction",xlim=c(-150,150))

plot(fitted.values(lm1),resid(lm1),main="lm",
     xlab="fitted",ylab="residuals")
abline(h=0)
plot(fitted.values(lm.i),resid(lm.i),main="lm-interaction",
     xlab="fitted",ylab="residuals")
abline(h=0)
plot(fitted.values(fm1),resid(fm1),main="lm-random int",
     xlab="fitted",ylab="residuals")
abline(h=0)
plot(fitted.values(fm2),resid(fm2),main="lm-random slope+int",
     xlab="fitted",ylab="residuals")
abline(h=0)

####################################################################

(v1 <- as.data.frame(VarCorr(fm1)))
(v2 <- as.data.frame(VarCorr(fm2)))

sum(vday$vcov)
sum(vsub$vcov)
sum(v1$vcov)
sum(v2$vcov)
