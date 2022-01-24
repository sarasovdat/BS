########################################
# 3. domača naloga : Bayesova statistika
# Sara Bizjak, 27202020
########################################


########################################

# 1. NALOGA

# Uvoz podatkov:
setwd("~/Documents/IŠRM magistrski/Bayesova statistika/BS/Domače naloge/DN3")
source("podatki_sole.R")
str(pod)

library(ggplot2)
ggplot(pod, aes(x = school, y = mathscore, group = school)) +
        stat_summary(fun.ymin = min, fun.ymax = max, fun.y = mean) + 
        labs(title = "Povprecja (pika) in razpon rezultatov po solah")

# Preureditev podatkov: 
library(dplyr)
pod.sole = pod %>%
  group_by(school) %>% 
  summarise(povprecje = mean(mathscore),
                        n=length(mathscore), 
                        varianca = var(mathscore))


# Nastavitev parametrov:

# Parametri (hiper)apriornih porazdelitev, isti kot na vajah:
sigma20 = 100
nu0 = 1
eta20 = 100
kappa0 = 1
mu0 = 50
tau20 = 25

# Parametri iz domače naloge: 
a = 2
b = 1 / 10
alpha = 2
k.max = 1000

### Pripravimo si kolicine, ki jih bomo potrebovali iz podatkov
x = pod
m = length(pod.sole$school) 
n = pod.sole$n
x.povpr = pod.sole$povprecje
x.var = pod.sole$varianca

### Dolocimo si zacetne vrednosti
muGroups = x.povpr
mu = mean(muGroups)
eta2 = var(muGroups)

# DODANO: shranimo sigmo za vsako skupino ločeno
sigma2Groups = x.var 

### Pripravimo si prostor za shranjevanje
n.iter = 5000

muGroups.all = matrix(nrow = n.iter, ncol = m) 
mu.all = rep(NA, n.iter)
eta2.all = rep(NA, n.iter)

# DODANO:
sigma2Groups.all = matrix(nrow=n.iter, ncol= m)
sigma20.all = rep(NA, n.iter)
nu0.all = rep(NA, n.iter)

### Na prvo mesto si shranimo zacetne vrednosti (nepotrebno)
muGroups.all[1, ] = muGroups 
mu.all[1] = mu
eta2.all[1] = eta2

# DODANO: shranimo zacetne vrednosti novih parametrov
sigma2Groups.all[1,] = sigma2Groups
sigma20.all[1] = sigma20 
nu0.all[1] = nu0

### Pozenemo Gibbsov vzorčevalnik

set.seed(1)
for (s in 1 : n.iter) {
  
  # Vzorcimo muGroups
  for (j in 1 : m) {
    # SPREMENJENO: sigma2 v enačbi je zamenjana s sigma2Groups[j], ostalo je enako
    muGroups[j] = rnorm(1,
                        mean = (x.povpr[j] * n[j] / sigma2Groups[j] + mu / eta2) / (n[j] / sigma2Groups[j] + 1 / eta2),
                        sd = sqrt(1 / (n[j] / sigma2Groups[j] + 1 / eta2)))
  }
  
  # DODANO: Vzorcimo sigma2Groups namesto sigma2, delamo po skupinah... po formuli iz domace naloge
  for(j in 1 : m){
    sigma2Groups[j] = 1 / rgamma(1,
                                 (nu0 + n[j]) / 2,
                                 (nu0 * sigma20 + sum((x[x[, 1] == j, 2] - muGroups[j])^2)) / 2)
  }
  
  # DODANO: Vzorcimo sigma20... po formuli iz domace
  sigma20 = rgamma(1,
                   a + m * nu0 / 2,
                   b + nu0 * sum(1 / sigma2Groups) / 2)
  
  # Vzorcimo mu, ne spreminjamo nic
  mu <- rnorm(1,
              mean = (mean(muGroups) * m / eta2 + mu0 / tau20) / (m / eta2 + 1 / tau20),
              sd = sqrt(1 / (m / eta2 + 1 / tau20)))
  
  # Vzorcimo eta2, ne spreminjamo nic
  ss   <- kappa0 * eta20 + sum((muGroups - mu)^2)
  eta2 <- 1 / rgamma(1, (kappa0 + m) / 2, ss / 2)
  
  # DODANO: Vzorcimo nu0, koda kopirana iz navodil za domaco nalogo
  k <- 1:k.max
  logp.nu0 <- m * (0.5 * k * log(k*sigma20/2) - lgamma(k/2)) +
    (k/2-1) * sum(log(1/sigma2Groups)) +
    - k * (alpha + 0.5 * sigma20 * sum(1/sigma2Groups))
  nu0 <- sample(k, 1, prob = exp(logp.nu0 - max(logp.nu0)))
  
  # Shranimo nove parametre
  muGroups.all[s,] = muGroups 
  mu.all[s] = mu
  eta2.all[s] = eta2
  
  # DODANO:
  sigma2Groups.all[s,] = sigma2Groups
  sigma20.all[s] = sigma20 
  nu0.all[s] = nu0
}

########################################################################################################################
########################################################################################################################

# 2. NALOGA

########################################

# Trace plots:

# Hiperparametri
par(mfrow=c(2,2))
plot(mu.all, type="l", main="mu") 
plot(eta2.all, type="l", main="eta2") 
plot(sigma20.all, type="l", main="sigma20")
plot(nu0.all, type="l", main="nu0") 


# Podobno še za samo prvih 500 iteracij
par(mfrow=c(2,2))
plot(mu.all[1:500], type="l", main="mu") 
plot(eta2.all[1:500], type="l", main="eta2") 
plot(sigma20.all[1:500], type="l", main="sigma20")
plot(nu0.all[1:500], type="l", main="nu0") 


# Izberemo j = 1, 10, 50, 80
# MUGROUPS
par(mfrow=c(2,2))
plot(muGroups.all[,1], type="l", main="muGroups, j = 1") 
plot(muGroups.all[,10], type="l", main="muGroups, j = 10") 
plot(muGroups.all[,50], type="l", main="muGroups, j = 50") 
plot(muGroups.all[,80], type="l", main="muGroups, j = 80")

# SIGMA2GROUPS
par(mfrow=c(2,2))
plot(sigma2Groups.all[,1], type="l", main="sigma2Groups, j = 1") 
plot(sigma2Groups.all[,10], type="l", main="sigma2Groups, j = 10") 
plot(sigma2Groups.all[,50], type="l", main="sigma2Groups, j = 50") 
plot(sigma2Groups.all[,80], type="l", main="sigma2Groups, j = 80")


# Podobno še za samo prvih 500 iteracij:
par(mfrow=c(2,2))
plot(muGroups.all[1:500, 1], type="l", main="muGroups,j = 1") 
plot(muGroups.all[1:500, 10], type="l", main="muGroups, j = 10") 
plot(muGroups.all[1:500, 50], type="l", main="muGroups, j = 50") 
plot(muGroups.all[1:500, 80], type="l", main="muGroups, j = 80")

par(mfrow=c(2,2))
plot(sigma2Groups.all[1:500, 1], type="l", main="sigma2Groups, j = 1") 
plot(sigma2Groups.all[1:500, 10], type="l", main="sigma2Groups, j = 10") 
plot(sigma2Groups.all[1:500, 50], type="l", main="sigma2Groups, j = 50") 
plot(sigma2Groups.all[1:500, 80], type="l", main="sigma2Groups, j = 80")


########################################

# Porazdelitve podvzorcev

library(gridExtra)


# Hiperparametri

mu.all2 = data.frame(sample = mu.all, podvzorec = factor(sort(rep(1:10,500)))) 
p1 = ggplot(mu.all2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "mu")

eta2.all2 = data.frame(sample = eta2.all, podvzorec = factor(sort(rep(1:10,500)))) 
p2 = ggplot(eta2.all2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "eta2")

sigma20.all2 = data.frame(sample = sigma20.all, podvzorec = factor(sort(rep(1:10,500)))) 
p3 = ggplot(sigma20.all2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "sigma20")

nu0.all2 = data.frame(sample = nu0.all, podvzorec = factor(sort(rep(1:10,500)))) 
p4 = ggplot(nu0.all2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "nu0")

grid.arrange(p1, p2, p3, p4, ncol = 2)


# Ostali, izberemo iste j kot prej

### MUGROUPS

mu1_1= data.frame(sample = muGroups.all[,1], podvzorec = factor(sort(rep(1:10,500)))) 
p5 = ggplot(mu1_1, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "muGroups, j = 1")

mu1_2 = data.frame(sample = muGroups.all[,10], podvzorec = factor(sort(rep(1:10,500)))) 
p6 = ggplot(mu1_2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "muGroups, j = 10")

mu1_3 = data.frame(sample = muGroups.all[,50], podvzorec = factor(sort(rep(1:10,500)))) 
p7 = ggplot(mu1_3, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "muGroups, j = 50")

mu1_4 = data.frame(sample = muGroups.all[,80], podvzorec = factor(sort(rep(1:10,500)))) 
p8 = ggplot(mu1_4, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "muGroups, j = 80")

grid.arrange(p5, p6, p7, p8, ncol = 2)

### SIGMA2GROUPS

sigma1_1 = data.frame(sample = sigma2Groups.all[,1], podvzorec = factor(sort(rep(1:10,500)))) 
p9 = ggplot(sigma1_1, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "sigma2Groups, j = 1")

sigma1_2 = data.frame(sample = sigma2Groups.all[,10], podvzorec = factor(sort(rep(1:10,500)))) 
p10 = ggplot(sigma1_2, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "sigma2Groups, j = 10")

sigma1_3 = data.frame(sample = sigma2Groups.all[,50], podvzorec = factor(sort(rep(1:10,500)))) 
p11 = ggplot(sigma1_3, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "sigma2Groups, j = 50")

sigma1_4 = data.frame(sample = sigma2Groups.all[,80], podvzorec = factor(sort(rep(1:10,500)))) 
p12 = ggplot(sigma1_4, aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "sigma2Groups, j = 80")

grid.arrange(p9, p10, p11, p12, ncol = 2)

########################################


# AVTOKORELACIJE:

###

# Hiperparametri
par(mfrow=c(2,2)) 
acf(mu.all, main = "mu")
acf(eta2.all, main = "eta2") 
acf(sigma20.all, main = "sigma20")
acf(nu0.all, main = "nu0")

# MUGROUPS
par(nfrow = c(2,2))
acf(muGroups.all[,1], main = "muGroups, j = 1") 
acf(muGroups.all[,10], main = "muGroups, j = 10") 
acf(muGroups.all[,50], main = "muGroups, j = 50") 
acf(muGroups.all[,80], main = "muGroups, j = 80") 

#SIGMA2GROUPS
par(nfrow = c(2,2))
acf(sigma2Groups.all[,1], main = "sigma2Groups, j = 1") 
acf(sigma2Groups.all[,10],main = "sigma2Groups, j = 10") 
acf(sigma2Groups.all[,50], main = "sigma2Groups, j = 50") 
acf(sigma2Groups.all[,80],main = "sigma2Groups, j = 80") 

###

# AVTOKORELACIJE se za samo prvih 100 iteracij

# Hiperparametri
par(mfrow = c(2,2)) 
acf(mu.all[1:100], main = "mu")
acf(eta2.all[1:100], main = "eta2") 
acf(sigma20.all[1:100], main = "sigma20")
acf(nu0.all[1:100], main = "nu0")

# MUGROUPS
par(nfrow = c(2,2))
acf(muGroups.all[1:100,1], main = "muGroups, j = 1") 
acf(muGroups.all[1:100,10], main = "muGroups, j = 10") 
acf(muGroups.all[1:100,50], main = "muGroups, j = 50") 
acf(muGroups.all[1:100,80], main = "muGroups, j = 80") 

#SIGMA2GROUPS
par(nfrow = c(2,2))
acf(sigma2Groups.all[1:100,1], main = "sigma2Groups, j = 1") 
acf(sigma2Groups.all[1:100,10], main = "sigma2Groups, j = 10") 
acf(sigma2Groups.all[1:100,50], main = "sigma2Groups, j = 50") 
acf(sigma2Groups.all[1:100,80], main = "sigma2Groups, j = 80") 

###

# AVTOKORELACIJE s thinningom, vsakega drugega izbrisemo

# Hiperparametri
par(mfrow = c(2,2)) 
acf(mu.all[seq(1, length(mu.all), by=2)], main = "mu")
acf(eta2.all[seq(1, length(mu.all), by=2)], main = "eta2") 
acf(sigma20.all[seq(1, length(mu.all), by=2)], main = "sigma20")
acf(nu0.all[seq(1, length(mu.all), by=2)], main = "nu0")

# MUGROUPS
par(nfrow = c(2,2))
acf(muGroups.all[seq(1, length(mu.all), by=2),1], main = "muGroups, j = 1") 
acf(muGroups.all[seq(1, length(mu.all), by=2),10], main = "muGroups, j = 10") 
acf(muGroups.all[seq(1, length(mu.all), by=2),50], main = "muGroups, j = 50") 
acf(muGroups.all[seq(1, length(mu.all), by=2),80], main = "muGroups, j = 80") 

#SIGMA2GROUPS
par(nfrow = c(2,2))
acf(sigma2Groups.all[seq(1, length(mu.all), by=2),1], main = "sigma2Groups, j = 1") 
acf(sigma2Groups.all[seq(1, length(mu.all), by=2),10], main = "sigma2Groups, j = 10") 
acf(sigma2Groups.all[seq(1, length(mu.all), by=2),50], main = "sigma2Groups, j = 50") 
acf(sigma2Groups.all[seq(1, length(mu.all), by=2),80], main = "sigma2Groups, j = 80") 


# Effective sample size:
library(coda) 

# Hiperparametri
effectiveSize(mu.all) 
effectiveSize(eta2.all)
effectiveSize(sigma20.all)
effectiveSize(nu0.all)

# MUGROUPS
effectiveSize(muGroups.all[,1])
effectiveSize(muGroups.all[,10])
effectiveSize(muGroups.all[,50])
effectiveSize(muGroups.all[,80])

#SIGMA2GROUPS
effectiveSize(sigma2Groups.all[,1])
effectiveSize(sigma2Groups.all[,10])
effectiveSize(sigma2Groups.all[,50])
effectiveSize(sigma2Groups.all[,80])

########################################

# 3. NALOGA

###

# Hiperparametri

par(mfrow=c(2,2))

plot(density(mu.all), type="l", main="mu")
abline(v = quantile(mu.all, prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(eta2.all), type="l", main="eta2")
abline(v = quantile(eta2.all, prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(sigma20.all), type="l", main="sigma20")
abline(v = quantile(sigma20.all, prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

hist(nu0.all, main = "nu0")
abline(v = quantile(nu0.all, prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red")

###

# MUGROUPS

par(mfrow=c(2,2))

plot(density(muGroups.all[,1]), type="l", main="muGroups, j = 1")
abline(v = quantile(muGroups.all[,1], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(muGroups.all[,10]), type="l", main="muGroups, j = 10")
abline(v = quantile(muGroups.all[,10], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(muGroups.all[,50]), type="l", main="muGroups, j = 50")
abline(v = quantile(muGroups.all[,50], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(muGroups.all[,80]), type="l", main="muGroups, j = 80")
abline(v = quantile(muGroups.all[,80], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red")

###

# SIGMA2GROUPS

par(mfrow=c(2,2))

plot(density(sigma2Groups.all[,1]), type="l", main="sigma2Groups, j = 1")
abline(v = quantile(sigma2Groups.all[,1], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(sigma2Groups.all[,10]), type="l", main="sigma2Groups, j = 10")
abline(v = quantile(sigma2Groups.all[,10], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(sigma2Groups.all[,50]), type="l", main="sigma2Groups, j = 50")
abline(v = quantile(sigma2Groups.all[,50], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red") 

plot(density(sigma2Groups.all[,80]), type="l", main="sigma2Groups, j = 80")
abline(v = quantile(sigma2Groups.all[,80], prob=c(0.025, 0.5, 0.975)), lty = 2, col = "red")

###

# Intervali:

# Hiperparametri
quantile(mu.all, prob=c(0.025, 0.5, 0.975))
quantile(eta2.all, prob=c(0.025, 0.5, 0.975))
quantile(sigma20.all, prob=c(0.025, 0.5, 0.975))
quantile(nu0.all, prob=c(0.025, 0.5, 0.975))

# MUGROUPS
quantile(muGroups.all[,1], prob=c(0.025, 0.5, 0.975))
quantile(muGroups.all[,10], prob=c(0.025, 0.5, 0.975))
quantile(muGroups.all[,50], prob=c(0.025, 0.5, 0.975))
quantile(muGroups.all[,80], prob=c(0.025, 0.5, 0.975))

# SIGMA2GROUPS
quantile(sigma2Groups.all[,1], prob=c(0.025, 0.5, 0.975))
quantile(sigma2Groups.all[,10], prob=c(0.025, 0.5, 0.975))
quantile(sigma2Groups.all[,50], prob=c(0.025, 0.5, 0.975))
quantile(sigma2Groups.all[,80], prob=c(0.025, 0.5, 0.975))

########################################

# 4. NALOGA

# MUGROUPS povprecje
pod.sole$EmuGroups = colMeans(muGroups.all)

par(mfrow=c(1,2))

plot(x = pod.sole$povprecje, 
     y = pod.sole$EmuGroups,
     xlab = "vzorcno povprecje", 
     ylab = expression(E(mu_j))) 
abline(a = 0, b = 1)

plot(x = pod.sole$n, 
     y = pod.sole$povprecje - pod.sole$EmuGroups,
     xlab = "velikost vzorca sole",
     ylab = expression(paste("vzorcno povprecje - "," ",E(mu_j), sep="")))
abline(h = 0)

###

# SIGMA2GROUPS varianca
pod.sole$Esigma2Groups = colMeans(sigma2Groups.all)

par(mfrow=c(1,2))

plot(x = pod.sole$varianca, 
     y = pod.sole$Esigma2Groups,
     xlab = "varianca", 
     ylab = expression(E(sigma2_j))) 
abline(a = 0, b = 1)

plot(x = pod.sole$n, 
     y = pod.sole$varianca - pod.sole$Esigma2Groups,
     xlab = "velikost vzorca sole",
     ylab = expression(paste("varianca - "," ",E(sigma2_j), sep="")))
abline(h = 0)


########################################




