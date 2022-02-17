# Bayesova statistika, seminarska naloga

library(Rlab)
library(BayesX)
library(R2BayesX)

# stevilo simulcij 
N = 1000

# velikost vzorca
n = 100

# porazdelitve suma (napake)
napaka.1 = rnorm(n, mean = 0, sd = 2) # z majhno varianco
napaka.2 = rnorm(n, mean = 0, sd = 100) # z veliko varianco
napaka.3 = rchisq(n, df = 1)

# pojasnjevalne spremenljivke, simulirane neodvisno ena od druge
X.1 = rnorm(n, mean = 10, sd = 2)
X.2 = rbern(n, prob = 0.8)
X.3 = rnorm(n, mean = 15, sd = 5)
X.4 = rbinom(n, size = 10, prob = 0.7)
X.5 = rbern(n, prob = 0.2)

# koeficienti
intercept = 25 # zac vrednost 
beta.1 = 0 # brez vpliva
beta.2 = 5
beta.3 = 50 # mocan vpliv
beta.4 = 4 # srednje velik vpliv
beta.5 = 10 # sibak vpliv
  
# Generiramo spremenljivko Y:
Y = intercept + beta.1 * X.1 + beta.2 * X.2 + beta.3 * X.3 + beta.4 * X.4 + beta.5 * X.5

# Scenariji, razlicni sumi:
Y.1 = Y + napaka.1
Y.2 = Y + napaka.2
Y.3 = Y + napaka.3

# Shranimo v dataframe:
bayes.1 = cbind(Y.1, X.1, X.2, X.3, X.4, X.5)
bayes.2 = cbind(Y.2, X.1, X.2, X.3, X.4, X.5)
bayes.3 = cbind(Y.3, X.1, X.2, X.3, X.4, X.5)

df.1 = data.frame(bayes.1)
df.2 = data.frame(bayes.2)
df.3 = data.frame(bayes.3)

### KONVERGENCA ###
# Preverimo, ali so privzeti parametri funkcije primerni (stevilo iteracij, burn-in, thinning)
# MCMC verige

# model.bayes.1 = model z normalno porazdeljeno napako z majhno varianco
# model.bayes.2 = model z normalno porazdeljeno napako z veliko varianco
# model.bayes.3 = model z asimetrično porazdeljeno napako


model.bayes.1 = bayesx(Y.1 ~ X.1 + X.2 + X.3 + X.4 + X.5,
                       data = df.1,
                       method = "MCMC",
                       family = "gaussian")

model.bayes.2 = bayesx(Y.2 ~ X.1 + X.2 + X.3 + X.4 + X.5,
                       data = df.2,
                       method = "MCMC",
                       family = "gaussian")

model.bayes.3 = bayesx(Y.3 ~ X.1 + X.2 + X.3 + X.4 + X.5,
                       data = df.3,
                       method = "MCMC",
                       family = "gaussian")


#############################################################
# Model 1: z normalno porazdeljeno napako z majhno varianco #
#############################################################

library(gridExtra)
library(ggplot2)

b1.majhna = attr(model.bayes.1$fixed.effects, "sample")[,1]
b2.majhna = attr(model.bayes.1$fixed.effects, "sample")[,2]
b3.majhna = attr(model.bayes.1$fixed.effects, "sample")[,3]
b4.majhna = attr(model.bayes.1$fixed.effects, "sample")[,4]
b5.majhna = attr(model.bayes.1$fixed.effects, "sample")[,5]

# Celotna zaporedja:
plot(b1.majhna, type = "l", main = "beta.1, celotno zaporedje")
plot(b2.majhna, type = "l", main = "beta.2, celotno zaporedje")
plot(b3.majhna, type = "l", main = "beta.3, celotno zaporedje")
plot(b4.majhna, type = "l", main = "beta.4, celotno zaporedje")
plot(b5.majhna, type = "l", main = "beta.5, celotno zaporedje")

# Prvih 100 clenov zaporedja
plot(b1.majhna[1:100], type = "l", main = "beta.1, prvih 100 clenov zaporedja")
plot(b2.majhna[1:100], type = "l", main = "beta.2, prvih 100 clenov zaporedja")
plot(b3.majhna[1:100], type = "l", main = "beta.3, prvih 100 clenov zaporedja")
plot(b4.majhna[1:100], type = "l", main = "beta.4, prvih 100 clenov zaporedja")
plot(b5.majhna[1:100], type = "l", main = "beta.5, prvih 100 clenov zaporedja")

# Podvzorci
ggplot(data.frame(sample = data.frame(b1.majhna)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.1")
ggplot(data.frame(sample = data.frame(b2.majhna)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.2")
ggplot(data.frame(sample = data.frame(b3.majhna)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.3")
ggplot(data.frame(sample = data.frame(b4.majhna)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.4")
ggplot(data.frame(sample = data.frame(b5.majhna)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.5")

# Avtokorelacije
plot(acf(b1.majhna, plot = FALSE), main = "Avtokorelacija za beta.1")
plot(acf(b2.majhna, plot = FALSE), main = "Avtokorelacija za beta.2")
plot(acf(b3.majhna, plot = FALSE), main = "Avtokorelacija za beta.3")
plot(acf(b4.majhna, plot = FALSE), main = "Avtokorelacija za beta.4")
plot(acf(b5.majhna, plot = FALSE), main = "Avtokorelacija za beta.5")

#############################################################
# Model 2: z normalno porazdeljeno napako z veliko varianco #
#############################################################

b1.velika = attr(model.bayes.2$fixed.effects, "sample")[,1]
b2.velika = attr(model.bayes.2$fixed.effects, "sample")[,2]
b3.velika = attr(model.bayes.2$fixed.effects, "sample")[,3]
b4.velika = attr(model.bayes.2$fixed.effects, "sample")[,4]
b5.velika = attr(model.bayes.2$fixed.effects, "sample")[,5]

# Celotna zaporedja:
plot(b1.velika, type = "l", main = "beta.1, celotno zaporedje")
plot(b2.velika, type = "l", main = "beta.2, celotno zaporedje")
plot(b3.velika, type = "l", main = "beta.3, celotno zaporedje")
plot(b4.velika, type = "l", main = "beta.4, celotno zaporedje")
plot(b5.velika, type = "l", main = "beta.5, celotno zaporedje")

# Prvih 100 clenov zaporedja
plot(b1.velika[1:100], type = "l", main = "beta.1, prvih 100 clenov zaporedja")
plot(b2.velika[1:100], type = "l", main = "beta.2, prvih 100 clenov zaporedja")
plot(b3.velika[1:100], type = "l", main = "beta.3, prvih 100 clenov zaporedja")
plot(b4.velika[1:100], type = "l", main = "beta.4, prvih 100 clenov zaporedja")
plot(b5.velika[1:100], type = "l", main = "beta.5, prvih 100 clenov zaporedja")

# Podvzorci
ggplot(data.frame(sample = data.frame(b1.velika)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.1")
ggplot(data.frame(sample = data.frame(b2.velika)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.2")
ggplot(data.frame(sample = data.frame(b3.velika)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.3")
ggplot(data.frame(sample = data.frame(b4.velika)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.4")
ggplot(data.frame(sample = data.frame(b5.velika)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.5")


# Avtokorelacije
plot(acf(b1.velika, plot = FALSE), main = "Avtokorelacija za beta.1")
plot(acf(b2.velika, plot = FALSE), main = "Avtokorelacija za beta.2")
plot(acf(b3.velika, plot = FALSE), main = "Avtokorelacija za beta.3")
plot(acf(b4.velika, plot = FALSE), main = "Avtokorelacija za beta.4")
plot(acf(b5.velika, plot = FALSE), main = "Avtokorelacija za beta.5")

##############################################
# Model 3: z asimetrično porazdeljeno napako #
##############################################

b1.asim = attr(model.bayes.3$fixed.effects, "sample")[,1]
b2.asim = attr(model.bayes.3$fixed.effects, "sample")[,2]
b3.asim = attr(model.bayes.3$fixed.effects, "sample")[,3]
b4.asim = attr(model.bayes.3$fixed.effects, "sample")[,4]
b5.asim = attr(model.bayes.3$fixed.effects, "sample")[,5]

# Celotna zaporedja:
plot(b1.asim, type = "l", main = "beta.1, celotno zaporedje")
plot(b2.asim, type = "l", main = "beta.2, celotno zaporedje")
plot(b3.asim, type = "l", main = "beta.3, celotno zaporedje")
plot(b4.asim, type = "l", main = "beta.4, celotno zaporedje")
plot(b5.asim, type = "l", main = "beta.5, celotno zaporedje")

# Prvih 100 clenov zaporedja
plot(b1.asim[1:100], type = "l", main = "beta.1, prvih 100 clenov zaporedja")
plot(b2.asim[1:100], type = "l", main = "beta.2, prvih 100 clenov zaporedja")
plot(b3.asim[1:100], type = "l", main = "beta.3, prvih 100 clenov zaporedja")
plot(b4.asim[1:100], type = "l", main = "beta.4, prvih 100 clenov zaporedja")
plot(b5.asim[1:100], type = "l", main = "beta.5, prvih 100 clenov zaporedja")

# Podvzorci
ggplot(data.frame(sample = data.frame(b1.asim)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.1")
ggplot(data.frame(sample = data.frame(b2.asim)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.2")
ggplot(data.frame(sample = data.frame(b3.asim)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.3")
ggplot(data.frame(sample = data.frame(b4.asim)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.4")
ggplot(data.frame(sample = data.frame(b5.asim)[,1], podvzorec = factor(sort(rep(1:10,100)))), aes(x = podvzorec, y = sample)) + geom_boxplot() + labs(title = "Podvzorci za beta.5")


# Avtokorelacije
plot(acf(b1.asim, plot = FALSE), main = "Avtokorelacija za beta.1")
plot(acf(b2.asim, plot = FALSE), main = "Avtokorelacija za beta.2")
plot(acf(b3.asim, plot = FALSE), main = "Avtokorelacija za beta.3")
plot(acf(b4.asim, plot = FALSE), main = "Avtokorelacija za beta.4")
plot(acf(b5.asim, plot = FALSE), main = "Avtokorelacija za beta.5")

##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################
##########################################################################################################################################################################################

### PRIMERJAVA ###
##################

# Pomozne funkcije

simuliraj = function(napaka){
  
  # pojasnjevalne spremenljivke, simulirane neodvisno ena od druge
  X.1 = rnorm(n, mean = 10, sd = 2)
  X.2 = rbern(n, prob = 0.8)
  X.3 = rnorm(n, mean = 15, sd = 5)
  X.4 = rbinom(n, size = 10, prob = 0.7)
  X.5 = rbern(n, prob = 0.2)
  
  # koeficienti
  intercept = 25 # zac vrednost 
  beta.1 = 0 # brez vpliva
  beta.2 = 5
  beta.3 = 50 # mocan vpliv
  beta.4 = 4 # srednje velik vpliv
  beta.5 = 10 # sibak vpliv
  
  # Generiramo spremenljivko Y:
  Y = intercept + beta.1 * X.1 + beta.2 * X.2 + beta.3 * X.3 + beta.4 * X.4 + beta.5 * X.5
  Y = Y + napaka
  
  # Shranimo v dataframe:
  bayes = cbind(Y, X.1, X.2, X.3, X.4, X.5)
  df = data.frame(bayes)
  
  # Modela:
  
  # Bayesov:
  model.bayes = bayesx(Y ~ X.1 + X.2 + X.3 + X.4 + X.5,
                       data = df,
                       method = "MCMC",
                       family = "gaussian")
  
  # Frekventisticen
  model.frekventisticno = lm(Y ~ X.1 + X.2 + X.3 + X.4 + X.5, 
                             data = df)
  
  # Shranjevanje
  attr.bayes = attr(model.bayes$fixed.effects, "sample")
  
  rezultat = list("Koeficienti_frekvenstisticni" = coef(summary(model.frekventisticno))[,1],
                  "Standardna_napaka_frekvenstisticni" = coef(summary(model.frekventisticno))[,2],
                  "Koeficienti_Bayesov" = apply(attr.bayes, 2, mean),
                  "Standardna_napaka_Bayesov" = apply(attr.bayes, 2, sd))
                  
  return(rezultat)
}


##############################
# Model 1: z majhno varianco #
##############################

### 1

shrani_podatke.1 = t(replicate(N, simuliraj(napaka.1)))

koeficienti_frekv.1 = do.call(rbind, shrani_podatke.1[,1])
standardna_napaka_frekv.1 = do.call(rbind, shrani_podatke.1[,2])
koeficienti_bayes.1 = do.call(rbind, shrani_podatke.1[,3])
standardna_napaka_bayes.1 = do.call(rbind, shrani_podatke.1[,4])

beta_koef = c(intercept, beta.1, beta.2, beta.3, beta.4, beta.5)

tabela.1 = data.frame("Beta" = c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Ocena frekvenstisticni" = colMeans(koeficienti_frekv.1),
                      "Pristranskost frekventist" = beta_koef - colMeans(koeficienti_frekv.1),
                      "Ocena Bayesov" = colMeans(koeficienti_bayes.1),
                      "Pristranskost Bayesov" = beta_koef - colMeans(koeficienti_bayes.1))

#install.packages("huxtable")
library(huxtable)
library(magrittr)
library(matrixStats)

tabela.1_hux = hux(tabela.1)

number_format(tabela.1_hux) = 6
number_format(tabela.1_hux)[1] = 1
bold(tabela.1_hux)[1,] = TRUE
bottom_border(tabela.1_hux)[1,] = 0.5
right_border(tabela.1_hux)[,2] = 2
right_border(tabela.1_hux)[,4] = 2
align(tabela.1_hux) = "centre"

print_screen(tabela.1_hux)

### 2

tabela.2 = data.frame("Beta" =c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Standardna napaka frekventisticni" = colSds(koeficienti_frekv.1),
                      "Koren standardne napake frekventisticni" = sqrt( colSds(koeficienti_frekv.1)**2 + (colMeans(koeficienti_frekv.1) - beta_koef)**2 ),
                      "Standardna napaka Bayesov" = colSds(koeficienti_bayes.1),
                      "Koren standardne napake Bayesov" = sqrt(colSds(koeficienti_bayes.1)**2 + (colMeans(koeficienti_bayes.1)- beta_koef)**2)
)


tabela.2_hux = hux(tabela.2)

number_format(tabela.2_hux) = 6
bold(tabela.2_hux)[1,] = TRUE
bottom_border(tabela.2_hux)[1,] = 0.5
right_border(tabela.2_hux)[,2] = 0.3
right_border(tabela.2_hux)[,4] = 2
align(tabela.2_hux) = "centre"

print_screen(tabela.2_hux)

### 3

tabela.3 = data.frame("Beta" = c("Intercept", "Ena", "Dba", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Povprečje standardne napake frekventisticni" = colMeans(standardna_napaka_frekv.1),
                      "Povprečje standardne napake Bayesov" = colMeans(standardna_napaka_bayes.1)
                      )


tabela.3_hux = hux(tabela.3)

number_format(tabela.3_hux) = 6
bold(tabela.3_hux)[1,] = TRUE
bottom_border(tabela.3_hux)[1,] = 0.5
right_border(tabela.3_hux)[,2] = 0.3
right_border(tabela.3_hux)[,3] = 2
align(tabela.3_hux) = "centre"

print_screen(tabela.3_hux)


##############################
# Model 2: z veliko varianco #
##############################

# 1

shrani_podatke.2 = t(replicate(N, simuliraj(napaka.2)))

koeficienti_frekv.2 = do.call(rbind, shrani_podatke.2[,1])
standardna_napaka_frekv.2 = do.call(rbind, shrani_podatke.2[,2])
koeficienti_bayes.2 = do.call(rbind, shrani_podatke.2[,3])
standardna_napaka_bayes.2 = do.call(rbind, shrani_podatke.2[,4])

beta_koef = c(intercept, beta.1, beta.2, beta.3, beta.4, beta.5)

tabela.1 = data.frame("Beta" = c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Ocena frekvenstisticni" = colMeans(koeficienti_frekv.2),
                      "Pristranskost frekventist" = beta_koef - colMeans(koeficienti_frekv.2),
                      "Ocena Bayesov" = colMeans(koeficienti_bayes.2),
                      "Pristranskost Bayesov" = beta_koef - colMeans(koeficienti_bayes.2))

tabela.1_hux = hux(tabela.1)

number_format(tabela.1_hux) = 6
number_format(tabela.1_hux)[1] = 1
bold(tabela.1_hux)[1,] = TRUE
bottom_border(tabela.1_hux)[1,] = 0.5
right_border(tabela.1_hux)[,2] = 2
right_border(tabela.1_hux)[,4] = 2
align(tabela.1_hux) = "centre"

print_screen(tabela.1_hux)

# 2

tabela.2 = data.frame("Beta" =c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Standardna napaka frekventisticni" = colSds(koeficienti_frekv.2),
                      "Koren standardne napake frekventisticni" = sqrt(colSds(koeficienti_frekv.2)**2 + (colMeans(koeficienti_frekv.2) - beta_koef)**2),
                      "Standardna napaka Bayesov" = colSds(koeficienti_bayes.2),
                      "Koren standardne napake Bayesov" = sqrt(colSds(koeficienti_bayes.2)**2 + (colMeans(koeficienti_bayes.2)- beta_koef)**2)
)


tabela.2_hux = hux(tabela.2)

number_format(tabela.2_hux) = 6
bold(tabela.2_hux)[1,] = TRUE
bottom_border(tabela.2_hux)[1,] = 0.5
right_border(tabela.2_hux)[,1] = 0.3
right_border(tabela.2_hux)[,2] = 0.3
right_border(tabela.2_hux)[,4] = 2
align(tabela.2_hux) = "centre"

print_screen(tabela.2_hux)


# 3

tabela.3 = data.frame("Beta" = c("Intercept", "Ena", "Dba", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Povprečje standardne napake frekventisticni" = colMeans(standardna_napaka_frekv.2),
                      "Povprečje standardne napake Bayesov" = colMeans(standardna_napaka_bayes.2)
)


tabela.3_hux = hux(tabela.3)

number_format(tabela.3_hux) = 6
bold(tabela.3_hux)[1,] = TRUE
bottom_border(tabela.3_hux)[1,] = 0.5
right_border(tabela.3_hux)[,1] = 0.3
right_border(tabela.3_hux)[,2] = 0.3
align(tabela.3_hux) = "centre"

print_screen(tabela.3_hux)


####################################################
# Model 3: model z asimetrično porazdeljeno napako #
####################################################

# 1

shrani_podatke.3 = t(replicate(N, simuliraj(napaka.3)))

koeficienti_frekv.3 = do.call(rbind, shrani_podatke.3[,1])
standardna_napaka_frekv.3 = do.call(rbind, shrani_podatke.3[,2])
koeficienti_bayes.3 = do.call(rbind, shrani_podatke.3[,3])
standardna_napaka_bayes.3 = do.call(rbind, shrani_podatke.3[,4])

beta_koef = c(intercept, beta.1, beta.2, beta.3, beta.4, beta.5)

tabela.1 = data.frame("Beta" = c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Ocena frekvenstisticni" = colMeans(koeficienti_frekv.3),
                      "Pristranskost frekventist" = beta_koef - colMeans(koeficienti_frekv.3),
                      "Ocena Bayesov" = colMeans(koeficienti_bayes.3),
                      "Pristranskost Bayesov" = beta_koef - colMeans(koeficienti_bayes.3))

tabela.1_hux = hux(tabela.1)

number_format(tabela.1_hux) = 6
number_format(tabela.1_hux)[1] = 1
bold(tabela.1_hux)[1,] = TRUE
bottom_border(tabela.1_hux)[1,] = 0.5
right_border(tabela.1_hux)[,2] = 2
right_border(tabela.1_hux)[,4] = 2
align(tabela.1_hux) = "centre"

print_screen(tabela.1_hux)

# 2

tabela.2 = data.frame("Beta" =c("Intercept", "Ena", "Dva", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Standardna napaka frekventisticni" = colSds(koeficienti_frekv.3),
                      "Koren standardne napake frekventisticni" = sqrt(colSds(koeficienti_frekv.3)**2 + (colMeans(koeficienti_frekv.3) - beta_koef)**2),
                      "Standardna napaka Bayesov" = colSds(koeficienti_bayes.3),
                      "Koren standardne napake Bayesov" = sqrt(colSds(koeficienti_bayes.3)**2 + (colMeans(koeficienti_bayes.3)- beta_koef)**2)
)


tabela.2_hux = hux(tabela.2)

number_format(tabela.2_hux) = 6
bold(tabela.2_hux)[1,] = TRUE
bottom_border(tabela.2_hux)[1,] = 0.5
right_border(tabela.2_hux)[,1] = 0.3
right_border(tabela.2_hux)[,2] = 0.3
right_border(tabela.2_hux)[,4] = 2
align(tabela.2_hux) = "centre"

print_screen(tabela.2_hux)


# 3

tabela.3 = data.frame("Beta" = c("Intercept", "Ena", "Dba", "Tri", "Stiri", "Pet"),
                      "Tocne vrednosti beta" = beta_koef,
                      "Povprečje standardne napake frekventisticni" = colMeans(standardna_napaka_frekv.3),
                      "Povprečje standardne napake Bayesov" = colMeans(standardna_napaka_bayes.3)
)


tabela.3_hux = hux(tabela.3)

number_format(tabela.3_hux) = 6
bold(tabela.3_hux)[1,] = TRUE
bottom_border(tabela.3_hux)[1,] = 0.5
right_border(tabela.3_hux)[,1] = 0.3
right_border(tabela.3_hux)[,2] = 0.3
align(tabela.3_hux) = "centre"

print_screen(tabela.3_hux)
