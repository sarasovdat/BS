# Bayesova statistika: 1. domača naloga
# Oktober, 2021
# Sara Bizjak, 27202020


###############################################################################################################################
### 1. naloga:
# Pri vsaki razlicici apriorne beta apriorne porazdleitve (spreminjamo alpha in beta) narisemo graf: 
#     - apriorna porazdelitev,
#     - aposteriorna porazdelitev.

narisi_graf <- function(alpha, beta, n = 26, k = 6) {
  theta <- seq(0, 1, 0.001)
  
  # Apriorna:
  apriorna <- dbeta (theta, alpha, beta)
  
  # Aposterirna: 
  alpha.apost <- alpha + k
  beta.apost <- beta + n - k
  
  aposteriorna <- dbeta(theta, alpha.apost, beta.apost)
  
  # Graf:
  plot(theta, aposteriorna, type='l', ylab='', xlab=expression(theta), main = paste0( c("alpha = ", " beta = "), c(alpha, beta)))
  lines(theta, apriorna, col='red')
  legend("topright", legend = c("apriorna","aposteriorna"), col = c("red","black"), lty = 1, bty = "n", cex = 1.3)
}

# Spreminjamo alphe in bete za razlicne grafe:
narisi_graf(1000, 500)
  
###############################################################################################################################
### 2. naloga:
# Pri vsaki razlicici apriorne porazdelitve narisemo graf:
#     - apriorna porazdelitev,
#     - aposteriorna porazdelitev.
# Note: alpha in beta sta v vseh primerih izbrani tako, da je pricakovana vrednost apriorne porazdelitve enaka 1/4.

# Spreminjamo alphe in bete za vsako razlicico (veljati mora beta = 3 * alpha)
alpha <- 0.3
beta  <- 0.9

narisi_graf(alpha, beta)

# Izracunajmo se oceno pricakovane vrednosti za izbrana alpha in beta:
ocena_pricakovane_vrednosti <- function(alpha, beta, n = 26, k = 6) {
  ocena <- (k + alpha) / (n + alpha + beta)
  ocena
}

ocena_pricakovane_vrednosti(alpha, beta)

###############################################################################################################################
### 3. naloga:
# Imamo nov vzorec studentov velikosti 30, izmed katerih jih je 21 pravilno odgovorila na zastavljeno vprasanje. Porazdelitev je Beta(1, 1). 
# Aposteriorna porazdelitev je:

izracunaj_aposteriorno <- function(alpha = 1, beta = 1, n = 30, k = 21) {
  theta <- seq(0, 1, 0.001)
  
  alpha.apost <- alpha + k
  beta.apost <- beta + n - k
  
  aposteriorna <- dbeta(theta, alpha.apost, beta.apost)
  aposteriorna
}

narisi_aposteriorno <- function(alpha = 1, beta = 1, n = 30, k = 21) {
  theta <- seq(0, 1, 0.001)
  
  alpha.apost <- alpha + k
  beta.apost <- beta + n - k
  
  aposteriorna <- dbeta(theta, alpha.apost, beta.apost)
  
  plot(theta, aposteriorna, type='l', ylab='', xlab=expression(theta), main = c("Graf aposteriorne porazdelitve."))
}

narisi_aposteriorno()
izracunaj_aposteriorno()

###############################################################################################################################
### 4. naloga:
# Primerjava aposteriornih porazdelitev
#       - Z1 : Aposteriorna porazdelitev za 
#         Beta(1 + 21, 1 + 30 - 21) = Beta(22, 10)
#       - Z2 : Aposteriorna poradelitev za 
#         Beta(7, 21)

Z1 <- rbeta(10000, 22, 10)
Z2 <- rbeta(10000, 7, 21)

# Primerjamo vektorja po komponentah. Za vsako komponento nas zanima, če je element iz Z2 < elementa iz Z1

stej <- 0
for (i in 1 : 10000) {
  if (Z2[i] < Z1[i]) {
    stej <- stej + 1
  }
}
ocena <- stej / 10000
ocena

# 95% interval zaupanja
# Najprej za vsako porazd posebej
quantile(Z1, probs = c(0.025, 0.975))
quantile(Z2, probs = c(0.025, 0.975))

# Skupno -- na podlagi simulacije
# Generiramo vektor razlik Z1 in Z2 (simulacija zgornje for zanke, lahko bi gledali, kdaj razlika pozitivna) in na tem vektorju izracunajmo interval zaupanja
vektor_razlik <- Z1 - Z2
quantile(vektor_razlik, probs = c(0.025, 0.975))


