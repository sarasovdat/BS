######################################################################################################################
############################### 1. naloga : implementacija Metropolis Hasting algoritma ############################## 
######################################################################################################################

# Dani podatki: 
x <- c(2.11, 9.75, 13.88, 11.3, 8.93, 15.66, 16.38, 4.54, 8.86, 11.94, 12.47, 11.11, 11.65, 14.53, 9.61, 7.38, 3.34, 9.06, 9.45, 5.98, 7.44, 8.5, 1.55, 11.45, 9.73)
mean(x)
plot(x, type = "l")

######################################################################################################################
# PRIVZET MODEL N(theta, sigma^2)
# Ocenjujemo theto
sigma <- 2

# APRIORNA = N(theta_0, tau_0^2)
theta_0 <- 6
tau_0 <- 3

# PREDLAGALNO JEDRO
sigma_predlagalno_jedro <- 0.1


######################################################################################################################

# UPORABIMO LOGARITEMSKO SKALO, DA NE PRIDE DO DELJENJA Z 0!!! 
# (zaokroževanje majhnih številk, ko npr vzamemo nesmiselno majhno vrednost... imamo počasno konvergenco in se za zelo spreminja)

# UPOŠTEVAMO SIMETRIČNOST NORMALNE PORAZDELITVE!


MetropolisHasting <- function(theta_0, koraki, sigma_jedro = sigma_predlagalno_jedro) {
  theta <- numeric(koraki)
  theta[1] = theta_0
  
  for (i in 2 : koraki) {
    predlagalno_jedro <- rnorm(1, mean = theta[i - 1], sd = sigma_jedro)
    
    # Trenutno
    theta_trenutna_apriorna <- dnorm(theta[i - 1], theta_0, tau_0, log = TRUE)
    theta_trenutna_verjetje <- sum(dnorm(x, theta[i - 1], sigma, log = TRUE))
    
    # Naslednje
    theta_naslednja_apriorna <- dnorm(predlagalno_jedro, theta_0, tau_0, log = TRUE)
    theta_naslednja_verjetje <- sum(dnorm(x, predlagalno_jedro, sigma, log = TRUE))
    
    # Kvocient
    kvocient <- exp(min(log(1), ((theta_naslednja_apriorna + theta_naslednja_verjetje) - (theta_trenutna_apriorna + theta_trenutna_verjetje))))
    ro <- runif(1, 0, 1)
    
    if (kvocient > ro) {
      # sprejmemo
      theta[i] <- predlagalno_jedro
    }
    else {
      # zavrnemo
      theta[i] <- theta[i - 1]
    }
  }
  theta
}
  
######################################################################################################################
############################### 2. naloga : Testiranje s smiselno začetno vrednostjo #################################
######################################################################################################################

# 2.1, S = 100k
vektor <- MetropolisHasting(5, 100000)
plot(vektor, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 5")

vektor <- MetropolisHasting(9, 100000)
plot(vektor, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 9")

######################################################################################################################

# 2.2
# S = 500
vektor1 <- MetropolisHasting(5, 500)
plot(vektor1, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 5")

vektor1 <- MetropolisHasting(9, 500)
plot(vektor1, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 9")

# S = 5000
vektor2 <- MetropolisHasting(5, 5000)
plot(vektor2, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 5")

vektor2 <- MetropolisHasting(9, 5000)
plot(vektor2, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = 9")

######################################################################################################################

# 2.3
vektor3 <- MetropolisHasting(5, 100000)
plot(vektor3[501:100000], type = "l", ylab = expression(theta), xlab = "Iteracija")

######################################################################################################################

# 2.4

n = length(x)

kvocient_sigma_n = sig^2 / n
theta_1 = (tau_0^2 / (kvocient_sigma_n + tau_0^2)) * mean(x) + (kvocient_sigma_n / (kvocient_sigma_n + tau_0^2)) * theta_0
tau_1 = sqrt((kvocient_sigma_n * tau_0^2) / (kvocient_sigma_n + tau_0^2))

theta <- seq(7, 11.5, 0.001)
aposteriorna_prava <- dnorm(theta, mean = theta_1, sd = tau_1)
aposteriorna_neprava <- density(vektor3[501:100000])

plot(theta, aposteriorna_prava, type = "l",
     xlab = expression(theta), ylab = "")
lines(aposteriorna_neprava, col = "red")
legend("topleft", legend = c("Aposteriorna porazdelitev (teoretična)", "Aposteriorna porazdelitev -- Metropolis Hasting"), col = c("black", "red"), lty = 1, bty = "n", cex = 0.7)

######################################################################################################################

# 2.5

# 95% interval zaupanja -- teoretično
qnorm(c(0.025, 0.975), mean = theta_1, sd = tau_1)

# 95% interval zaupanja -- M-H
# install.packages("HDInterval")
library(HDInterval)
hdi(vektor3[501:100000], credMass = 0.95)

######################################################################################################################
############################### 3. naloga : Testiranje z nesmiselno začetno vrednostjo ###############################
######################################################################################################################

# 3.1
vektor4 <- MetropolisHasting(-1, 100000)
plot(vektor4, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = -1")

######################################################################################################################

# 3.2
# S = 500
vektor5 <- MetropolisHasting(-1, 500)
plot(vektor5, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = -1")

# S = 5000
vektor6 <- MetropolisHasting(-1, 5000)
plot(vektor6, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Začetna vrednost = -1")

######################################################################################################################

# 3.3
vektor7 <- MetropolisHasting(-1, 100000)
plot(vektor7[1001:100000], type = "l", ylab = expression(theta), xlab = "Iteracija")


######################################################################################################################
################################ 4. naloga: Preiskovanje variance predlagalnega jedra ################################
######################################################################################################################

# 4.1

S = 5000
celotno = 100000

# Izberemo si začetno vrednost, kje začnemo risati spodnje grafe.
zv = 5

# a)
v1 <- MetropolisHasting(zv, S, sigma_jedro = 0.0001)
plot(v1, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.0001")

v2 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.0001)
plot(v2, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.0001")

# b)
v3 <- MetropolisHasting(zv, S, sigma_jedro = 0.0005)
plot(v3, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.0005")

v4 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.0005)
plot(v4, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.0005")

# c)
v5 <- MetropolisHasting(zv, S, sigma_jedro = 0.001)
plot(v5, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.001")

v6 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.001)
plot(v6, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.001")

# d)
v7 <- MetropolisHasting(zv, S, sigma_jedro = 0.005)
plot(v7, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.005")

v8 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.005)
plot(v8, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.005")

# e)
v9 <- MetropolisHasting(zv, S, sigma_jedro = 0.01)
plot(v9, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.01")

v10 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.01)
plot(v10, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.01")

# f)
v11<- MetropolisHasting(zv, S, sigma_jedro = 0.05)
plot(v11, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.05")

v12 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.05)
plot(v12, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.05")

# g)
v11<- MetropolisHasting(zv, S, sigma_jedro = 0.5)
plot(v11, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.5")

v12 <- MetropolisHasting(zv, celotno, sigma_jedro = 0.5)
plot(v12, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 0.5")

# h)
v13 <- MetropolisHasting(zv, S, sigma_jedro = 1)
plot(v13, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 1")

v14 <- MetropolisHasting(zv, celotno, sigma_jedro = 1)
plot(v14, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 1")

# i)
v15 <- MetropolisHasting(zv, S, sigma_jedro = 5)
plot(v15, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 5")

v16 <- MetropolisHasting(zv, celotno, sigma_jedro = 5)
plot(v16, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 5")

# j)
v17 <- MetropolisHasting(zv, S, sigma_jedro = 10)
plot(v17, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 10")

v18 <- MetropolisHasting(zv, celotno, sigma_jedro = 10)
plot(v18, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 10")

# k)
v19 <- MetropolisHasting(zv, S, sigma_jedro = 15)
plot(v19, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 15")

v20 <- MetropolisHasting(zv, celotno, sigma_jedro = 15)
plot(v20, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 15")

# l)
v21 <- MetropolisHasting(zv, S, sigma_jedro = 50)
plot(v21, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 50")

v22 <- MetropolisHasting(zv, celotno, sigma_jedro = 50)
plot(v22, type = "l", ylab = expression(theta), xlab = "Iteracija", main = "Predlagalno jedro = 50")


