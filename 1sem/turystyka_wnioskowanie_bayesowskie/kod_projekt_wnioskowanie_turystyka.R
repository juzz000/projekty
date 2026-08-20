install.packages("openxlsx")
install.packages("MASS")
install.packages("coda")
install.packages("MCMCpack")
install.packages("writexl")

library(openxlsx)
library(MASS)
library(MCMCpack)
library(coda)
library(writexl)

dane <- read.xlsx("dane_projekt_wnioskowanie.xlsx", sheet=1, startRow = 1, colNames = TRUE, rowNames = FALSE,
                  detectDates = FALSE)

colnames(dane)
typeof(dane)
ncol(dane)
nrow(dane)

y <- as.matrix(dane[,1]) 
X <- as.matrix(dane[, -1])

T_obs <- nrow(y)
k <- ncol(X)

N_spal <- 20000 
N_post <- 100000
T_total <- N_spal + N_post

a <- matrix(0, nrow = k, ncol = 1) 
C <- 0.0001 * diag(k)

n0 <- 0.001 
s0 <- 0.001

macierz_losowan <- matrix(0, nrow = T_total, ncol = k + 1)
colnames(macierz_losowan) <- c(colnames(X), "tau")

beta_0 <- matrix(0, nrow = k, ncol = 1)
tau_0 <- 1

beta_biezace <- beta_0
tau_biezace <- 1

beta_daszek <- solve(t(X)%*%X) %*% t(X) %*% y

n_falka <- n0 + T_obs
wektor_reszt <- y - X %*% beta_biezace
S_beta <- t(wektor_reszt) %*% wektor_reszt 
s_falka <- s0 + S_beta

C_falka <- C + tau_biezace * t(X) %*% X
a_falka <- solve(C_falka) %*% (C %*% a + tau_biezace * t(X) %*% X %*% beta_daszek)

for(q in 1:T_total){
  
  tau_biezace <- rgamma(1, shape = n_falka/2, rate = s_falka/2)

  C_falka <- C + tau_biezace * t(X) %*% X
  a_falka <- solve(C_falka) %*% (C %*% a + tau_biezace * t(X) %*% X %*% beta_daszek)
  
  beta_biezace <- MASS::mvrnorm(1, mu = a_falka, Sigma = solve(C_falka))
  beta_biezace <- matrix(beta_biezace, k, 1)
  
  wektor_reszt <- y - X %*% beta_biezace
  S_beta <- t(wektor_reszt) %*% wektor_reszt 
  s_falka <- s0 + S_beta
  
  macierz_losowan[q,] <- c(beta_biezace, tau_biezace)
  
  if(q %% 5000 == 0) cat("Obieg:", q, "\n")
}

posterior <- as.mcmc(macierz_losowan[(N_spal + 1):T_total, ])

jpeg("Wykresy_posteriori_priori.jpg", width = 1600, height = 1200, res = 150)

par(mfrow=c(4,5))

for(i in 1:(k+1)){
    histogram <- hist(posterior[,i], probability = TRUE, main=paste(colnames(posterior)[i]), 
                      col="lightgray", border="white")
    if(colnames(posterior)[i] != "tau"){
    lines(histogram$mids, dnorm(histogram$mids, mean=a[i], sd=(solve(C[i,i])^0.5)), col="red", lwd=2)
    
  }
  
  if(colnames(posterior)[i] == "tau"){
    
    lines(histogram$mids, dgamma(histogram$mids, shape=n0/2, rate=s0/2),
          col="blue", lwd=2)
  }
}

dev.off()

E_priori = c(rep(0, k), (n0/2)/(s0/2))
Me_priori = c(rep(0, k), qgamma(0.5, shape=n0/2, rate=s0/2))
Mo_priori = c(rep(0, k), ifelse((n0/2) > 1, ((n0/2) - 1) / (s0/2), 0))
Sd_priori = c(sqrt(diag(solve(C))), sqrt((n0/2))/(s0/2))

E_posteriori <- apply(posterior, 2, mean)
Me_posteriori <- apply(posterior, 2, quantile, probs=0.5)
Mo_posteriori <- numeric(ncol(posterior))
for(i in 1:ncol(posterior)){
  d <- density(posterior[,i])
  Mo_posteriori[i] <- d$x[which.max(d$y)]
}
Sd_posteriori <- apply(posterior, 2, sd)
CI_posteriori <- t(apply(posterior, 2, quantile, probs=c(0.025,0.975)))
HPD <- HPDinterval(posterior)

characteristics <- cbind(E_priori, E_posteriori, Me_priori, Me_posteriori, Mo_priori, Mo_posteriori, Sd_priori, Sd_posteriori, CI_posteriori, HPD)      
print(characteristics)

write_xlsx(as.data.frame(characteristics), "Charakterystyki.xlsx")

posterior_matrix <- as.matrix(posterior)
names <- colnames(posterior_matrix)

jpeg("Analiza_zbieznosci.jpg", width = 5000, height = 9000, res = 150)
par(mfrow=c((k+1),3))

for(i in 1:(k+1)){
  
  plot(posterior_matrix[,i], type="l",
       main=paste("Trace plot:",names[i]),
       xlab="Iteracja", ylab="Wartość")
  
  acf(posterior_matrix[,i],
      main=paste("ACF:",names[i]))
  
  param <- posterior_matrix[,i]
  erg_mean <- cumsum(param)/(1:length(param))
  
  plot(erg_mean, type="l",
       main=paste("Średnie ergodyczne:",names[i]),
       xlab="Iteracja", ylab="Średnia")
  
  abline(h=E_posteriori[i], col="red", lwd=2)
}

dev.off()

jpeg(filename = "Wykres_CUSUM.jpg", width = 1200, height = 800, res = 150)

plot(NULL, xlim=c(1,nrow(posterior_matrix)), ylim=c(-2,2),
     xlab="Iteracja", ylab="CUSUM",
     main="CUSUM") 

kolory <- rainbow(k+1)

for(i in 1:(k+1)){
  
  param <- posterior_matrix[, i]
  N <- length(param)
  
  m_N <- mean(param)
  s_N <- sd(param) 
  m_t <- cumsum(param)/(1:N)
  
  cs_t <- (m_t-m_N)/s_N
  
  lines(cs_t, col = kolory[i], lwd = 1)
}

abline(h = 0, col = "red", lwd = 2, lty = 2)

legend("topright", legend = colnames(posterior_matrix), 
       col = kolory, lty = 1, cex = 0.6, ncol = 2)
dev.off()



posterior_funkcja <- MCMCregress(Turystyka ~ trend +  M2 + M3 + M4 + M5 + M6 + M7 + M8 + M9 + M10 + M11 + M12 
                                  + CPI + Bezrobocie + Płaca + Temperatura, 
                                  data = dane, burnin = N_spal, mcmc = N_post,
                                  thin = 1, verbose = 200, beta.start = 0, 
                                  b0 =a, B0 = C, c0 = n0, d0 = s0,
                                  marginal.likelihood = c("none"))

posterior_funkcja[,k+1] <- 1/posterior_funkcja[,k+1]
colnames(posterior_funkcja)[k+1] <- c("tau")
E_funkcja <- apply(posterior_funkcja, 2, mean)
SD_funkcja <- apply(posterior_funkcja, 2, sd) 

options("scipen" = 100, "digits" = 3)
comparison <- cbind(E_posteriori, E_funkcja, Sd_posteriori, SD_funkcja)
print(comparison)      
write_xlsx(as.data.frame(comparison), "Porownanie_z_funkcja.xlsx")

form_full <- Turystyka ~ trend + M2+M3+M4+M5+M6+M7+M8+M9+M10+M11+M12 +
  CPI + Bezrobocie + Płaca + Temperatura

form_no_cpi <- update(form_full, . ~ . - CPI)
form_no_bez <- update(form_full, . ~ . - Bezrobocie)
form_no_pla <- update(form_full, . ~ . - Płaca)
form_no_tem <- update(form_full, . ~ . - Temperatura)

forms <- list(full=form_full,
              no_cpi=form_no_cpi,
              no_bez=form_no_bez,
              no_pla=form_no_pla,
              no_tem=form_no_tem)

logml <- numeric(length(forms))


for(i in 1:length(forms)){
  
  form <- forms[[i]]
  
  z <- model.matrix(form, dane)
  k <- ncol(z)
  
  a <- matrix(0, nrow = k, ncol = 1)
  C <- 0.0001 * diag(k) 
  beta_0 <- matrix(0, nrow = k, ncol = 1)
  
  fit <- MCMCregress(form,
                     data = dane,
                     burnin = N_spal,
                     mcmc = N_post,
                     thin = 1,
                     verbose = 0,
                     beta.start = beta_0,
                     b0 = a,
                     B0 = C,
                     c0 = n0,
                     d0 = s0,
                     marginal.likelihood = "Chib95")
  
  logml[i] <- attr(fit, "logmarglike")
  
  cat("Model", i, "logML =", logml[i], "\n")
}

options("scipen" = 100, "digits" = 3)
exp(logml)
write_xlsx(as.data.frame(exp(logml)), "Porownanie_modeli.xlsx")
