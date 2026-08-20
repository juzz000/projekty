library(quantmod)
library(zoo)
library(rugarch)
library(rmgarch)
library(psych)
library(TSA)
library(MTS)

MXX  <- read.csv("BMVIPC_1991_2001.csv",  stringsAsFactors = FALSE)
MERV <- read.csv("Merval_1991_2001.csv",   stringsAsFactors = FALSE)
IBOV <- read.csv("Bovespa_1991_2001.csv",  stringsAsFactors = FALSE)

MXX  <- MXX[, 1:2]
MERV <- MERV[, 1:2]
IBOV <- IBOV[, 1:2]

MXX$Date  <- as.Date(MXX$Date,  format = "%m/%d/%Y")
MERV$Date <- as.Date(MERV$Date, format = "%m/%d/%Y")
IBOV$Date <- as.Date(IBOV$Date, format = "%m/%d/%Y")

MXX[, 2]  <- as.numeric(gsub(",", "", MXX[, 2]))
MERV[, 2] <- as.numeric(gsub(",", "", MERV[, 2]))
IBOV[, 2] <- as.numeric(gsub(",", "", IBOV[, 2]))

MXX  <- MXX[order(MXX$Date), ]
MERV <- MERV[order(MERV$Date), ]
IBOV <- IBOV[order(IBOV$Date), ]

msci_data <- getSymbols("^892000-USD-STRD", src = "yahoo",
                        from = "1991-01-01", to = "2001-12-31",
                        auto.assign = FALSE)
ts_msci     <- zoo(coredata(msci_data[, 4]), order.by = as.Date(index(msci_data)))
colnames(ts_msci) <- "MSCI_LatAm"

dane_wstepne <- merge(MXX, MERV, by = "Date", all = TRUE)
dane_wstepne <- merge(dane_wstepne, IBOV, by = "Date", all = TRUE)
colnames(dane_wstepne) <- c("Date", "MXX", "MERV", "IBOV")

ts_csv <- zoo(dane_wstepne[, -1], order.by = dane_wstepne$Date)
dane   <- merge(ts_csv, ts_msci, all = TRUE)
dane   <- as.matrix(dane)
dane   <- dane[-c(1:627), ]
dane   <- na.omit(dane)

dim(dane)
sum(is.na(dane))

lnrdane <- diff(log(dane))
colnames(lnrdane) <- c("MXX", "MERV", "IBOV", "MSCI_LatAm")
daty  <- as.Date(rownames(lnrdane))
daty_dane  <- as.Date(rownames(dane))

par(mfrow = c(2, 1))
plot(daty_dane, dane[,1],    type = "l", main = colnames(dane)[1],    ylab = "Ceny")
plot(daty, lnrdane[,1], type = "l", main = colnames(lnrdane)[1], ylab = "Stopay zwrotu")
plot(daty_dane, dane[,2],    type = "l", main = colnames(dane)[2],    ylab = "Ceny")
plot(daty, lnrdane[,2], type = "l", main = colnames(lnrdane)[2], ylab = "Stopay zwrotu")
plot(daty_dane, dane[,3],    type = "l", main = colnames(dane)[3],    ylab = "Ceny")
plot(daty, lnrdane[,3], type = "l", main = colnames(lnrdane)[3], ylab = "Stopay zwrotu")

par(mfrow = c(3, 1))
TSA::acf(dane[, 1],       main = "MXX - poziom")  #brak stacjonarności
TSA::acf(lnrdane[, 1],    main = "MXX - stopy zwrotu")  #dosyć szybko zgiegają się do zera, choć pojawiają się póżne opóźnienia statystycznie istotne, szereg stacjonarny
TSA::acf(lnrdane[, 1]^2,  main = "MXX - kwadraty")  #zbiegają się do 0, choć w wolniejszym tempie niż logarytmy, efekt ARCH
par(mfrow = c(1, 1))
Box.test(lnrdane[, 1],    lag = 20, type = "Ljung-Box") #wskazuje na występowanie autokorelacji
Box.test(lnrdane[, 1]^2,  lag = 20, type = "Ljung-Box")
archTest(lnrdane[, 1])
McLeod.Li.test(y = lnrdane[, 1])
#wszystkie testy mówią o występowaniu efektu ARCH
options(scipen=999)
par(mfrow = c(3, 1))
TSA::acf(dane[, 2],       main = "MERV - poziom") #analogicznie
TSA::acf(lnrdane[, 2],    main = "MERV - stopy zwrotu") #analogicznie
TSA::acf(lnrdane[, 2]^2,  main = "MERV - kwadraty") #analogicznie
par(mfrow = c(1, 1))
Box.test(lnrdane[, 2],    lag = 20, type = "Ljung-Box") #wskazuje na występowanie autokorelacji
Box.test(lnrdane[, 2]^2,  lag = 20, type = "Ljung-Box")
archTest(lnrdane[, 2])
McLeod.Li.test(y = lnrdane[, 2])
#wszystkie testy mówią o występowaniu efektu ARCH

par(mfrow = c(3, 1))
TSA::acf(dane[, 3],       main = "IBOV - poziom") #analogicznie
TSA::acf(lnrdane[, 3],    main = "IBOV - stopy zwrotu") #analogicznie
TSA::acf(lnrdane[, 3]^2,  main = "IBOV - kwadraty") #analogicznie
par(mfrow = c(1, 1))
Box.test(lnrdane[, 3],    lag = 20, type = "Ljung-Box") #wskazuje na występowanie autokorelacji
Box.test(lnrdane[, 3]^2,  lag = 20, type = "Ljung-Box")
archTest(lnrdane[, 3])
McLeod.Li.test(y = lnrdane[, 3])
#wszystkie testy mówią o występowaniu efektu ARCH

msci <- lnrdane[, 4]

y_resztyMSCI <- matrix(0, length(msci), 3)
colnames(y_resztyMSCI) <- c("MXX", "MERV", "IBOV")

pom1 <- lm(lnrdane[, 1] ~ 1 + msci)
y_resztyMSCI[, 1] <- pom1$residuals
pom2 <- lm(lnrdane[, 2] ~ 1 + msci)
y_resztyMSCI[, 2] <- pom2$residuals
pom3 <- lm(lnrdane[, 3] ~ 1 + msci)
y_resztyMSCI[, 3] <- pom3$residuals

spec_sGarch_norm = ugarchspec(variance.model=list(model= "sGARCH", garchOrder= c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                              distribution.model="norm")
spec_gjrGARCH_norm = ugarchspec(variance.model=list(model= "gjrGARCH", garchOrder= c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                                distribution.model="norm")
spec_gjrGARCH_std = ugarchspec(variance.model=list(model= "gjrGARCH", garchOrder= c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                               distribution.model="std")
spec_gjrGARCH_sstd = ugarchspec(variance.model=list(model= "gjrGARCH", garchOrder= c(1, 1)), mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                                distribution.model="sstd")


#szereg [,1]
fit_sGarch_norm<-ugarchfit(data = lnrdane[, 1], spec = spec_sGarch_norm, solver = "hybrid")
fit_gjrGARCH_norm<-ugarchfit(data = lnrdane[, 1], spec = spec_gjrGARCH_norm, solver = "hybrid")
fit_gjrGARCH_std<-ugarchfit(data = lnrdane[, 1], spec = spec_gjrGARCH_std, solver = "hybrid")
fit_gjrGARCH_sstd<-ugarchfit(data = lnrdane[, 1], spec = spec_gjrGARCH_sstd, solver = "hybrid")
infocriteria(fit_sGarch_norm)
infocriteria(fit_gjrGARCH_norm)
infocriteria(fit_gjrGARCH_std)
infocriteria(fit_gjrGARCH_sstd)
#najlepszy model fit_gjrGARCH_std (-5.319307)

#szereg [,2]
fit_sGarch_norm<-ugarchfit(data = lnrdane[,2], spec = spec_sGarch_norm, solver = "hybrid")
fit_gjrGARCH_norm<-ugarchfit(data = lnrdane[,2], spec = spec_gjrGARCH_norm, solver = "hybrid")
fit_gjrGARCH_std<-ugarchfit(data = lnrdane[,2], spec = spec_gjrGARCH_std, solver = "hybrid")
fit_gjrGARCH_sstd<-ugarchfit(data = lnrdane[,2], spec = spec_gjrGARCH_sstd, solver = "hybrid")
infocriteria(fit_sGarch_norm)
infocriteria(fit_gjrGARCH_norm)
infocriteria(fit_gjrGARCH_std)
infocriteria(fit_gjrGARCH_sstd)
#najlepszy model fit_gjrGARCH_sstd (-4.938041)

#szereg [,3]
fit_sGarch_norm<-ugarchfit(data = lnrdane[,3], spec = spec_sGarch_norm, solver = "hybrid")
fit_gjrGARCH_norm<-ugarchfit(data = lnrdane[,3], spec = spec_gjrGARCH_norm, solver = "hybrid")
fit_gjrGARCH_std<-ugarchfit(data = lnrdane[,3], spec = spec_gjrGARCH_std, solver = "hybrid")
fit_gjrGARCH_sstd<-ugarchfit(data = lnrdane[,3], spec = spec_gjrGARCH_sstd, solver = "hybrid")
infocriteria(fit_sGarch_norm)
infocriteria(fit_gjrGARCH_norm)
infocriteria(fit_gjrGARCH_std)
infocriteria(fit_gjrGARCH_sstd)
#najlepszy model fit_gjrGARCH_sstd (-4.477138)

y <- lnrdane[, 1:3]

mspec  <- multispec(replicate(3, spec_gjrGARCH_sstd))
multf  <- multifit(mspec, data = y)
spec_dcc_y <- dccspec(uspec = mspec, VAR = FALSE, dccOrder = c(1, 1), distribution = "mvt")
fit_dcc_y  <- dccfit(spec_dcc_y, data = y, fit.control = list(eval.se = TRUE), fit = multf)

print(fit_dcc_y@mfit$matcoef)

print(DCCtest(y, garchOrder = c(1, 1), n.lags = 2, solver = "solnp", solver.control = list()))  #H0 odrzucone

par(mfrow = c(3, 1))
plot(daty, rcor(fit_dcc_y)[1, 2, ], type = "l", main = "MXX-MERV")
plot(daty, rcor(fit_dcc_y)[1, 3, ], type = "l", main = "MXX-IBOV")
plot(daty, rcor(fit_dcc_y)[2, 3, ], type = "l", main = "MERV-IBOV")
par(mfrow = c(1, 1))

multf_bez  <- multifit(mspec, data = y_resztyMSCI)
spec_dcc_bez_regional <- dccspec(uspec = mspec, VAR = FALSE, dccOrder = c(1, 1), distribution = "mvt")
fit_DCC_bez_regional  <- dccfit(spec_dcc_bez_regional, data = y_resztyMSCI, fit.control = list(eval.se = TRUE), fit = multf_bez)

print(fit_DCC_bez_regional@mfit$matcoef)

print(DCCtest(y_resztyMSCI, garchOrder = c(1, 1), n.lags = 2, solver = "solnp", solver.control = list())) #H0 odrzucone

par(mfrow = c(3, 1))
plot(daty, rcor(fit_DCC_bez_regional)[1, 2, ], type = "l", main = "MXX-MERV")
plot(daty, rcor(fit_DCC_bez_regional)[1, 3, ], type = "l", main = "MXX-IBOV")
plot(daty, rcor(fit_DCC_bez_regional)[2, 3, ], type = "l", main = "MERV-IBOV")
