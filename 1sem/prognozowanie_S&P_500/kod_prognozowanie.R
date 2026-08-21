library(tsibble)
library(dplyr)
library(ggplot2)
library(forecast)
library(fable)
library(fabletools)
library(tseries)
library(fBasics)
library(feasts)
library(readr)
library(tidyr)
library(ggtime)
library(lubridate)
library(astsa)

#uaktualnienie wczytywania danych
plik <- "dane_S&P500.csv"

dane_surowe <- read_csv(
  plik,
  locale = locale(encoding = "UTF-8", decimal_mark = ",", grouping_mark = "."),
  col_types = cols(.default = "c")
)

# Wybór kolumn: data oraz kurs zamknięcia
dane <- dane_surowe %>%
  transmute(
    Data = as.Date(Data, format = "%d.%m.%Y"),
    `S&P500` = parse_number(
      Ostatnio,
      locale = locale(decimal_mark = ",", grouping_mark = ".")
    )
  ) %>%
  arrange(Data)
###############################################

# Filtrowanie zakresu dat
dane <- subset(dane, Data >= as.Date("2010-01-01") & Data <= as.Date("2026-01-31"))

#wykres surowych danych, do prezentacji
ggplot(dane, aes(x = Data, y = `S&P500`)) +
  geom_line(color = "black") +
  labs(title = "S&P 500", y = "Cena", x = "Data") +
  theme_minimal()

sp_ts <- ts(dane$`S&P500`, frequency = 4044)

dane_tsibble <- dane %>%
  as_tsibble(index = Data) %>%
  arrange(Data) %>%
  fill_gaps() %>%
  fill(`S&P500`, .direction = "down")  #uzupełnienie w weekendy i święta danymi danymi z ostaniego dnia sesji

dane_tsibble2 <- dane_tsibble %>%
  mutate(logdane = log(`S&P500`))

ggplot(dane_tsibble, aes(x = Data, y = `S&P500`)) +
  geom_line(color = "black") +
  labs(title = "S&P 500", y = "Cena", x = "Data") +
  theme_minimal()

#statystyki opisowe
cat("\n--- STATYSTYKI OPISOWE ---\n")
print(basicStats(dane_tsibble$`S&P500`))

################################################
plot(density(dane_tsibble$`S&P500`))

run_stationarity_analysis <- function(series, name) {
  x <- na.omit(series)
  dx <- na.omit(diff(x))
  cat("\n--- TESTY DLA:", name, "---\n")
  data.frame(
    Ujecie = c("Poziomy", "Poziomy", "1. Różnice", "1. Różnice"),
    Test = c("ADF", "KPSS", "ADF", "KPSS"),
    Stat = c(adf.test(x)$statistic, kpss.test(x)$statistic, adf.test(dx)$statistic, kpss.test(dx)$statistic),
    P_Val = c(adf.test(x)$p.value, kpss.test(x)$p.value, adf.test(dx)$p.value, kpss.test(dx)$p.value)
  )
}

print(run_stationarity_analysis(dane_tsibble2$logdane, "S&P 500")) #testy wskazują, że dane są stacjonarne po różnicowaniu
#adf.test(sp_ts) #dla testu ADF: p-value > 0.05 więc nie można odrzucić hipotezy zerowej zakładającej, że dane są niestacjonarne

options("scipen" = 100, "digits" = 5)

model_stl_mult <- dane_tsibble %>%
  model(
    STL(log(`S&P500`) ~ trend(window = 91) +
          season(period = 7, window = "periodic"))
  )

model_stl_mult %>%
  components() %>%
  autoplot() +
  labs(title = "Dekompozycja STL indeksu S&P 500 (dane zlogarytmowane)")

# Ostatnie 252 obserwacje (rok sesyjny)
model_stl_mult %>%
  components() %>%
  tail(252) %>%
  # Rysujemy tylko komponent sezonowy
  autoplot(season_7) +
  labs(
    title = "Sezonowość tygodniowa S&P 500 (dane zlogarytmowane)",
    subtitle = "Ostatnie 252 dni sesyjne",
    y = "Wartość",
    x = "Data"
  ) +
  theme_minimal()

model_stl_mult %>%
  components() %>%
  as_tibble() %>%
  mutate(dzien_tyg = wday(Data, label = TRUE, abbr = FALSE)) %>%
  ggplot(aes(x = dzien_tyg, y = season_7, fill = dzien_tyg)) +
  geom_boxplot() +
  labs(
    title = "Rozkład sezonowości S&P 500 wg dni tygodnia (dane zlogarytmowane)",
    y = "Wartość",
    x = "Dzień tygodnia"
  ) +
  theme_minimal() +
  guides(fill = "none")
model_stl <- dane_tsibble %>%
  model(
    STL(`S&P500` ~ trend(window = 63) +
          season(period = 7, window = 13),
        robust = TRUE)
  )

p_acf <- dane_tsibble %>%
  gg_tsdisplay(log(`S&P500`), plot_type = 'partial') +
  labs(title = "Korelogramy ACF/PACF dla S&P 500 (dane zlogarytmowane)")
print(p_acf)

dane_log_diff <- dane
dane_log_diff$`S&P500` <- c(NA, diff(log(dane$`S&P500`)))
dane_log_diff <- dane_log_diff[-1, ]

#tsibble
dane_tsibble3 <- dane_log_diff %>%
  arrange(Data) %>%
  as_tsibble(index = Data) %>%
  fill_gaps() %>%
  fill(`S&P500`, .direction = "down")

p_acf <- dane_tsibble3 %>%
  gg_tsdisplay(`S&P500`, plot_type = 'partial') +
  labs(title = "Korelogramy ACF/PACF dla S&P 500 (pierwsze różnice logarytmów)")
print(p_acf)

# Wykresy diagnostyczne
#Na podstawie wykresu ACF można stwierdzić, że dane są niestacjonarne
acf2(dane$`S&P500`, max.lag = 50)

sp_ts_diff <- diff(sp_ts)
sp_ts_logdiff <- diff(log(sp_ts))

#wykres danych różnicowanych
plot(sp_ts_diff)
plot(sp_ts_logdiff)

acf2(sp_ts_diff)
acf2(sp_ts_logdiff)

# Podział na zbiór uczący i testowy (zgodnie z najnowszymi danymi)
train_data <- dane_tsibble %>% filter_index(. ~ "2024-12-31")
test_data  <- dane_tsibble %>% filter_index("2025-01-01" ~ .)

fit <- train_data %>%
  model(
    arima    = ARIMA(log(`S&P500`)),
    ets      = ETS(log(`S&P500`)),
    naive    = NAIVE(log(`S&P500`)),
    rw_drift = RW(log(`S&P500`) ~ drift())
  )

# Generowanie prognoz na okres testowy (horyzont testowy)
fc <- fit %>% forecast(new_data = test_data)

# Obliczanie miar błędów Ex-post (MAE, RMSE, MAPE)
acc <- fc %>% accuracy(dane_tsibble) %>% dplyr::select(.model, RMSE, MAE, MAPE)

# Tabela porównawcza dla obu szeregów
comparison_table <- fit %>% glance() %>% dplyr::select(.model, AIC, AICc, BIC) %>% arrange(AICc)

print("--- KRYTERIA INFORMACYJNE: S&P 500 ---")
print(comparison_table)

##########################################################
fit %>% dplyr::select(rw_drift) %>% report()

# Pełna diagnostyka reszt dla najlepszych modeli
# Diagnostyka
fit %>%
  dplyr::select(arima) %>%
  gg_tsresiduals() +
  labs(title = "Diagnostyka reszt: ARIMA(1,1,1)")

calculate_full_metrics <- function(forecast_obj, actual_data, var_name) {

  # Przygotowanie danych rzeczywistych - wybieramy tylko Date i zmienną celu
  actual_df <- actual_data %>%
    as_tibble() %>%
    dplyr::select(Data, y = !!sym(var_name))

  # Przygotowanie prognoz i złączenie
  comp_data <- forecast_obj %>%
    as_tibble() %>%
    dplyr::select(Data, .model, .mean) %>%
    dplyr::rename(y_hat = .mean) %>%
    dplyr::inner_join(actual_df, by = "Data") %>%
    dplyr::filter(!is.na(y))

  # Obliczenia mierników
  comp_data %>%
    dplyr::group_by(.model) %>%
    dplyr::summarise(
      MSE  = mean((y_hat - y)^2),
      RMSE = sqrt(MSE),
      MAE  = mean(abs(y_hat - y)),
      MPE  = mean((y - y_hat) / y, na.rm = TRUE) * 100,
      # Dekompozycja współczynnika Theila:
      I1_Bias = (mean(y_hat) - mean(y))^2 / MSE,
      I2_Var  = (sd(y_hat) - sd(y))^2 / MSE,
      I3_Cov  = 2 * (1 - cor(y_hat, y)) * sd(y_hat) * sd(y) / MSE,
      Theil_I2 = sqrt(sum((y_hat - y)^2) / sum(y^2)),
      .groups = "drop"
    )
}

fc <- fc %>%
  mutate(y_hat = as.numeric(.mean))

# Wykonanie obliczeń
metrics <- calculate_full_metrics(fc, dane_tsibble, "S&P500")

cat("\n--- PEŁNE STATYSTYKI BŁĘDÓW: S&P 500 ---\n")
print(metrics)

# Ponowna estymacja na pełnych danych
final_fit_arima <- dane_tsibble %>%
  model(ARIMA_Final = ARIMA(log(`S&P500`)))

# Prognoza na 12 miesięcy w przód
final_fc_arima <- final_fit_arima %>% forecast(h = "12 months")

# Wykres wachlarzowy
final_fc_arima %>%
  autoplot(dane_tsibble %>% filter_index("2023 Jan" ~ .), level = c(80, 95)) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "darkgreen") +
  labs(title = "Prognoza Ex-ante",
       subtitle = "Model ARIMA | Przedziały ufności: 80% i 95%",
       y = "% r/r", x = "Data") +
  theme_minimal()

final_fit_rw <- dane_tsibble %>% model(ARIMA_Final = RW(log(`S&P500`) ~ drift()))

# Prognoza na 12 miesięcy w przód
final_fc_rw <- final_fit_rw %>% forecast(h = "12 months")

# Wykres wachlarzowy
final_fc_rw %>%
  autoplot(dane_tsibble %>% filter_index("2023 Jan" ~ .), level = c(80, 95)) +
  geom_hline(yintercept = 2.5, linetype = "dashed", color = "darkgreen") +
  labs(title = "Prognoza Ex-ante",
       subtitle = "Model Random Walk with drift | Przedziały ufności: 80% i 95%",
       y = "% r/r", x = "Data") +
  theme_minimal()

