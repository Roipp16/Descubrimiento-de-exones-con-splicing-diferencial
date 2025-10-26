library(ggplot2)
library(dplyr)

### Cargar datos
low_vs_interm_exons <- read.delim("C:/Users/roger/Downloads/Telegram Desktop/resultados_Giancarlo/Compare_results/Significant_LOW_vs_INTERM", header = TRUE, sep = "\t")
interm_vs_high_exons <- read.delim("C:/Users/roger/Downloads/Telegram Desktop/resultados_Giancarlo/Compare_results/Significant_INTERM_vs_HIGH", header = TRUE, sep = "\t")
low_vs_high_exons <- read.delim("C:/Users/roger/Downloads/Telegram Desktop/resultados_Giancarlo/Compare_results/Significant_LOW_vs_HIGH", header = TRUE, sep = "\t")
all_exons <-  read.delim("C:/Users/roger/Downloads/Telegram Desktop/resultados_Giancarlo/Compare_results/All_exons.tab", header = TRUE, sep = "\t")


### Estadísticas descriptivas para cada comparación
get_stats <- function(data, label){
  data.frame(
    Comparison = label,
    mean_LENGTH = mean(data$LENGTH),
    median_LENGTH = median(data$LENGTH),
    sd_LENGTH = sd(data$LENGTH)
  )
}

stats_low_high   <- get_stats(low_vs_high_exons,   "LOW vs HIGH")
stats_low_interm <- get_stats(low_vs_interm_exons, "LOW vs INTERM")
stats_interm_high<- get_stats(interm_vs_high_exons,"INTERM vs HIGH")
stats_bg         <- get_stats(all_exons,          "Background")

de_stats <- rbind(stats_low_high, stats_low_interm, stats_interm_high)
print(de_stats)
print(stats_bg)


### Tests estadísticos
test_low_high <- wilcox.test(low_vs_high_exons$LENGTH, all_exons$LENGTH)
test_low_interm <- wilcox.test(low_vs_interm_exons$LENGTH, all_exons$LENGTH)
test_interm_high <- wilcox.test(interm_vs_high_exons$LENGTH, all_exons$LENGTH)

cat("Wilcoxon test LOW vs HIGH vs Background:\n")
print(test_low_high)
cat("\nWilcoxon test LOW vs INTERM vs Background:\n")
print(test_low_interm)
cat("\nWilcoxon test INTERM vs HIGH vs Background:\n")
print(test_interm_high)


### Preparar datos para gráfico combinado
low_high_df <- data.frame(LENGTH = log10(low_vs_high_exons$LENGTH), Group = "LOW vs HIGH")
low_interm_df <- data.frame(LENGTH = log10(low_vs_interm_exons$LENGTH), Group = "LOW vs INTERM")
interm_high_df <- data.frame(LENGTH = log10(interm_vs_high_exons$LENGTH), Group = "INTERM vs HIGH")
bg_df <- data.frame(LENGTH = log10(all_exons$LENGTH), Group = "Background")

plot_data <- rbind(low_high_df, low_interm_df, interm_high_df, bg_df)

### Gráfico combinado
p_all <- ggplot(plot_data, aes(x = LENGTH, color = Group)) +
  geom_density(alpha = 0.7) +
  labs(title = "Distribución de la longitud de exones (log10)",
       x = "log10(Exon LENGTH)", y = "Densidad") +
  theme_minimal() +
  scale_color_manual(values = c(
    "Background" = "grey40",
    "LOW vs HIGH" = "#E41A1C",
    "LOW vs INTERM" = "#377EB8",
    "INTERM vs HIGH" = "#4DAF4A"
  ))

print(p_all)


############ filtrado de significancia ###########
### se lee la tabla con eventos
diff <- read.delim(file = "C:\\Users\\roger\\Downloads\\Telegram Desktop\\resultados_Giancarlo\\Differential_tests/LH.tab", sep = "\t")
diff <- na.omit(diff)

### se aplican filtros para encontrar los eventos más importantes en el grupo High
events_in_high <- diff[diff$E.dPsi. > 0.4 & diff$MV.dPsi._at_0.95 > 0.4, ]
events_in_high <- events_in_high[order(events_in_high$E.dPsi., decreasing = TRUE), ]

### se aplican filtros para encontrar los eventos más importantes en el grupo Low
events_in_low <- diff[diff$E.dPsi. < -0.4 & diff$MV.dPsi._at_0.95 > 0.4, ]
events_in_low <- events_in_low[order(events_in_low$E.dPsi., decreasing = FALSE), ]
