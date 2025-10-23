sig <- read.delim(file = "C:\\Users\\roger\\Downloads\\Telegram Desktop\\resultados_Giancarlo\\Compare_results\\Significant_LOW_vs_HIGH", sep = "\t")
mat <- read.delim(file = "C:\\Users\\roger\\Downloads\\Telegram Desktop\\resultados_Giancarlo\\Compare_results\\All_exons.tab", sep = "\t")
diff <- read.delim(file = "C:\\Users\\roger\\Downloads\\Telegram Desktop\\resultados_Giancarlo\\Differential_tests/LH.tab", sep = "\t")
diff <- na.omit(diff)
events_in_high <- diff[diff$E.dPsi. > 0.4 & diff$MV.dPsi._at_0.95 > 0.4, ]
events_in_high <- events_in_high[order(events_in_high$E.dPsi., decreasing = TRUE), ]

events_in_low <- diff[diff$E.dPsi. < -0.4 & diff$MV.dPsi._at_0.95 > 0.4, ]
events_in_low <- events_in_low[order(events_in_low$E.dPsi., decreasing = FALSE), ]
