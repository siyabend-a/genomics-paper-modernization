total <- read.table("total2.txt", header = TRUE, sep = "\t")
tfas = total[, c(2:13)]
rpkm = total[, c(14)]
head(tfas)
#rpkm <- read.table("final_rpkm.csv", header = TRUE, sep = ",")
#tfas_c <- tfas[, -c(1, 2)]
#tfas <- log(tfas + 0.1)
svd_result <- svd(tfas, nu = min(nrow(tfas), ncol(tfas)), nv = min(nrow(tfas), ncol(tfas)))
X <- svd_result$u %*% diag(svd_result$d)
#logY <- rpkm
logY <- log(rpkm + 0.1) # To avoid taking the logarithm of zero, a small positive constant is added to Yi. 
pcr_model <- lm(logY ~ X)

library(ggplot2)
plot_data <- data.frame(RNA_seq_ESC_predicted = predict(pcr_model), RNA_seq_ESC_exp = logY)
ggplot(plot_data, aes(y = RNA_seq_ESC_predicted, x = RNA_seq_ESC_exp)) +
                  geom_point(size=0.000001) +
                  annotate("text", x = Inf, y = -Inf, label = paste0("R^2 == ", summary(pcr_model)$r.squared), 
                  hjust = 1.1, vjust = -1.1, size = 5, parse = TRUE, fontface = "bold")
