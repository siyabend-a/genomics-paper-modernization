# Predictors and target from supplementary information (https://www.pnas.org/doi/10.1073/pnas.0904863106)
tfas <- read.table("sd2.txt", header = TRUE, sep = "\t")
rpkm <- read.table("sd1.txt", header = TRUE, sep = "\t")
tfas_c <- tfas[, -c(1, 2)]
svd_result <- svd(tfas_c, nu = min(nrow(tfas_c), ncol(tfas_c)), nv = min(nrow(tfas_c), ncol(tfas_c)))
X <- svd_result$u %*% diag(svd_result$d)
logY <- log(rpkm[, 3] + 0.1) # To avoid taking the logarithm of zero, a small positive constant is added to Yi. 
pcr_model <- lm(logY ~ X)

library(ggplot2)
plot_data <- data.frame(RNA_seq_ESC_predicted = predict(pcr_model), RNA_seq_ESC_exp = logY)
ggplot(plot_data, aes(y = RNA_seq_ESC_predicted, x = RNA_seq_ESC_exp)) +
                  geom_point(size=0.001) +
                  annotate("text", x = Inf, y = -Inf, label = paste0("R^2 == ", summary(pcr_model)$r.squared), 
                  hjust = 1.1, vjust = -1.1, size = 5, parse = TRUE, fontface = "bold")