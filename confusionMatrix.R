library(mlearning);
library(mlbench);
data("Glass", package = "mlbench")
Glass$Type <- as.factor(paste("Glass", Glass$Type))

summary(glassLvq <- mlLvq(Type ~ ., data = Glass));
(glassConf <- confusion(predict(glassLvq, Glass, type = "class"), Glass$Type))

plot(glassConf) # Image by default
