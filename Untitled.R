library("ggplot2") 
library("ggfortify") 
library("cluster") 
library("Biobase") 
library("ape") 
library('RColorBrewer') 
library("caTools") 
library("randomForest") 
library("caret") 
library("ROCR") 
library("lattice")

d <- read.csv("GEOG0610.csv") 
protein.id <- d[1:2320,1] 
country <- rep(c("Kenya", "Mali", "PNG", "Solomon", "Naive"), 
               c(110, 32, 50, 26, 20)) 
continent <- rep(c("Africa", "Asia", "AU"), c(142, 76, 20))

#################### Negative Control ########################## 
noDNA <- d[2322:2369,] 
boxplot(noDNA) 
#Remove samples with high varibility based on boxplot 
rem <- c("Kenya_20", "Kenya_21", "Kenya_42", "Kenya_112", "Kenya_115", "Mali_19", 
         "Mali_25", "Mali_25", "Mali_27", "Mali_27", "Mali_30", "Mali_31", "PNG_10", 
         "PNG_30", "PNG_40", "PNG_41", "PNG_44", "PNG_47", "PNG_49", "SOLOMON_1", 
         "SOLOMON_7", "SOLOMON_10", "SOLOMON_21") 
noDNA <- noDNA[, !names(noDNA) %in% rem]

#################### Continuous variables ###################### 
####### Data Transformation 
data <- d[1:2320, !names(d) %in% rem] 
data <- data[,-1]

##Subtract noDNA mean vlaue from sample data points for each sample 
centercolmeans <- function(x, y) { 
  xcenter = colMeans(x) 
  y - rep(xcenter, rep.int(nrow(y), ncol(y))) 
} 
data <- centercolmeans(noDNA[,-1], data)

##Normalization: log2 
prtn.norm <- log2(data) 
prtn.norm[is.na(prtn.norm)] = 0

####### Data analysis 
## 1.Hierarchical clustering (pearson) 
#Protein based 
prtn.cor <- cor(prtn.norm, method = "pearson") 
prtn.cor[is.na(prtn.cor)] = 0 p
rtn.cor.dist <- as.dist(1-prtn.cor) 
prtn.tree <- hclust(prtn.cor.dist, method = "average") 
#Country based (samples) clust <- t(prtn.norm) 
cntry.cor <- cor(clust, method = "pearson") 
cntry.cor[is.na(cntry.cor)] = 0 
cntry.dist <- as.dist(1-cntry.cor) 
cntry.tree <- hclust(cntry.dist, method = "average") 
#Heatmap 
col <- colorRampPalette(brewer.pal(10, 'RdBu'))(256) 
col <- rev(col)#reverse colors as indicated above 
heatmap_data <- as.matrix(prtn.norm) 
colnames(heatmap_data) <- country 
heatmap_factor <- as.factor(colnames(heatmap_data)) 
color <- rep(c("red", "blue", "yellow", "green", "black"), 
             c(110, 32, 50, 26, 20)) 
heatmap(as.matrix(prtn.norm), scale='row', col = col,
        Rowv=as.dendrogram(cntry.tree),
        Colv=as.dendrogram(prtn.tree),
        main = "Sample countries data",
        ColSideColors=color)

## 2.PCoA pcoa(prtn.cor.dist) 
biplot.pcoa(pcoa(prtn.cor.dist))

## 3.PCA 
s <- scale(t(prtn.norm)) 
s[is.na(s)] <- 0

prtn.pc <- prcomp(s, center = TRUE) 
s <- cbind(country, s) 
s <- cbind(continent, s) 
autoplot(prtn.pc, data = s, colour = "country") 
autoplot(prtn.pc, data = s, colour = "continent")

## 4.Differential analysis (t-test) 
kenya <- t(prtn.norm[,1:110]) 
mali <- t(prtn.norm[,111:142]) 
PNG <- t(prtn.norm[,143:192]) 
solomon <- t(prtn.norm[,193:218]) 
native <- t(prtn.norm[,219:238]) 
#Variance test first, "not equal" for all combo 
ttest <- function(dataset1, dataset2){ 
  p_value <- c(NULL) 
  for (i in 1:2320){ 
    p_value[i] <- t.test(dataset1[,i], 
                         dataset2[,i], 
                         alternative = "two.sided", 
                         var.equal = FALSE)$p.value 
  } 
  return (p_value) 
} 
#FDR adjust 
prnt_d <- c(which(p.adjust(p_value, method = "fdr") < 0.05)) 
#T-tests: differentially expressed proteins 
kenya_naive <- ttest(kenya, naive)

## 5.Classifier (randomForest) 
classify <- function(data, factor){
  set.seed(1234) 
  split <- sample.split(data, SplitRatio = 0.8) 
  data.training <- data[split,] 
  data.testing <- data[!split,] 
  data_model <- randomForest(factor ~ ., 
                             data = data.training, 
                             na.action = na.exclude) 
  data_testing_pred <- predict(data_model, 
                               newdata = data.testing) 
  data_testing_matrix <- confusionMatrix(data_testing_pred, 
                                         reference = data.testing$factor) 
  print(data3_testing_matrix)

#ROC curve 
data.pred <- predict(data_model, type = "prob", 
                     newdata = data.testing)[,2] 
#ROCR package 
data_pred <- prediction(data.pred, data.testing$factor)#labels 
data.perf <- performance(data_pred, "tpr", "fpr") 
#plot ROC curve 
plot(data4.perf,
     main = "ROC curve for Malaria exposed continents",
     col = 2, lwd = 2) 
abline(a=0,b=1,lwd=2,lty=2,col="gray")
return (data_testing_matrix) 
} 
#Five Countries 
class <- data.frame(t(data)) 
class <- cbind(country, class) 
colnames(class)[1] <- "Country" 
class$Country <- as.factor(class$Country) 
country_matrix <- classify(class, "Country") 
#Without Naive
data1 <- data.frame(t(data[,1:218])) 
data1 <- cbind(country[1:218], data1) 
colnames(data1)[1] <- "Malaria" 
data1$Malaria <- as.factor(data1$Malaria) 
malaria_matrix <- classify(data1, "Malaria") 
#Kenya vs. Mali d
ata2 <- data.frame(t(data[,1:142])) 
data2 <- cbind(country[1:142], data2) 
colnames(data2)[1] <- "africa" 
data2$africa <- as.factor(data2$africa) 
africa_matrix <- classify(data2, "africa") 
#PNG vs. Solomon 
data3 <- data.frame(t(data[,143:218])) 
data3 <- cbind(country[143:218], data3) 
colnames(data3)[1] <- "asia" 
data3$asia <- as.factor(data3$asia) 
aisa_matrix <- classify(data3, "asia") 
#Afica vs. Asia 
data4 <- data.frame(t(data[,1:218])) 
data4 <- cbind(continent[1:218], data4) 
colnames(data4)[1] <- "Countries" 
data4$Countries <- as.factor(data4$Countries) 
continent_matrix <- classify(data4, "Countries")
