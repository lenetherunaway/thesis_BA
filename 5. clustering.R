library(tidyr)
library(dtwclust)
library(dplyr)
library(ggplot2)
library(reshape)
co2 <- read.csv('/Users/elena/Documents/Диплом/Кластеризация/co2_to_cluster.csv')
#http://rstudio-pubs-static.s3.amazonaws.com/398402_abe1a0343a4e4e03977de8f3791e96bb.html
co2_long <- gather(co2)
co2_long$time<-rep(1992:2014,92)
#co2_long  %>% 
#  ggplot(aes(x= time, y= value, color= key)) +
#  geom_line( size=0.2) +
#  ggtitle("Control chart sequences") + 
#  facet_wrap(~ key , scales = 'free_x', nrow= 10) 
co2_list <- as.list(utils::unstack(co2_long, value ~ key))

co2_list_z <- dtwclust::zscore(co2_list)

#hierarchical clustering with 10% window size for up to k=10 clusters
cluster_dtw_h <-list()
for (i in 2:10)
{
  cluster_dtw_h[[i]] <- tsclust(co2_list_z, type = "h", k = i, distance = "dtw", control = hierarchical_control(method = "ward.D"), seed = 390, preproc = NULL, args = tsclust_args(dist = list(window.size =5L)))
}
#method  should be one of “ward.D”, “ward.D2”, “single”, “complete”, “average”, “mcquitty”, “median”, “centroid”, “all”
# take a look at the object
cluster_dtw_h[[3]]

cluster_dtw_h[[3]]@clusinfo

plot(cluster_dtw_h[[2]])
plot(cluster_dtw_h[[3]])
plot(cluster_dtw_h[[2]], type='sc') #Cluster members
plot(cluster_dtw_h[[3]], type='sc') #Cluster members
plot(cluster_dtw_h[[4]], type='sc') #Cluster members
plot(cluster_dtw_h[[5]], type='sc') #Cluster members

cluster_dtw_h <- tsclust(co2, k=2L:7L, distance = "dtw", centroid = "dba", seed=94L)
names <- paste0("k_", 2L:7L)
sapply(cluster_dtw_h, cvi, type="internal")
hc.l2 <- tsclust(co2, type = "hierarchical",
                 k = 4L, trace = TRUE,
                 control = hierarchical_control(method = "all"))

# Plot the best dendrogram according to variation of information
plot(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]])


#https://petolau.github.io/TSrepr-clustering-time-series-representations/
library(TSrepr)
library(ggplot2)
library(data.table)
library(cluster)
library(clusterCrit)

clusterings <- lapply(c(2:7), function(x)
  pam(co2, x))
DB_values <- sapply(seq_along(clusterings), function(x) 
  intCriteria(as.matrix(t(co2)), as.integer(clusterings[[x]]$clustering),
              c("Davies_Bouldin")))

ggplot(data.table(Clusters = 2:7, DBindex = unlist(DB_values)),
       aes(Clusters, DBindex)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw()

#2 кластера оптимально
# prepare data for plotting
data_plot <- data.table(melt(data.table(class = as.factor(clusterings[[1]]$clustering),
                                        co2)))
data_plot[, Time := rep(1:ncol(co2), each = nrow(co2))]
data_plot[, ID := rep(1:nrow(co2), ncol(co2))]

# prepare medoids
centers <- data.table(melt(clusterings[[6]]$medoids))
setnames(centers, c("Var1", "Var2"), c("class", "Time"))
centers[, ID := class]

# plot the results
ggplot(data_plot, aes(Time, value, group = ID)) +
  facet_wrap(~class, ncol = 2, scales = "free_y") +
  geom_line(color = "grey10", alpha = 0.65) +
  geom_line(data = centers, aes(Time, value),
            color = "firebrick1", alpha = 0.80, size = 1.2) +
  labs(x = "Time", y = "Load (normalised)") +
  theme_bw()


#https://uc-r.github.io/hc_clustering#kmeans
#https://uc-r.github.io/kmeans_clustering choose number of clusters
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend)

d <- dist(co2, method='dtw') 
#hc1 <- hclust(d, method = "complete" )
hc1 <- agnes(d, method = "kmeans" )
plot(hc1)
hc1$ac

hc2 <- agnes(d, method = "average" )
plot(hc2, hang = -1)
hc2$ac

hc3 <- agnes(d, method = "pam" )
plot(hc3, hang = -1)
hc3$ac

hc4 <- agnes(d, method = "ward" )
plot(hc4, hang = -1)
hc4$ac

# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(co2, method = x)$ac
}

map_dbl(m, ac)


#https://medium.com/@sametgirgin/hierarchical-clustering-model-in-5-steps-with-python-6c45087d4318
#https://rpubs.com/esobolewska/dtw-time-series
#https://cran.r-project.org/web/packages/dtwclust/vignettes/dtwclust.pdf
#https://towardsdatascience.com/10-tips-for-choosing-the-optimal-number-of-clusters-277e93d72d92
#https://stats.stackexchange.com/questions/131281/dynamic-time-warping-clustering
library(tidyverse)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(GGally)
library(ggiraphExtra)
library(knitr)
library(kableExtra)

d <- dist(t(co2), method = "dtw") 
d_eu <- dist(t(co2), method = "euclidean") 
km <- kmeans(t(co2), centers=3, nstart=30)
p <- pam(t(co2), k=3)
clar <- clara(t(co2), k=3)
set.seed(42)

fviz_nbclust(t(co2), kmeans, method = "wss", diss = d, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
fviz_nbclust(t(co2), kmeans, method = "wss", diss = d_eu, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
par(mfrow=c(1,2))
fviz_nbclust(t(co2), pam, method = "wss", diss = d, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
fviz_nbclust(t(co2), pam, method = "wss", diss = d_eu, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))

fviz_nbclust(t(co2), clara, method = "wss", diss = d, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
fviz_nbclust(t(co2), clara, method = "wss", diss = d_eu, k.max = 10) + theme_minimal() + ggtitle("Метод локтя") +ylab("Общая внутренняя сумма квадратов")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))

gap_stat <- clusGap(t(co2), FUN = kmeans, dist=d, nstart = 30, K.max = 10, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("Статистика разрыва")+ylab("Статистика разрыва (k)")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))

gap_stat <- clusGap(t(co2), FUN = pam, nstart = 30, K.max = 10, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("Статистика разрыва")+ylab("Статистика разрыва (k)")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))

fviz_nbclust(t(co2), kmeans, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("Силуэтный метод")+ylab("Средняя ширина ")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
fviz_nbclust(t(co2), pam, method = "silhouette", k.max = 10) + theme_minimal() + ggtitle("Силуэтный метод")+ylab("Средняя ширина ")+xlab("Число кластеров")+theme(plot.title = element_text(hjust = 0.5))
intern <- clValid(t(co2), nClust = 2:24, 
                  clMethods = c("hierarchical","kmeans","pam"), metric = "euclidean", validation = "internal")


# Compute dissimilarity matrix with euclidean distances

# Hierarchical clustering using Ward's method
res.hc <- hclust(d, method = "complete" )
# Cut tree into 5 groups
grp <- cutree(res.hc, k = 2)
# Visualize
plot(res.hc, cex = 0.6) # plot tree
rect.hclust(res.hc, k = 2, border = 2:2) # add rectangle
