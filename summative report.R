## Read in data
library(nclSLR)
library(readr)
gexpr = read.csv("/Volumes/KINGSTON/Ch10Ex11.csv", header = FALSE)
dim(gexpr)

#Transposing the data set
# The genes are the rows and the columns are the individuals/tissues  
gexpr_t= t(gexpr)
dim(gexpr_t) 


# Question 1 
# a)
# First we need to compute the distance matrix 

d = dist(gexpr_t) 
hc_c = hclust(d, method="single")
plot(hc_c, cex=0.5, main="", sub="", xlab="")

#In general we observe two different clusters
#In order to decide where to cut the tree we take a look at the heights of consecutive fusions
# and we search for a big jump. It is easier to observe that if we plot the heights
hc_c$height
plot(hc_c$height, xlab="Fusion number", ylab="Fusion height", main = "Plot of fusion number against fusion height")

#We observe a bigger jump from around 45 to 47, so we can decide to cut the tree there 
cutree(hc_c, h=45.5)

# So we observe that the genes separate the two groups of healthy and diseased patients. 
#This is expected, as we would expect that the healthy patients would have more similar 
# gene expressions among them and quite different gene expressions with the ones that are 
# diseased 

# b)
#Now we will repeat part a using complete and average linkage in order to see if the result
# changes. First we do complete linkage:

hc_c_complete = hclust(d, method="complete")
plot(hc_c_complete, cex=0.5, main="", sub="", xlab="")
hc_c_complete$height
plot(hc_c_complete$height, xlab="Fusion number", ylab="Fusion height")

#Again we observe a jump between 48.5 and 52.5 so again we cut the dentogram at 49 and we end
# up with two clusters 
cutree(hc_c_complete, h=49)

#Similarly we can use the average linkage method 
hc_c_average = hclust(d, method="average")
plot(hc_c_average, cex=0.5, main="", sub="", xlab="") 
hc_c_average$height
plot(hc_c_average$height, xlab="Fusion number", ylab="Fusion height")

# Using average linkage and looking at the heights we see that there is a jump between 47 
# and 50, so we cut the dentogram at 47. 
cutree(hc_c_average, h = 47)

#So again we observe that the genes separate the sample into two groups of healthy and deseased
#So for this data set we see that we get the same result regardless of the linkage method 
# that we use 

# c)
# Now what we would like to know is which genes are the ones that differ across the healthy 
# and diseased groups, as this would give an insight of the genes causing the problem 
par(mfrow = c(1, 2))
colMeans(gexpr_t) 
plot(1:1000, colMeans(gexpr_t), col = ifelse(colMeans(gexpr_t)>0.5 ,"red", "black"))
# Just be looking at the column means of the genes in the 40 patients, we observe that some 
#specific genes have quite a higher mean compared to the most other genes. Now what we would 
#like to do is find which are these genes. So the red colored genes are the ones who probably
# affect the health of the patients 

healthy = gexpr_t[1:20,]
dim(healthy)

diseased = gexpr_t[21:40,]
dim(diseased) 

colMeans(healthy)
colMeans(diseased)

plot(colMeans(healthy), colMeans(diseased), col = ifelse(colMeans(diseased)>1 ,"red", "black"))

which(colMeans(diseased)>1)
which(colMeans(gexpr_t)>0.5)

setequal(which(colMeans(diseased)>1),which(colMeans(gexpr_t)>0.5))

pca1 = prcomp(x=gexpr)
summary(pca1)
# Plot the first PC against the second PC
plot(pca1$x[,1], pca1$x[,2], xlab="First PC", ylab="Second PC")
# Add labels representing the cancer types
text(pca1$x[,1], pca1$x[,2], labels=rownames(gexpr), cex=0.7, pos=3)
# We can check if the points are the same as the ones found above: 
pca1vector = pca1$x[,1]
pca1vector[pca1vector < -4]
problematic_genes = which(pca1vector < -4)
problematic_genes

#So the genes found above are the ones who we have found above and this means that we have 
# identify which genes are the ones that are more expressed in the diseased people, so they 
# could be the ones that affect their health


# Question 2
# a)
kmeans_function = function(x,k){
#making sure that the data entered are in matrix form
x = as.matrix(x)
#choosing the initial means randomly
initial_means = x[sample(nrow(x),size=k,replace=FALSE),]
old_clusters = rep(1,nrow(x))
clusters = rep(2, nrow(x))
# Create an empty matrix which will later be populated with the distance of the matrix 
#points to the means
distance_from_means = matrix(nrow = nrow(x), ncol = k)
while(all(old_clusters != clusters)){
  # Set the old clusters equal to the clusters found
  old_clusters = clusters
  for(i in 1:nrow(x)){
    for(j in 1:k){
      # Find the euclidian distance for each each row to the k randomly selected initial means
      distance = sqrt(sum((x[i,] - initial_means[j,])^2))
      # populate the columns of the matrix with the distances
      distance_from_means[i,j] = distance
    }
}
# Assign clusters based on which distance is closer to the means 
clusters = apply(distance_from_means, 1, which.min)
# Find the means of the different clusters
centers = apply(x, 2, tapply, clusters, mean)
# Set the new means as the initial means in order to run the for loop with the new means
initial_means = centers
    
}
#Create a for loop, that goes over the rows of the matrix finds the total within variation
within_clusters_sum_of_squares = 0
for(g in 1:nrow(x)){
  within_clusters_sum_of_squares =  within_clusters_sum_of_squares + sum((x[g,] - centers[clusters[g],])^2)
}
#Create a list that contains the cluster partition, the centers and the witihn cluster sum of squares  
y = list(clusters,centers, within_clusters_sum_of_squares) 
# Put headers on the list elements 
names(y) = c("cluster partition", "cluster means", "within clusters sum of squares")
# Return the list
return(y)
}

km = kmeans_function(gexpr_t,2)
km[[1]]


# b)
## Load nclSLR package
library(nclSLR)
## Read in data
data(USArrests, package="nclSLR")
## Check size
dim(USArrests)
kmeans_function(USArrests[,1:4],3)
kmeans(USArrests[,1:4],3)$tot.withinss
kmeans(USArrests[,1:4],3)$centers

#run the function 50 times 
run = replicate(50, kmeans_function(USArrests[,1:4],3), simplify = FALSE) 

#run[[3]][3]
#is.list(run)

# Extract the third element 
ssw = sapply(run, "[[", 3)
ssw
# We find which run has generated the smallest within clusters sum of squares
which.min(ssw)
#We get that the first run generated the smallest ssw so we extract that run to compare it with the 
#built in algorith
min_run = run[which.min(ssw)]
min_run

#Now we run the algoritm and extract the cluster partition, centers and total within variation in
#order to compare if we get the same results 
km = kmeans(USArrests[,1:4],3)
km$cluster
km$centers
km$tot.withinss

#We end up with the same result so everything is fine

# c)
# No, it is possible that the procedure will not give the same result. This happens because the 
# k means algorithm is guaranteed to converge to a minumun value of SSW, but this may be the local
# and not the global minimum. So depending on the starting point the kmeans algorithm may find a 
#partition with a local mimimun SSW abd stop there, while the procedure in part 2, which may have 
# a different starting point could find the global maximums and thus the results could be different. 
# Even if the algorithm is run with different initial starting points, this may find the global 
#minimum, but again there is no guarantee that this will always happen. 

# Question e
## Load the nclSLR package
library(nclSLR)
## Load the data
data(diabetes)
## Check size
dim(diabetes)
head(diabetes , 3) 
diabetes[,1:10] = scale(diabetes[,1:10])

# a)
# Now we need to define the training and test sets 

train_set = diabetes[1:350,]
test_set = diabetes[351:442, ]

# b) 

lsq_train = lm(dis ~ ., data= train_set)

## Compute fitted values for test data:
yhat_test = predict(lsq_train, test_set)
head(yhat_test)

## Compute test error
test_error = mean((test_set$dis - yhat_test)^2)
test_error

#The test error is 2842.256

# c) 
# From practical 5 we see that the model selected by subset selection includes sex, bmi, map,
# tc, ldl and ltg. So we fit the model with just these variables 

lsq_train_ss = lm(dis ~ sex + bmi + map + tc + ldl + ltg , data = train_set)
## Compute fitted values for test data:
yhat_test_ss = predict(lsq_train_ss, test_set)
head(yhat_test_ss)

## Compute test error
test_error = mean((test_set$dis - yhat_test_ss)^2)  
test_error

# The test error is 2800.964. As expected it has reduced 

# d) 
# i)
# I)
library(pls)
## Fit model using PLS
plsr_fit = plsr(dis ~ ., data = train_set, scale=FALSE)
summary(plsr_fit)

# We see that the first 4 transformed variables explain approximatelly 72% of the
# variation in X and 50% of the variation in Y. We also observe that in general the variation 
#explained by the principal components does not really changes after we have included the third pc.

## Fit model, applying 10-fold cross validation:
plsr_cv_fit = plsr(dis ~ ., data=train_set, scale=FALSE, validation = "CV")
## Plot the cross-validation scores:
plot(plsr_cv_fit, plottype="validation", legend="topright", val.type="MSEP")
# The optimal number of variables are 2

# II)
yhat_m2 = predict(plsr_fit, test_set, ncomp = 2)
test_error_1 = mean((test_set$dis - yhat_m2)^2)
test_error_1
# The test error is 2929.477


# ii) 
plsr_fit_full = plsr(dis ~ ., data = diabetes, scale=FALSE)
## Examine the directions defined by PLS
plsr_load_full = loadings(plsr_fit_full)
## Print the directions to the screen
C_full = unclass(plsr_load_full)
C_full[,1]

plsr_cv_fit_full = plsr(dis ~ ., data= diabetes, scale=FALSE, validation = "CV")
plot(plsr_cv_fit_full, plottype="validation", legend="topright", val.type="MSEP")

# Now we can find the test error. First we need to fit the model with the 2 variables 
yhat_full_m2 = predict(plsr_fit_full, diabetes, ncomp = 2) 
test_error_2 = mean((diabetes$dis - yhat_full_m2)^2)
test_error_2
# the test error is 2915.573

# The optimal number of variables is 2 
coef(plsr_fit_full, intercept=TRUE, ncomp=2)

# Fit using least squares
lsq_fit = lm(dis ~ ., data=diabetes)
lsq_fit$coefficients

# We can think the first axis as a contrast between hdl and all the other variables 

# e) We observe that we get the smallest test error when we use just the 7 variables identified 
# by subset selection 



