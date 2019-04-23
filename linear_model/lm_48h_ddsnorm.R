## read norm-Table with normalized histon-modification counts in promoter regions (already normalized by control)
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/linear_regression/")
norm <- read.table("norm_48h_new2.txt",header = T)
m <- read.table("med_RNA_48h.txt")
names <- m[,1]
m <- as.numeric(m[,2])
names(m) <- names



## randomly select "training" (1), "validation" (2), "test" (3) for every index/row
random <- sample(x= 1:3, size=22369, prob=c(0.5,0.25,0.25), replace=T)

training <- c()
val <- c()
test <- c()

training <- m[(random == 1)]
val <- m[(random == 2)]
test <- m[(random == 3)]




# generalized linear model (glm) fitting for expression values dependent on histon modifications 
norm <- as.matrix(norm)
norm_glm <- norm[(random == 1),]
data  <- data.frame(cbind(training,norm_glm))

model <- lm(training ~ H3.3 + H3K27ac + H3K27me3 + H3K36me3 + H3K4me1 +  H3K4me2 +  
               H3K4me3 + H3K79me2 + H3K9ac + H3K9me3 + H3 + Hdac1 + 0,data = data)
#model <- glm(training ~ norm_glm + 0, family = gaussian) 
#m2 <- lm(training ~ norm_glm + 0) 
#beta <- model$residuals
coefs <- model$coefficients

summary(model)
q <- summary(model)
rsq <- q[["r.squared"]]
#0.8304184

## greatest beta-coefficients and corresponding Histon modifications
max.coefs <- tail(sort(coefs),n=3)
# norm_glmH3K4me1 norm_glmH3K79me2  norm_glmH3K4me3

## smallest beta-coefficients and corresponding Histon modifications
min.coefs <- head(sort(coefs),n=3)
# norm_glmH3K9ac norm_glmH3K4me2 norm_glmH3K9me3  

library(affy)
scatter1 <- smoothScatter((norm[,12]),(m),xlab = "Hdac1", ylab = "normalized RNA-Counts")
scatter1 <- smoothScatter((norm[,7]),(m),xlab = "H3K4me3", ylab = "normalized RNA-Counts")

scatter1 <- smoothScatter((norm[,2]),(m),xlab = "H3K27ac", ylab = "normalized RNA-Counts")
cor(norm[,7],m)
#cor((counts.f[,1]),(norm[,2]))

library(boot)
model2 <- glm(training ~ norm_glm + 0, family = gaussian) 
diag <- glm.diag(model2)
glm.diag.plots(model2)


# Multiple R-squared:  0.8304,	Adjusted R-squared:  0.8302 
# F-statistic:  4522 on 12 and 11081 DF,  p-value: < 2.2e-16

norm_test = norm[(random == 3),]
pr <- predict(model, newdata = data.frame(H3.3 = norm_test[,1],
                                           H3K27ac = norm_test[,2],
                                           H3K27me3 = norm_test[,3],
                                           H3K36me3 = norm_test[,4],
                                           H3K4me1 = norm_test[,5],
                                           H3K4me2 = norm_test[,6],
                                           H3K4me3 = norm_test[,7],
                                           H3K79me2 = norm_test[,8],
                                           H3K9ac = norm_test[,9],
                                           H3K9me3 = norm_test[,10],
                                           H3 = norm_test[,11],
                                           Hdac1 = norm_test[,12]))


plot((pr-test), pch = 20)

mean(abs(pr-test))
#1.683918

hist((pr-test), breaks=seq(-15,16, l = 50))


# down <- which(m<3)
# up <- which(m>=3)
# library(rtracklayer)
# BiocManager::install("rtracklayer", version = "3.8")
# export(prot.gene.prom[up],"up.bed")
# export(prot.gene.prom[down],"down.bed")

