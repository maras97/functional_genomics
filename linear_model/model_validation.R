## read norm-Table of 48h with normalized histon-modification counts in promoter regions (already normalized by control)
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/linear_regression/")
norm_48h <- read.table("norm_48h_new2.txt",header = T)
m_48h <- read.table("med_RNA_48h.txt")
names_48h <- m_48h[,1]
m_48h <- as.numeric(m_48h[,2])
names(m_48h) <- names_48h


## randomly select "training" (1), "validation" (2), "test" (3) for every index/row
random <- sample(x= 1:3, size=22369, prob=c(0.5,0.25,0.25), replace=T)

training_48h <- c()
val_48h <- c()
test_48h <- c()

training_48h <- m_48h[(random == 1)]
val_48h <- m_48h[(random == 2)]
test_48h <- m_48h[(random == 3)]


## generalized linear model (glm) fitting for expression values dependent on histon modifications 
norm_48h <- as.matrix(norm_48h)
norm_glm_48h <- norm_48h[(random == 1),]
data_48h  <- data.frame(cbind(training_48h,norm_glm_48h))

model_48h <- lm(training ~ H3.3 + H3K27ac + H3K27me3 + H3K36me3 + H3K4me1 +  H3K4me2 +  
               H3K4me3 + H3K79me2 + H3K9ac + H3K9me3 + H3 + Hdac1 + 0,data = data)

## read norm-Table of ESC with normalized histon-modification counts in promoter regions (already normalized by control)
setwd("/Users/marasteiger/Documents/Uni/Functional Genomics/linear_regression/")
norm_ESC <- read.table("norm_ESC_new2.txt",header = T)
m_ESC <- read.table("med_RNA_ESC.txt")
names_ESC <- m_ESC[,1]
m_ESC <- as.numeric(m_ESC[,2])
names(m_ESC) <- names_ESC

training_ESC <- c()
val_ESC <- c()
test_ESC <- c()

training_ESC <- m_ESC[(random == 1)]
val_ESC <- m_ESC[(random == 2)]
test_ESC <- m_ESC[(random == 3)]


## predict ESC-RNA Counts by using the model trained on 48h expression data with ESC-histone modifications as input 
norm_test_ESC = norm_ESC[(random == 3),]
pr_ESC <- predict(model_48h, newdata = data.frame(H3.3 = norm_test_ESC[,1],
                                           H3K27ac = norm_test_ESC[,2],
                                           H3K27me3 = norm_test_ESC[,3],
                                           H3K36me3 = norm_test_ESC[,4],
                                           H3K4me1 = norm_test_ESC[,5],
                                           H3K4me2 = norm_test_ESC[,6],
                                           H3K4me3 = norm_test_ESC[,7],
                                           H3K79me2 = norm_test_ESC[,8],
                                           H3K9ac = norm_test_ESC[,9],
                                           H3K9me3 = norm_test_ESC[,10],
                                           H3 = norm_test_ESC[,11],
                                           Hdac1 = norm_test_ESC[,12]))
                                           
                                           
plot((pr_ESC-test_ESC))

mean(abs(pr-test))
#1.703608


hist((pr_ESC-test_ESC), breaks=seq(-13,18, l = 50))