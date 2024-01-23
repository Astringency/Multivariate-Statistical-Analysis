---
title: "多元正态分布、参数估计与假设检验"
author: "Zhang Xifeng"
date: "`r Sys.Date()`"
documentclass: ctexart
output: rticles::ctex
---

# 1 多元正态分布及参数的估计

## 1.1导入数据

```{r}
data <- iris
head(data)
```

## 1.2 二元正态分布

```{r}
# 导入程序包
library(MASS)
```

```{r fig.height=4, fig.width=4}
# 参数
n <- 100
mu <- c(0, 0)
rho <- 0.75
sigma1 <- 1
sigma2 <- 1
Sigma <- as.matrix(rbind(c(sigma1, rho * sigma1 * sigma2), c(rho * sigma1 * sigma2, sigma2)), nrow = 2, ncol = 2)

# 生成二元正态分布
x <- mvrnorm(n = n, mu = mu, Sigma = Sigma)

# 绘制散点图
plot(x)
```

## 1.3 数据基本描述性统计

```{r}
# 将数据转化为矩阵
X <- as.matrix(data[1:4])
X
```

```{r}
# 均值向量
mu <- apply(X, 2, mean)
mu
```

```{r}
# X自身协方差阵
Sigma <- cov(X)
Sigma
```

```{r}
# 相关阵
R <- cor(X)
R
```

# 2 多元正态总体参数的假设检验

## **2.1 多元正态分布的均值向量的检验**

<https://blog.csdn.net/cassie_note/article/details/80210070>

### **2.1.1 单个总体协方差阵已知时均值向量的检验**

```{r}
mu.test.known <- function(data, mu0, Sigma0, alpha=0.05)   
###########################################################
## H0: mu=mu0 when Sigma0 is known
## this is a Chisq testing
##############  Input  ####################################
## data = design matrix with the ith sample in the ith line
## mu0 = mu0 for null hypothesis
## Sigma0 = the known variance matrix
## alpha = the significant level, default value = 0.05
############## Output  ####################################
## Reject.area = reject region
## p.value     = p value
###########################################################

{
  data <- as.matrix(data) #将数据框转化为矩阵#
  n <- nrow(data) #n行#
  p <- ncol(data) #p列#

  X.bar <- apply(data, 2, mean) #按列求均值#
  T1 <- n * t(X.bar - mu0) %*% solve(Sigma0) %*% (X.bar - mu0)
  
  #求下侧分位点，上侧：lower.tail = FALSE#
  a2 <- qchisq(1 - alpha, p) 
  
  #按行排#
  reject <- matrix(c(T1, a2), nrow=1) 
  rownames(reject) <- c("Reject") #行名#
  colnames(reject) <- c("Obs", ">1-alpha") #列名#
  
  #右半累积概率，T越大，P越小，越拒绝#
  pv <- 1-pchisq(T1, p)
  return(list(Reject.area = reject, p.value = pv))
}
```

```{r}
mu.test.known(X[10:20,], mu0 = mu, Sigma0 = Sigma, alpha=0.05)
```

### **2.1.2 单个总体协方差阵未知时均值向量的检验**

```{r}
mu.test <- function(data, mu0)   
############################################################
## H0: mu=mu0 when Sigma is unknown
## this is an F testing
##############  Input  #####################################
## data = design matrix with the ith sample in the ith line
## mu0 = mu0 for null hypothesis
############## Output  #####################################
## p.value = p value
############################################################
{
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)

  X.bar <- apply(data, 2, mean)
  A <- (n-1) * var(data)

  T2 <- (n - 1) * n * t(X.bar - mu0) %*% solve(A) %*% (X.bar - mu0)
  F <- (n-p)/((n-1)*p)*T2

  p.two <- 1-pf(F, p, n-p)
  return(list(p.value = p.two))
}
```

### 2.1.3 协**时**方差阵**未知但相等两个总体均值向量的检验**

```{r}
#### two independent normal distribution  ########################
two.mu.test <- function(data1, data2)   
##################################################################
## H0: mu1=mu2 
## this is an F testing
##############  Input  ###########################################
## data1 = design matrix for X with the ith sample in the ith line
## data2 = design matrix for X with the ith sample in the ith line
############## Output  ###########################################
## p.value = p value
##################################################################

{
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  p <- ncol(data1)

  X.bar <- apply(data1, 2, mean) 
  A1 <- (n1-1)*var(data1)
  Y.bar <- apply(data2, 2, mean)
  A2 <- (n2-1)*var(data2) 
  A <- (A1+A2)/(n1+n2-2)

  T2 <- (n1*n2/(n1+n2))*t(X.bar-Y.bar)%*%solve(A)%*%(X.bar-Y.bar)
  F <- (n1+n2-2-p+1)/((n1+n2-2)*p)*T2

  p.two <- 1-pf(F, p, (n1+n2-p-1))
  return(list(p.value=p.two))
}
```

### 2.1.4 **多个正态总体均值向量的检验：多元方差分析**

```{r}
#### k independent normal distribution ######### 
###################################################################
## H0: mu1=mu2=...=muk
## this is asymptotically a Chisq testing
##############  Input  ############################################
## data  = design matrix with a group index ind
############## Output  ############################################
## p.value = p value
###################################################################
multi.mu.test <- function(data, k)            

{
  ind <- data$ind

  n <- nrow(data)
  p <- ncol(data)-1

  data <- data[ ,1:p]
  T <- (n-1)*var(data)

  A <- 0
  for (i in 1:k)                                
  {
    datai <- data[ind == i, ]
    ni <- nrow(datai)                                 
    A <- A+(ni-1)*var(datai)
  }

  Lambda <- det(A)/det(T)
  n1 <- n-k
  n2 <- k-1
  r <- n1-(p-n2+1)/2
  Chi <- (-1)*r*log(Lambda)

  p.value <- 1-pchisq(Chi, p*n2)
  return(p.value = p.value)
}
```

## 2.2 **多元正态分布协方差阵的检验**

<https://blog.csdn.net/cassie_note/article/details/80212572>

### 2.2.1 单个p元正态总体协方差阵的检验

```{r}
#### variance testing ######### 
var.test <- function(data, Sigma0)
###############################################################
## H0: Sigma=Sigma0
## this is aymptotically a Chisq testing
##############  Input  ########################################
## data  = design matrix with the ith sample in the ith line
## Sigma0= Simga0 for null hypothesis
############## Output  ########################################
## p.value     = p value
###############################################################
{
  n=nrow(data)
  p=ncol(data)

  A=(n-1)*var(data)
  S=A%*%solve(Sigma0)

  lambda=exp(sum(diag((-1)*S/2)))*(det(S))^(n/2)*(exp(1)/n)^(n*p/2)
  T5=-2*log(lambda)

  p.value=1-pchisq(T5, p*(p+1)/2)
  return(p.value)
}
```

### 2.2.2 两个p元正态总体协方差阵的检验

```{r}
#### k independent normal distribution  ######### 
multi.var.test <- function(data, k)
 ###################################################################
## H0: Sigma1=Sigma2=...=Sigmak
## this is asymptotically a Chisq testing
##############  Input  ############################################
## data  = design matrix with a group index ind
############## Output  ############################################
## p.value     = p value
###################################################################
{
  ind=data$ind

  n=nrow(data)
  p=ncol(data)-1
  data=data[ ,1:p]

  A=0
  for (i in 1:k)                                
  {
    datai=data[ind==i, ]
    ni=nrow(datai)                                 
    A=A+(ni-1)*var(datai)
  }

  det.A=0
  for (i in 1:k)                                
  {
    datai=data[ind==i, ]
    ni=nrow(datai)                                 
    det.A=det.A+(ni-1)*log(det(var(datai)))
  }

  M=(n-k)*log(det(A/(n-k)))-det.A
  d=(2*p^2+3*p-1)*(k+1)/(6*(p+1)*(n-k))
  f=p*(p+1)*(k-1)/2

  T6=(1-d)*M

  p.value=1-pchisq(T6, f)
  return(p.value=p.value)
}
```

## 2.3 两总体独立性检验

<https://blog.csdn.net/cassie_note/article/details/80212572>

```{r}
######### testing for independent  ######### 
multi.ind.test <- function(data, k)
###################################################################
## H0: sigma(ij)=0,i!=j
## this is asymptotically a Chisq testing
##############  Input  ############################################
## data  = design matrix with a group index ind
############## Output  ############################################
## p.value     = p value
###################################################################
{
  ind=data$ind #不同总体分组#

  n=nrow(data)-1 #因为最后一行为ind#
  p=ncol(data)
  data=data[1:n , ]
  X.bar=apply(data, 2, mean)
  S=(n-1)*var(data)

  det.B=1
  P1=0
  P2=0
  P3=0
  for (i in 1:k)                                
  {
    datai=data[ ,ind==i]
    p_i=ncol(datai)
    P1=P1+p_i**3
    P2=P2+p_i**2
    P3=P3+p_i*(p_i+1)
    Si=(n-1)*var(datai)
    det.B=det.Si*det.B
  }

  V=det.S/det.B
  b=n-3/2-(p^3-P1)/(3*p^2-3*P2)
  f=0.5*(p(p+1)-P3)

  T6=(-b)*log(V)

  p.value=1-pchisq(T6, f)
  return(p.value=p.value)
}
```

```{r}
source("D:/R/TEST/MV_function/mulit_test.R")
```

## 2.4 多元正态性检验

```{r}
# 检验包
library(mvnormtest) # 基于Shapiro-Wilk
library(MVN) # 有多种方法
```

```{r}
# mvnormtest
data(EuStockMarkets) # 欧洲股票市场数据
C <- t(EuStockMarkets[15:29,1:4])
mshapiro.test(C)
```

```{r}
# MVN
mvn(data = X)
```
