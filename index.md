---
title       : Disease Propagation
subtitle    : using CIM model
author      : Ehsan Siavashi
job         : 
framework   : io2012        # {io2012, html5slides, shower, dzslides, ...}
highlighter : highlight.js  # {highlight.js, prettify, highlight}
hitheme     : tomorrow      # 
widgets     : []            # {mathjax, quiz, bootstrap}
mode        : selfcontained # {standalone, draft}
knit        : slidify::knit2slides
---

## What the app does?
This is an application which models disease propagation is a network of 31 people using the Conditional Influence Model. It has the following properties:

- The disease is contagous.
- People influence each other due to physical contact.
- The disease is fatal.
- There is immunization. If someone gets imminization, he/she will not get the disease again.

The following code shows the 

--- .class #id 
## How to find who gets the disease?
I have applied the following code for calculating the probabilities for getting sick in the next step:

```r
C <- matrix(0, N, N)
for (i in 1:N) {
        for (j in 1:N) {
                if(i==j){
                        C[i,j] <- 1
                }else
                        # Node i gets influenced by node j iff i is H and j is I.
                        C[i,j] <- S[3*(i-1)+1]*S[3*(j-1)+2]
        }
}
# transpose of E (E_t): Refer to the paper for the formula.
E_t <- t(D) * t(C) + diag(1,N,N) * (D %*% (matrix(1,N,N)-t(C)))

# Calculating A_0  
A_0 <- cbind(MC2,matrix(0,3,3*N-3))
O <- matrix(0,3,3)
for (i in 2:N) {
        A <- O
        for (j in 2:N) {
                if(i==j && (i %in% Inf.list)){
                        A<-cbind(A, E_t[i,i]*MC2)
                }
                else if(i==j){
                        A<-cbind(A, E_t[i,i]*MC1) 
                }else{
                        A<-cbind(A,O)
                }
        }
        A_0 <- rbind(A_0,A)
}

H <- A_0 + kronecker(E_t - (diag(1, N,N)*E_t), input$DI)

P <- S %*% H
```


--- .class #id 

## The input/output

Users can input the number of people infected and also the amount of influences between nodes by changing the Markov chain matrices. Select the parameter and the app will provide a table and a graph of the number of infected, healthy and dead prople for a period of 30 time steps. Copy right researved for the author.

# Thank you for your attention!
