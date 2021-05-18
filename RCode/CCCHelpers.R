library(DescTools)

getCCC = function(x, y) {
  CCC(x,y)$rho.c$est
}

getMSE = function(x, y) {
  mean((x - y)^2)
}

#test: 
#getMSE(c(1,1.5,1), c(1,1,1)) #0.833333. ok