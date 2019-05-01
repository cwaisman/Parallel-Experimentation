library(np)

set.seed(290419)

rm(list = ls())

### Number of observations, advertisers and simulations

j = 1 # advertiser of interest (whose ATEs we want to recover)
J = 3 # number of advertisers
S = 100 # number of simulation rounds
n = 120 # sample size

### Possible treatment assignments

grid_tot = as.matrix(expand.grid(data.frame(matrix(rep((0:1),J),nrow=2)))) # 2^J including j
grid = unique(grid_tot[,-j]) # 2^(J-1) from j's perspective
G1 = length(grid_tot[,1])
G2 = length(grid[,1])

### Underlying parameters

# Each coordinate of theta corresponds to the expected value of the dependent variable
# conditional on the treatment assignment from the corresponding row in grid_tot.

theta = runif(2^J)

### Implied ATEs

# Fix a row in grid. The ATE for j at this row consists of the difference between theta
# at the row where the j-th coordinate equals 1 and the remaining coordinates correspond
# to the ones in grid and the theta at the row where the j-th coordinate equals 1 and 
# the remaining coordinates correspond to the ones in grid.

ATEs = matrix(rep(NA,G2))
for (g in 1:G2) {

  ATEs[g] = theta[rowSums(1*(grid_tot[,-j]==kronecker(matrix(rep(1,G1),nrow=G1),t(grid[g,]))))==J-1 & grid_tot[,j]==1] - 
    theta[rowSums(1*(grid_tot[,-j]==kronecker(matrix(rep(1,G1),nrow=G1),t(grid[g,]))))==J-1 & grid_tot[,j]==0]

}
ATE = mean(ATEs)

### Start running simulations

hat_ATEs = matrix(NA,S,2^(J-1))
hat_ATEsk = matrix(NA,S,2^(J-1))

for (s in 1:S) {
  
  ### Generating equiprobable treatment assignments
  
  # D contains the treatment assignments for the n observations in simulation s.
  # p contains the corresponding expected value of the dependent variable given
  # the treatment assignment in D.
  # The loop over g pre-multiplies all assignments by g/G1 so that they all happen
  # with the same probability.
  
  aux = runif(n)
  for (g in 1:G1){
    
    if (g==1) {
      
      aux2 = g*(aux<g/G1)
      D = (aux2==g)*matrix(rep(grid_tot[g,],each=n),nrow=n)
      p = (aux2==g)*theta[g]
      
    } else if (g>1 & g<G1) {
      
      aux2 = aux2 + g*(aux>=((g-1)/G1))*(aux<(g/G1))
      D = D + (aux2==g)*matrix(rep(grid_tot[g,],each=n),nrow=n)
      p = p + (aux2==g)*theta[g]
      
    } else {
      
      aux2 = aux2 + g*(aux>=((g-1)/G1))
      D = D + (aux2==g)*matrix(rep(grid_tot[g,],each=n),nrow=n)
      p = p + (aux2==g)*theta[g]
      
    }
    
  }
  colnames(D) = c(sprintf("D%02d", seq(1,J)))
  
  ### Generating outcomes (visit (1) vs not visit (0))

  Y = 1*(runif(n)<p)

  ### Estimated ATEs for advertiser of interest
  
  ## Unconditional (just a regression of Y on a constant and D[,j])
  
  X = cbind(rep(1,n),D[,j])
  hat = solve(t(X)%*%X)%*%(t(X)%*%Y)
  hat_ATE = hat[2]
  
  ## Conditional ATEs via linear regression
  
  # The matrix Xt gathers interactions between the matrix X and indicators for 
  # all the 2^(J-1) possible treatment assignments from j's perspective.
  
  for (g in 1:G2){
    
    if (g==1){
      
      Xt = (rowSums(D[,-j]==kronecker(matrix(rep(1,n),nrow=n),t(grid[g,])))==J-1)*X
      
    } else {
      
      Xt = cbind(Xt,(rowSums(D[,-j]==kronecker(matrix(rep(1,n),nrow=n),t(grid[g,])))==J-1)*X)
      
    }
    
  }
  
  hats = solve(t(Xt)%*%Xt)%*%(t(Xt)%*%Y)
  hat_ATEs[s,] = hats[seq(2,length(hats),by=2)]
  
  ## Conditional ATEs via kernel-based method
  
  # The matrix Z contains the 2^(J-1) possible treatment assignments from j's 
  # perspective. They are converted to factors so that the np package treats
  # them as categorical variables instead of continuous.
  
  Z = data.frame(D[,-j])
  
  for (z in colnames(Z)) {
    
    Z[[z]] = as.factor(Z[[z]])
    
  }

  # npscoefbw computes the bandwidths via leave-one-out cross-validation as recommended by 
  # Li et al (2013). We specify the kernel type to be 'liracine', which corresponds to the 
  # kernel given in Li et al (2013).

  sets<-npscoefbw(xdat=D[,-j],ydat=Y,zdat=Z,ukertype='liracine',optim.method='CG')
  summary(sets)
  lambdas = t(matrix(unlist(sets['bw'])))

  Lb = matrix(rep(lambdas,n),nrow=n,byrow=TRUE)

  # We then compute the estimator given by Li et al (2013) at each of the 2^(J-1)
  # treatment assignments from j's perspective.
  
  for (i in 1:G2) {

    print(i)
    auxil1 = Lb^(Z!=matrix(unlist(rep(grid[i,],n,nrow=n,byrow=TRUE))))
    auxil2 = apply(auxil1,1,prod)
    bes = solve(t(X*cbind(auxil2,auxil2))%*%X)%*%(t(X*cbind(auxil2,auxil2))%*%Y)
    hat_ATEsk[s,i] = bes[2]

  }

}
