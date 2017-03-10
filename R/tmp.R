require(MASS)
require(glmnet)
require(pcaPP)
require(glasso)
require(corpcor)
require(gglasso)
require(expm)
require(ramify)
require(matrixStats)


library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)

options(digits=3)

wd <- "C:/Users/newtonii/Documents/Yue writeup/MATLAB/MVR"
setwd(wd)
write_filename <- paste(wd, "/", "R_running_test.txt",sep = "")



#08: use the approximate MRCE method

rm(list = ls())


#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------


#Compute the Kendall's tau matrix
#The input data_matrix is of size n by p, with each row of 1 by p vector for a single observation
#Output is p by p matirx of Kendall's tau (before the sine function transform) 
#Calculation is done element-wise, with fast Kendall's tau implementation in  package pcaPP
compute_Kendall_matrix <- function(data_matrix) {
	p <- ncol(data_matrix)
	Kendall_matrix <- mat.or.vec(p,p)
	for (k in 2:p){
		for (l in 1:(k-1)){
			Kendall_matrix[k,l] = cor.fk( data_matrix[,k], data_matrix[,l])
		}
	}
	Kendall_matrix <- Kendall_matrix + 0.5 * diag(p)
	Kendall_matrix <- Kendall_matrix + t(Kendall_matrix)
	return(Kendall_matrix)
}


#--------------------------------------------------------------------------------

#From given design matrix S_hat_xx = S_hat[1:p,1:p], cross-correlation matrix S_hat_xy = S_hat[1:p,(p+1):(p+q)],
#and precision matrix Omega, returns equivalent matrices X_c and Y_c
#such that X_c^T X_c = kronecker( Omega, S_hat_xx ) and Y_c^T X_c = gamma (needs more explanation)
return_equivalent_X_Y <- function(S_hat, Omega) {
S_hat_xx <- S_hat[1:p,1:p]
S_hat_xy <- S_hat[1:p,(p+1):(p+q)]
XTX <- kronecker( Omega, S_hat_xx )
Y_equiv <- cbind(rep(0,p*q))
SxOmega <- cbind(c( S_hat_xy %*% Omega ))
X_equiv <- chol(XTX)
Y_equiv <- solve( t(X_equiv), SxOmega)
#print(Y_equiv)

##Legacy code for computing Y_equiv
##----------------------------------------
#XTY <- mat.or.vec(p*q,1)
#for (l in 1:q){
#	Sigma_diag <- diag(Omega[l,])
#	SdKxy <- S_hat_xy %*% Sigma_diag
#	XTY[((l-1)*p+1):(l*p)] <- rowSums (SdKxy, na.rm = FALSE)
#}
#Y_equiv_legacy <- solve( t(X_equiv), XTY)
##----------------------------------------

#newList <- list("X_c" = X_equiv, "Y_c" = Y_equiv, "Y_c_legacy" = Y_equiv_legacy)
newList <- list("X_c" = X_equiv, "Y_c" = Y_equiv, "SxOmega" = SxOmega)
return(newList)
}



#--------------------------------------------------------------------------------

### (Yue:) This is the algorithm from the CoCoLasso paper
### admm algo to create the positive definite estimate of the gram matrix ###
maxproj.cov<-function(mat, epsilon=1e-4, mu=10, nitr.max=1e3, etol=1e-4){

	p<-nrow(mat)
	
	# Initialization
	R<-diag(mat)
	S<-matrix(0,p,p)
	L<-matrix(0,p,p)
	
	itr<-0
	while (itr<nitr.max) {
	Rp<-R
	Sp<-S
	
	# Subproblem I: R step
	W<-mat+S+mu*L
	W.eigdec<-eigen(W, symmetric=TRUE)	
	W.V<-W.eigdec$vectors
	W.D<-W.eigdec$values
	R<-W.V%*%diag(pmax(W.D,epsilon))%*%t(W.V)
	
	# Subproblem II: S step
	M<-R-mat-mu*L	
	S[lower.tri(S, diag = TRUE)]<-M[lower.tri(M, diag = TRUE)]-l1proj(v=M[lower.tri(M, diag = TRUE)],b=mu/2)	
	for (i in 2:p){
		for (j in 1:(i-1)){
			S[j,i]<-S[i,j]
		}
	}
	
	# L step: update the Lagrange parameter
	L<-L-(R-S-mat)/mu
	
    # Stopping Rule                        
	#cat("check the stopping criterion:",max(abs(R-S-mat)),"\n")
    if ((max(abs(R-Rp))<etol) && (max(abs(S-Sp))<etol) && (max(abs(R-S-mat))<etol)){
	  itr<-nitr.max
    } else {
		itr<-itr+1
	}
	
	if (itr%%20==0) {
		mu<-mu/2
	}
	}
	
	return(R)

}

# Efficient projection onto L1 ball of specified radius (i.e. b), used by the admm algo
# Ref. Duchi et al. (2008). Efficient Projections onto the L1-Ball for Learning in High Dimensions, ICML
l1proj<-function(v, b){

	stopifnot(b>0)
	
	u <- sort(abs(v),decreasing=TRUE)
	sv <- cumsum(u)
	rho <- max(which(u>(sv-b)/1:length(u)))
	theta <- max(0, (sv[rho]-b)/rho)
	w <-sign(v) * pmax(abs(v)-theta,0)

	return(w)
}



#--------------------------------------------------------------------------------
#This tries to solve for B under a given Omega
#The implementation is Algorithm 1 in Adam et al (page 951)
#But I don't think I made it to work properly
#(Compute B adjusted for correlation)
compute_B <- function(S_mat, S_xy_mat, Omega_mat, lambda, epsilon = 10^(-8)) {
	H_mat <- S_xy_mat %*% Omega_mat
	B_mat_ridge <- solve(S_mat + lambda * diag(p)) %*% S_xy_mat
	B_ridge_sum <- sum( abs(B_mat_ridge) )
	B_mat_old <- matrix(0,nrow(S_xy_mat),ncol(S_xy_mat))
	B_mat_new <- B_mat_old

	repeat {
		U_mat <- S_mat %*% B_mat_old %*% Omega_mat
		for (r in 1:p) {
			for (c in 1:q) {
				im_prod <- S_mat[r,r] * Omega_mat[c,c]
				im <- B_mat_old[r,c] + (H_mat[r,c]-U_mat[r,c])/im_prod
				B_mat_new[r,c] <- sign( im ) * 
					max(0, abs( im ) - lambda / im_prod )
			}
		}
		#print(sum( abs(B_mat_old-B_mat_new) ))
		if ( sum( abs(B_mat_old-B_mat_new) ) > 10^(10) ) break
		if ( sum( abs(B_mat_old-B_mat_new) ) < epsilon * B_ridge_sum ) break
		B_mat_old <- B_mat_new
	}
	return(B_mat_new)
}

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------




#Some initial setup
N <- 500
N_fold <- 5
fold_size <- N/N_fold
cv_starts <- seq(1,N,fold_size)
cv_ends <- seq(fold_size,N,fold_size)


p <- 10
q <- 20

#This generates an initial coefficient matrix B
repeat {
B <- matrix(0, nrow=p, ncol=q)
mu <- 0
W <- matrix( mvrnorm(n = p*q, mu, 1, tol = 1e-6, empirical = FALSE, EISPACK = FALSE) , nrow=p, ncol=q )
s_1 <- 0.1
K <- matrix( rbinom(n = p*q, 1, s_1) , nrow=p, ncol=q )
s_2 <- 1
Q <- matrix( rep( rbinom(n = p, 1, s_2), q), nrow=p, ncol=q )
#Q_row <- rep(0,p) #The following three lines are for group Lasso, to set a few rows of B to be all zero; currently disabled
#Q_row[1:(floor(p/1))] <- 1
#Q <- matrix( rep( Q_row, q), nrow=p, ncol=q )
B <- W * K * Q
if ( sum( abs( B ) ) != 0 ) break
} #This checks if B happens to be the zero matrix, which is not allowed
#print(B)


B_hat_1 <- mat.or.vec(p,q)


Sigma_e_diag = 0.5; #Set the largest value of the matrix Sigma_{\bvare\bvare}
Sigma_xx <- toeplitz(0.7 ^ seq(0,p-1) ) #Generates the population design matrix
max_element <- max(diag( t(B) %*% Sigma_xx %*% B))
B <- B/sqrt(max_element)*sqrt(1-Sigma_e_diag) #Adjust the coefficient matrix B so that the matrix Sigma_yy indeed has ones on the diagonal

rho_e <- 0.8;
Sigma_ee <- Sigma_e_diag * toeplitz( rho_e ^ seq(0,q-1) ); #Generates Sigma_{\bvare\bvare}
Omega_ee <- solve(Sigma_ee)

Sigma_yy <- t(B) %*% Sigma_xx %*% B + Sigma_ee

Sigma_xz <- mat.or.vec(p+q,p+q)
Sigma_xz[1:p,1:p] <- Sigma_xx
Sigma_xz[(p+1):(p+q),(p+1):(p+q)] <- Sigma_yy
Sigma_xz[1:p,(p+1):(p+q)] <- Sigma_xx %*% B
Sigma_xz[(p+1):(p+q),1:p] <- t(B) %*% Sigma_xx
Sigma_xy <- Sigma_xx %*% B

data_matrix <- mvrnorm(n = N, rep(0,p+q), Sigma_xz, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

T_xz_hat <- mat.or.vec(p+q,p+q)
S_xz_hat <- mat.or.vec(p+q,p+q)

T_xz_hat <- compute_Kendall_matrix(data_matrix)
S_xz_hat <- sin(pi/2*T_xz_hat)
S_xx_hat <- S_xz_hat[1:p,1:p]
S_xy_hat <- S_xz_hat[1:p,(p+1):(p+q)]
S_yy_hat <- S_xz_hat[(p+1):(p+q),(p+1):(p+q)]



if(min(eigen(S_xx_hat)$value) < 1e-3) {
	S_xx_hat_p <- maxproj.cov(mat=S_xx_hat, epsilon=1e-3)
} else {
	S_xx_hat_p <- S_xx_hat
}



data_matrix_cv <- mat.or.vec(fold_size*(N_fold-1),p)
data_matrix_validation <- mat.or.vec(fold_size,p)

#--------------------------------------------------------------------------------
#Compute cv matrices (used to estimate parameters) and validation matrices (used to evaluate estimated parameters).  We have N_folds number of each matrices.  For a single index cv_index, the samples used to generate the cv matrix and the samples used to generate the validation matrix form a partition of the entire data sample

Kendall_matrix_cv_all <- array(0,dim=c(p+q,p+q,N_fold))
Kendall_matrix_validation_all <- array(0,dim=c(p+q,p+q,N_fold))

for (cv_index in 1:N_fold){
	sample_indices_cv <- ( (seq(cv_starts[cv_index]+fold_size,cv_starts[cv_index]+fold_size*(N_fold-1)-1,1)-1) %% N)+1
	data_matrix_cv <- data_matrix[sample_indices_cv,]
	Kendall_matrix_cv_all[,,cv_index] <- compute_Kendall_matrix(data_matrix_cv)

	sample_indices_validation <- seq(cv_starts[cv_index],cv_starts[cv_index]+fold_size-1,1)
	data_matrix_validation <- data_matrix[sample_indices_validation,]
	Kendall_matrix_validation_all[,,cv_index] <- compute_Kendall_matrix(data_matrix_validation)
}


#--------------------------------------------------------------------------------
#Get preliminary estimate of B


#lambda_init = 16*sqrt(log(p+q)/N) * (2* norm(B,type = c("O")) + 1)
#N_lambda <- 100
#lambda_sq <- lambda_init * sqrt(2)^(N_lambda/10-1:N_lambda)
#lambda_sq_all <- matrix( t(rep(lambda_sq,q)), nrow=N_lambda, ncol=q)

N_lambda = 20
lambda_sq_all <- matrix(0, nrow=N_lambda, ncol=q)
beta_sq_all <- array(0,dim=c(N_lambda,p,q))

#Produce equivalent matrix X (covariates) and Y (response) from Kendall's tau matrices so that we can use the standard LASSO estimator
#I guess I could have called the function return_equivalent_X_Y to do this
X_c <- chol(S_xx_hat_p)
Y_c <- solve( t(X_c), S_xy_hat)

for (k in 1:q){
	prelim_fitinfo <- glmnet(X_c, Y_c[,k],intercept=FALSE,nlambda=N_lambda)
	n_lambda <- length(prelim_fitinfo$lambda)
	lambda_sq_all[1:n_lambda,k] <- prelim_fitinfo$lambda
	for (lambda_index in 1:n_lambda){
		beta_sq_all[lambda_index,,k] <- prelim_fitinfo$beta[,lambda_index]
	}
}

cv_loss_all <- mat.or.vec(N_lambda,q)

for (cv_index in 1:N_fold){

	T_xz_hat_cv <- Kendall_matrix_cv_all[,,cv_index]
	S_xz_hat_cv <- sin(pi/2*T_xz_hat_cv)
	S_xx_hat_cv <- S_xz_hat_cv[1:p,1:p]
	S_xy_hat_cv <- S_xz_hat_cv[1:p,(p+1):(p+q)]

	if(min(eigen(S_xx_hat_cv)$value) < 1e-3) {
		S_xx_hat_p_cv <- maxproj.cov(mat=S_xx_hat_cv, epsilon=1e-3)
	} else {
		S_xx_hat_p_cv <- S_xx_hat_cv
	}


	X_c_cv <- chol(S_xx_hat_p_cv)
	Y_c_cv <- solve( t(X_c_cv), S_xy_hat_cv)

	sample_indices_validation <- seq(cv_starts[cv_index],cv_starts[cv_index]+fold_size-1,1)
	data_matrix_validation <- data_matrix[sample_indices_validation,]

	T_xz_hat_validation <- Kendall_matrix_validation_all[,,cv_index]
	S_xz_hat_validation <- sin(pi/2*T_xz_hat_validation)
	S_xx_hat_validation <- S_xz_hat_validation[1:p,1:p]
	S_xy_hat_validation <- S_xz_hat_validation[1:p,(p+1):(p+q)]

	if(min(eigen(S_xx_hat_validation)$value) < 1e-3) {
		S_xx_hat_p_validation <- maxproj.cov(mat=S_xx_hat_validation, epsilon=1e-3)
	} else {
		S_xx_hat_p_validation <- S_xx_hat_validation
	}


	for (k in 1:q){
		lambda_index <- which(lambda_sq_all[,k] > 0)
		lambda_sq <- lambda_sq_all[lambda_index,k]
		n_lambda <- length(lambda_sq)
		fitinfo_cv <- glmnet(X_c_cv, Y_c_cv[,k],intercept=FALSE,lambda=lambda_sq)
		
		for (lambda_index in 1:n_lambda) {
			beta <- fitinfo_cv$beta[,lambda_index]
			cv_loss_all[lambda_index,k] <- cv_loss_all[lambda_index,k] + t(beta) %*% S_xx_hat_p_validation %*% beta - 2 * t( S_xy_hat_validation[,k] )  %*% beta
		}

	}
}


for (k in 1:q){
	min_index <- which.min( cv_loss_all[,k] )
	B_hat_1[,k] <- beta_sq_all[min_index,,k]
}



#--------------------------------------------------------------------------------
#Compute S_ee_hat

S_ee_hat <- S_yy_hat - t(B_hat_1) %*% S_xx_hat %*% B_hat_1
if(min(eigen(S_ee_hat)$value) < 1e-3) {
	S_ee_hat_p <- maxproj.cov(mat=S_ee_hat, epsilon=1e-3)
} else {
	S_ee_hat_p <- S_ee_hat
}


#--------------------------------------------------------------------------------
#Parameter sequence for estimating Omega_ee
#Compute Omega_ee_hat

n_rholist <- 16
rho_init <- sqrt( ( log(p^2) + log(q^2) ) / N)
rho_sq <- rho_init * ( 4^( ( seq(0,n_rholist-1,1) - n_rholist / 2 ) ) )
glasso_output <- glassopath(S_ee_hat_p,rholist = rho_sq)



#--------------------------------------------------------------------------------
#Parameter sequence for ordinary (that is, "_ord") and Group Lasso (that is, "_grop")

N_lambda <- 32
lambda_init_group <- 1
lambda_sq_group <- lambda_init_group * 2 ^( 0.5 * ( seq(N_lambda-1,0,-1) - N_lambda/2 ) ) * sqrt(q/N)
lambda_sq_ord <- lambda_sq_group / sqrt(q)


cv_loss_ord_all_folds <- array(0,dim=c(N_lambda,n_rholist,N_fold))



#--------------------------------------------------------------------------------
#Iteration over folds for cross validation

for (cv_index in 1:N_fold) {

	T_xz_hat_cv <- Kendall_matrix_cv_all[,,cv_index]
	S_xz_hat_cv <- sin(pi/2*T_xz_hat_cv)
	S_xx_hat_cv <- S_xz_hat_cv[1:p,1:p]
	S_xy_hat_cv <- S_xz_hat_cv[1:p,(p+1):(p+q)]
	S_yy_hat_cv <- S_xz_hat_cv[(p+1):(p+q),(p+1):(p+q)]

	S_ee_hat_cv <- S_yy_hat_cv - t(B_hat_1) %*% S_xx_hat_cv %*% B_hat_1
	if(min(eigen(S_ee_hat_cv)$value) < 1e-3) {
		S_ee_hat_p_cv <- maxproj.cov(mat=S_ee_hat_cv, epsilon=1e-3)
	} else {
		S_ee_hat_p_cv <- S_ee_hat_cv
	}

	glasso_output_cv <- glassopath(S_ee_hat_p_cv, rholist = rho_sq)

	if(min(eigen(S_xx_hat_cv)$value) < 1e-3) {
		S_xx_hat_p_cv <- maxproj.cov(mat=S_xx_hat_cv, epsilon=1e-3)
	} else {
		S_xx_hat_p_cv <- S_xx_hat_cv
	}

	S_xz_hat_cv[1:p,1:p] <- S_xx_hat_p_cv


	T_xz_hat_validation <- Kendall_matrix_validation_all[,,cv_index]
	S_xz_hat_validation <- sin(pi/2*T_xz_hat_validation)
	S_xx_hat_validation <- S_xz_hat_validation[1:p,1:p]
	S_xy_hat_validation <- S_xz_hat_validation[1:p,(p+1):(p+q)]


	if(min(eigen(S_xx_hat_validation)$value) < 1e-3) {
		S_xx_hat_p_validation <- maxproj.cov(mat=S_xx_hat_validation, epsilon=1e-3)
	} else {
		S_xx_hat_p_validation <- S_xx_hat_validation
	}


	#--------------------------------------------------------------------------------
	#(inner) Loop over rholist, for ordinary Lasso

	cv_loss_ord_all <- matrix(0, nrow = N_lambda, ncol = n_rholist)
	group_structure <- 1:(p*q)

	cv_loss_ord_all <- foreach (rho_index=1:n_rholist, .packages=c('expm', 'glmnet', 'gglasso', 'ramify') ,  .combine='cbind' ) %dopar% {
	#for (rho_index in 1:n_rholist) {

		Omega_ee_hat_cv <- glasso_output_cv$wi[,,rho_index]
		newList <- return_equivalent_X_Y( S_xz_hat_cv, Omega_ee_hat_cv )
		X_c_cv <- newList$X_c
		Y_c_cv <- newList$Y_c

		#gfitinfo_cv <- gglasso(X_c_cv, Y_c_cv, group = group_structure, lambda=lambda_sq_ord )
		fitinfo_ord_cv <- glmnet(X_c_cv, Y_c_cv, intercept=FALSE, lambda=lambda_sq_ord )
		
		cv_loss <- matrix(0, nrow = N_lambda, ncol = 1)

		for (lambda_index in 1:N_lambda) {
			#beta <- gfitinfo_cv$beta[,lambda_index]
			beta <- fitinfo_ord_cv$beta[,lambda_index]
			B_cv <- matrix(beta, nrow=p, ncol=q)
			#print( lambda_index )
			#print( B_cv )
			cv_loss[lambda_index] <- tr( t(B_cv) %*% S_xx_hat_p_validation %*% B_cv %*% Omega_ee_hat_cv -
				2 * t( S_xy_hat_validation )  %*% B_cv  %*% Omega_ee_hat_cv )
		}
		#cv_loss_ord_all[,rho_index] = cv_loss
		cv_loss
	}
	cv_loss_ord_all_folds[,,cv_index] <- cv_loss_ord_all
}


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------




cv_loss_ord_all_folds_median <- apply(cv_loss_ord_all_folds, c(1,2), "median")
min_index_ord <- which( cv_loss_ord_all_folds_median == min( cv_loss_ord_all_folds_median ), arr.ind = TRUE )
#print(min_index_ord)


#min_index_ord[2] is the index of the best tuning parameter for the graphical LASSO estimator
#Therefore the following line find the best graphical LASSO tuning parameter, and do a graphical LASSO estimator over the ENTIRE (that is, non-split) data
Omega_ee_hat_ord <- glasso_output$wi[,,min_index_ord[2]]
#We then use this estimator of Omega to do a final, second estimation of \bB^* using the ENTIRE data
newList <- return_equivalent_X_Y(S_xz_hat, Omega_ee_hat_ord)
X_c <- newList$X_c
Y_c <- newList$Y_c
group_structure <- rep(1:(p*q))
gfitinfo <- glmnet(X_c, Y_c, intercept=FALSE, lambda=lambda_sq_ord )
beta_ord <- gfitinfo$beta[,min_index_ord[1]]
#B_hat_ord is the final second estimation of \bB^*
B_hat_ord <- matrix( beta_ord, nrow = p, ncol = q )


#This line prints out the index of minimizer of the two-dimensional grid of the tuning parameters for Graphical LASSO and the second estimation of \bB^*
print(min_index_ord)


#The following prints out the prediction error and the error of the coefficient matrix of the second estimation of \bB^*, using the regular (i.e., non-group) LASSO method
dev_ord <- tr( t(B-B_hat_ord) %*% Sigma_xx %*% (B-B_hat_ord) )
print( norm(B-B_hat_ord), type = c("F") )
print( dev_ord )


#The following prints out the prediction error and the error of the coefficient matrix of the FIRST estimation of \bB^*, that is, the one down with a column-by-column estimation
dev_1 <- tr( t(B-B_hat_1) %*% Sigma_xx %*% (B-B_hat_1) )
print( norm(B-B_hat_1), type = c("F") )
print( dev_1 )


#print(B)
#print(B_hat_1)
#print(B_hat_ord)

#stopCluster(cl)
#showConnections(all = TRUE)
closeAllConnections()