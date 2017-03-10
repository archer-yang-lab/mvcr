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

