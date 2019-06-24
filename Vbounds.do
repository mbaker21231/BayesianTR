/* This function calculates bounds on the for each parameter, given a */
/* vector of ancilliary parameters */
/* Order v11 v22 v12 v1t v0t, say  */


/* sample guess as a parameter */
/* Actually that won't quite work */
/* We need to keep a running value of Sigma, and then from this */
/* work - returns upper and lower bounds on possible values */

real matrix newSigma (real scalar realno, real entry, real matrix v) {
	if (entry == 1) {
	    lb = v[1,3]^2*v[2,2] + v[1,2]^2 - 2*v[1,2]*v[1,3]*v[2,3]
		lb = lb / (v[2,2] - v[2,3]^2)
		nv = lb*(1 + exp(realno))
	}
	if (entry == 2) {
	    lb = v[2,3]^2*v[1,1] + v[1,2]^2 - 2*v[1,2]*v[1,3]*v[2,3]
		lb = lb / (v[1,1] - v[2,3]^2)	    
		return(lb*(1 + exp(realno)))
	if (entry == 3) {
	    lb = v[1,3]*v[2,3] - sqrt(v[1,3]^2*v[2,3]^2 - v[1,1]*v[2,3]^2-
		                          v[1,3]^2*v[2,2]   + v[1,1]*v[2,2])
		ub = v[1,3]*v[2,3] + sqrt(v[1,3]^2*v[2,3]^2 - v[1,1]*v[2,3]^2-
		                          v[1,3]^2*v[2,2]   + v[1,1]*v[2,2])
		return(lb*(1/(1 + exp(realno)) + exp(realno)/(1 + exp(realno)))
	if (entry == 4) {
	    lb = v[1,2]*v[2,3] - sqrt(v[1,2]^2*v[2,3]^2 + v[1,1]*v[2,2]^2-
		                          v[1,2]^2*v[2,2]   - v[1,1]*v[2,2]*v[2,3]^2)
	    ub = v[1,2]*v[2,3] - sqrt(v[1,2]^2*v[2,3]^2 + v[1,1]*v[2,2]^2-
		                          v[1,2]^2*v[2,2]   - v[1,1]*v[2,2]*v[2,3]^2)
	    return(lb*(1/1 + exp(realno)) + exp(realno)/(1 + exp(realno)))
	if (entry == 5) {
	    lb = v[1,2]*v[1,3] - sqrt(v[1,2]^2*v[1,3]^2 + v[1,1]^2*v[2,2] -
		                          v[1,1]^2*v[2,2]   - v[1,1]*v[2,2]*v[1,3]^2)
	    ub = v[1,2]*v[1,3] - sqrt(v[1,2]^2*v[1,3]^2 + v[1,1]^2*v[2,2] -
		                          v[1,1]^2*v[2,2]   - v[1,1]*v[2,2]*v[1,3]^2)
	    return(lb*(1/1 + exp(realno)) + exp(realno)/(1 + exp(realno)))	
}
