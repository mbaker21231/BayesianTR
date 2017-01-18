/* In this worksheet, we will describe Bayesian 0,1 endogenous treatment effects */
/* This gives us personalized treatment effects and also gives us variances, etc. */
/* Highly parametric, but good in small samples. */

/* After having experimented with doing this with latent variables, we decided      */
/* it would work better without doing latent variables. This is because the process */
/* just never seemed to converge. In this version, we have five things to draw:     */
/* variance in outcome 1, variance in outcome 0, covariance between 1 and 0,        */
/* covariance between 1 and treatment, and covariance between 2 and treatment       */

clear all
use http://www.stata-press.com/data/r13/union3
set more off
keep if union != .
keep if tenure != . 
describe

regress wage age grade smsa black tenure union

etregress wage age grade smsa black tenure, treat(union = south black tenure)
ivregress 2sls wage age grade smsa black tenure (union = south black tenure)

mata:
	st_view(y=., ., "wage")
	st_view(tr=., ., "union")
	st_view(X=., ., "age grade smsa black tenure")
	st_view(W=., ., "south black tenure")
	X = X, J(rows(y), 1, 1)
	W = W, J(rows(y), 1, 1)
end

reg wage age grade smsa black tenure if union == 1
mat binit1 = e(b)

reg wage age grade smsa black tenure if union == 0
mat binit0 = e(b) 

probit union south black tenure
mat einit = e(b)

/* New likelihood to include more stuff...    */
/* We pass along six total parameters...      */
/* We also along two latent random effects    */
/* Everything else is handled before this     */
/* Prior_piv should contain four terms        */
/* T is an indicator as to whether or not     */
/* the treatment was received                 */
/* Note that we also add a "Flag" variable    */
/* which lets us know if we should return     */
/* the sum or the individual row elements     */
/* the latter will help drawing individual    */
/* level parameters                           */

mata: 
	b1 = st_matrix("binit1")
	b0 = st_matrix("binit0")
	e  = st_matrix("einit")
end

mata:
	lnsd1 = 0
	lnsd0 = 0
	v10 = 0
	v1t = 0
	v0t = 0	
end

/* A stable inverse-normal calculator */
mata:
	real matrix invnormstab(X) {
		XHat = editvalue(X, 0, 1e-323)
		XHat = editvalue(XHat, 1, 1e-16)
		return(invnormal(XHat))
	}
end
mata:
	real scalar ln_L(real matrix errs, real matrix Sigma) {
        part1 = -cols(errs)*rows(errs)/2*ln(2*pi())
	    part2 = -1/2*colsum(rowsum( (errs*invsym(Sigma):*errs)))
	    part3 = -rows(errs)/2*ln(det(Sigma))
        return(part1 + part2 + part3)
    }
end

/* Initial values for all the latent variables */

mata:
	y0Hat = X*b0' + rnormal(rows(y),1,0,1)*exp(lnsd0)
	y1Hat = X*b1' + rnormal(rows(y),1,0,1)*exp(lnsd1)
	y1 = tr:*y + (1 :- tr):*y1Hat
	y0 = (1 :- tr):*y + tr:*y0Hat
	
	muz = W*e'
	et  = invnormstab( normal(-muz) + (1 :- normal(-muz)):*runiform(rows(muz),1) )
	ent = invnormstab( normal(-muz):*runiform(rows(muz),1) )
	z = muz + et:*tr + ent:*(1 :- tr)
	
	m0 = X*b0'
	m1 = X*b1'
	mt = W*e'
	ey1 = (y1 - m1)
	ey0 = (y0 - m0)
	et  = (z - mt)	
	
end

/* Throughout the following does not change, and we use it alot, so we might */
/* as well just compute it now */

mata:
    XX = invsym(X'X)
	WW = invsym(W'W)
end

mata: 
	/* First, given parameters and everything else, we see that if y is missing */
	/* it is simply normal with variance vy1 */
	
	/* Vectors for collecting stuff */
	eHold = e
	b1Hold = b1
	b0Hold = b0
	sd0Hold = lnsd1
	sd1Hold = lnsd0
	v10Hold = v10
	v1tHold = v1t
	v0tHold = v0t
end	

/* Placeholders for proposal distributions and all that: */

mata:
    prosd1 = 1
	prosd0 = 1
	prov10 = 1
	prov1t = 1
	prov0t = 1
	
	delta =  1
	asta  = .4
end
	
mata:

for (i=1;i<=200000;i++) {

/* First part - regression parameters and latent values for y1 */
	
	Sig12 = (v10, v1t)
	Sig22 = exp(lnsd0)^2, v0t \ v0t, 1
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey0, et) )
	CV = exp(lnsd1)^2 - Sig12*Sig22m1*Sig12'
	
	mc1   = m1 + CM
	y1Hat = mc1 + rnormal(rows(y),1,0,1)*sqrt(CV)
	y1    = tr:*y + (1 :-tr):*y1Hat
	
	mb1 = XX*X'(y1 - CM)
	vb1 = exp(lnsd1)^2*XX
	b1 = mb1 + cholesky(vb1)*rnormal(cols(b1), 1, 0, 1)
	b1 = b1'	

	m1 = X*b1'
	ey1 = (y1 - m1)
	
/* Second part - regression parameters and latent values for y0 */

	Sig12 = (v10, v0t)
	Sig22 = exp(lnsd1)^2, v1t \ v1t, 1
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey1, et) )
	CV = exp(lnsd0)^2 - Sig12*Sig22m1*Sig12'
	
	mc0   = m0 + CM
	y0Hat = mc0 + rnormal(rows(y),1,0,1)*sqrt(CV)
	y0 = (1 :- tr):*y + tr:*y0Hat
	
	mb0 = XX*X'(y0 - CM)
	vb0 = exp(lnsd0)^2*XX
	b0 = mb0 + cholesky(vb0)*rnormal(cols(b0), 1, 0, 1)
	b0 = b0'	

	m0 = X*b0'
	ey0 = (y0 - m0)	
	
/* Third part - regression parameters and latent values for z */
	
	Sig12 = (v1t, v0t)
	Sig22 = exp(lnsd1)^2, v10 \ v10, exp(lnsd0)^2
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey1, ey0) )
	CV = 1 - Sig12*Sig22m1*Sig12'	
	
	mct = mt + CM
	et  = invnormstab( normal(-mct) + (1 :- normal(-mct)):*runiform(rows(mct),1) )
	ent = invnormstab( normal(-mct):*runiform(rows(mct),1) )
	z = mct + et:*tr + ent:*(1 :- tr)
	
	/* Draw eta coefficients given unit variance, common values, etc. */
	
	meane = WW*W'(z - CM)
	vare  = WW
	e = meane + cholesky(vare)*rnormal(cols(e), 1, 0, 1)
	e = e'
	
	mt = W*e'
	et = (z - mt)
	
	/* Next: variance parameters */

	gam = 1/i^delta
	
	/* lnsd1 */
	
	lnsd1Hat = lnsd1 + rnormal(1,1,0,1)*prosd1
	Sigma    = exp(lnsd1)^2,   v10, v1t \ v10, exp(lnsd0)^2, v0t \ v1t, v0t, 1
	SigmaHat = exp(lnsd1Hat)^2,v10, v1t \ v10, exp(lnsd0)^2, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) lnsd1 = lnsd1Hat
	    prosd1 = exp(gam*(alpha - asta))*prosd1
	}
	else {
	    prosd1 = exp(-asta*gam)*prosd1
	}
	
	/* lnsd0 */
	
	lnsd0Hat = lnsd0 + rnormal(1,1,0,1)*prosd0
	Sigma    = exp(lnsd1)^2, v10, v1t \ v10, exp(lnsd0)^2,    v0t \ v1t, v0t, 1
	SigmaHat = exp(lnsd1)^2, v10, v1t \ v10, exp(lnsd0Hat)^2, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) lnsd0 = lnsd0Hat
	    prosd0 = exp(gam*(alpha - asta))*prosd0
	}
	else {
	    prosd0 = exp(-asta*gam)*prosd0	
	}
	
	/* v10 */
	
	v10Hat = v10 + rnormal(1,1,0,1)*prov10
	Sigma    = exp(lnsd1)^2, v10,    v1t \ v10,    exp(lnsd0)^2, v0t \ v1t, v0t, 1
	SigmaHat = exp(lnsd1)^2, v10Hat, v1t \ v10Hat, exp(lnsd0)^2, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) v10 = v10Hat
	    prov10 = exp(gam*(alpha - asta))*prov10
	}
	else {
	    prov10 = exp(-asta*gam)*prov10
	}
	
	/* v1t */
	
	v1tHat = v1t + rnormal(1,1,0,1)*prov1t
	Sigma    = exp(lnsd1)^2, v10, v1t    \ v10, exp(lnsd0)^2, v0t \ v1t,    v0t, 1
	SigmaHat = exp(lnsd1)^2, v10, v1tHat \ v10, exp(lnsd0)^2, v0t \ v1tHat, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) v1t = v1tHat
	    prov1t = exp(gam*(alpha - asta))*prov1t
	}	
	else {
	    prov1t = exp(-asta*gam)*prov1t
	}
	
	v0tHat = v0t + rnormal(1,1,0,1)*prov0t
	Sigma    = exp(lnsd1)^2, v10, v1t \ v10, exp(lnsd0)^2, v0t    \ v1t, v0t,    1
	SigmaHat = exp(lnsd1)^2, v10, v1t \ v10, exp(lnsd0)^2, v0tHat \ v1t, v0tHat, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) v0t = v0tHat
	    prov0t = exp(gam*(alpha - asta))*prov0t
	}
	else {
	    prov0t = exp(-asta*gam)*prov0t
	}	
			
	/* Collection */
	
	b1Hold = b1Hold \ b1
	b0Hold = b0Hold \ b0
	eHold = eHold \ e
	
	sd1Hold = sd1Hold \ lnsd1
	sd0Hold = sd0Hold \ lnsd0
	v10Hold = v10Hold \ v10
	v1tHold = v1tHold \ v1t
	v0tHold = v0tHold \ v0t	
	
}	
end	
