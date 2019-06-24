/* In this worksheet, we will describe Bayesian 0,1 endogenous treatment effects */
/* This gives us personalized treatment effects and also gives us variances, etc. */
/* Highly parametric, but good in small samples. */

/* After having experimented with doing this with latent variables, we decided      */
/* it would work better without doing latent variables. This is because the process */
/* just never seemed to converge. In this version, we have five things to draw:     */
/* variance in outcome 1, variance in outcome 0, covariance between 1 and 0,        */
/* covariance between 1 and treatment, and covariance between 2 and treatment       */

set seed 5150
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

/* Initial values for variance-matrix parameters */

mata:
	v11 = 1
	v00 = 1
	v10 = 0
	v1t = 0
	v0t = 0	
end

/* Initial ancilliary parameters */

mata:    
    av11 = 0
	av00 = 0
	av10 = 0
	av1t = 0
	av0t = 0
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
		((part3,part1,part2))
        return(part1 + part2 + part3)
    }
end

/* Initial values for all the latent variables */

mata:
	y0Hat = X*b0' + rnormal(rows(y),1,0,1)*v11
	y1Hat = X*b1' + rnormal(rows(y),1,0,1)*v00
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
	vHold  = v11, v10, v00, v1t, v0t
	avHold = av11, av10, av00, av1t, av0t

end	

/* Placeholders for proposal distributions and all that: */

mata:
    proav11 = 1
	proav00 = 1
	proav10 = 1
	proav1t = 1
	proav0t = 1
	
	delta =  1
	asta  = .4
end
	
mata:

for (i=1;i<=1000;i++) {

/* First part - regression parameters and latent values for y1 */
	
	Sig12 = (v10, v1t)
	Sig22 = v00, v0t \ v0t, 1
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey0, et) )
	CV = v11 - Sig12*Sig22m1*Sig12'
	
	mc1   = m1 + CM
	y1Hat = mc1 + rnormal(rows(y),1,0,1)*sqrt(CV)
	y1    = tr:*y + (1 :-tr):*y1Hat
	
	mb1 = XX*X'(y1 - CM)
	vb1 = CV*XX
	b1 = mb1 + cholesky(vb1)*rnormal(cols(b1), 1, 0, 1)
	b1 = b1'	

	m1 = X*b1'
	ey1 = (y1 - m1)
	
/* Second part - regression parameters and latent values for y0 */

	Sig12 = (v10, v0t)
	Sig22 = v11, v1t \ v1t, 1
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey1, et) )
	CV = v00 - Sig12*Sig22m1*Sig12'
	
	mc0   = m0 + CM
	y0Hat = mc0 + rnormal(rows(y),1,0,1)*sqrt(CV)
	y0 = (1 :- tr):*y + tr:*y0Hat
	
	mb0 = XX*X'(y0 - CM)
	vb0 = CV*XX
	b0 = mb0 + cholesky(vb0)*rnormal(cols(b0), 1, 0, 1)
	b0 = b0'	

	m0 = X*b0'
	ey0 = (y0 - m0)	
	
/* Third part - regression parameters and latent values for z */
	
	Sig12 = (v1t, v0t)
	Sig22 = v11, v10 \ v10, v00
	Sig22m1 = invsym(Sig22)
	
	CM = rowsum( (Sig12*Sig22m1):*(ey1, ey0) )
	CV = 1 - Sig12*Sig22m1*Sig12'	
	
	mct = mt + CM
	et  = CV*invnormstab( normal(-mct/CV) + (1 :- normal(-mct/CV)):*runiform(rows(mct),1) )
	ent = CV*invnormstab( normal(-mct/CV):*runiform(rows(mct),1) )
	z = mct + et:*tr + ent:*(1 :- tr)
	
	/* Draw eta coefficients given unit variance, common values, etc. */
	
	meane = WW*W'(z - CM)
	vare  = CV*WW
	e = meane + cholesky(vare)*rnormal(cols(e), 1, 0, 1)
	e = e'
	
	mt = W*e'
	et = (z - mt)
	
	/* Next: variance parameters */

	gam = 1/i^delta
	
	/* ancparm1 */
	
	lb = (v00*v1t^2 - 2*v10*v0t*v1t + v10^2) / (v00 - v0t^2)

	av11Hat = av11 + rnormal(1,1,0,1)*proav11
	v11Hat  = lb + exp(av11Hat)
	
	Sigma    = v11,    v10, v1t \ v10, v00, v0t \ v1t, v0t, 1
	SigmaHat = v11Hat, v10, v1t \ v10, v00, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) {
		    av11 = av11Hat
			v11  = v11Hat
		}
	    proav11 = exp(gam*(alpha - asta))*proav11
	}
	else {
	    proav11 = exp(-asta*gam)*proav11
	}
	
	/* v00 */
	
	lb = (v11*v0t^2 - 2*v10*v0t*v1t + v10^2) / (v11 - v1t^2)

	av00Hat = av00 + rnormal(1,1,0,1)*proav00
	v00Hat  = lb + exp(av00Hat)
	Sigma    = v11, v10, v1t \ v10, v00,    v0t \ v1t, v0t, 1
	SigmaHat = v11, v10, v1t \ v10, v00Hat, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) {
		    av00 = av00Hat
			v00  = v00Hat
		}
	    proav00 = exp(gam*(alpha - asta))*proav00
	}
	else {
	    proav00 = exp(-asta*gam)*proav00	
	}
	
	/* v10 */
	
	st = sqrt(v0t^2*v1t^2 - v0t^2*v11 - v1t^2*v00 + v00*v11)
	lb = v0t*v1t - st
	ub = v0t*v1t + st

	av10Hat = av10 + rnormal(1,1,0,1)*proav10
	v10Hat = 1/(1 + exp(av10Hat))*lb + exp(av10Hat)/(1 + exp(av10Hat))*ub
	
	Sigma    = v11, v10,    v1t \ v10,    v00, v0t \ v1t, v0t, 1
	SigmaHat = v11, v10Hat, v1t \ v10Hat, v00, v0t \ v1t, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) {
		    av10 = av10Hat
             v10 = v10Hat		
		}
	    proav10 = exp(gam*(alpha - asta))*proav10
	}
	else {
	    proav10 = exp(-asta*gam)*proav10
	}
	
	/* v1t */
	
	st = sqrt(v0t^2*v10^2 - v00*v0t^2*v11 + v00^2*v11 - v00*v10^2)
	lb = (v10*v0t - st)/v00
	ub = (v10*v0t + st)/v00

	av1tHat = av1t + rnormal(1,1,0,1)*proav1t
	v1tHat = 1/(1 + exp(av1tHat))*lb + exp(av1tHat)/(1 + exp(av1tHat))*ub
	
	Sigma    = v11, v10, v1t    \ v10, v00, v0t \ v1t,    v0t, 1
	SigmaHat = v11, v10, v1tHat \ v10, v00, v0t \ v1tHat, v0t, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) {
		    av1t = av1tHat
             v1t = v1tHat		
		}
	    proav1t = exp(gam*(alpha - asta))*proav1t
	}	
	else {
	    proav1t = exp(-asta*gam)*proav1t
	}
	
	/* v0t */
	
	st = sqrt(v1t^2*v10^2 - v00*v1t^2*v11 + v00*v11^2 - v11*v10^2)
	lb = (v10*v1t - st)/v11
	ub = (v10*v1t + st)/v11

	av0tHat = av0t + rnormal(1,1,0,1)*proav0t
	 v0tHat = 1/(1 + exp(av0tHat))*lb + exp(av0tHat)/(1 + exp(av0tHat))*ub
	
	Sigma    = v11, v10, v1t \ v10, v00, v0t    \ v1t, v0t,    1
	SigmaHat = v11, v10, v1t \ v10, v00, v0tHat \ v1t, v0tHat, 1 

	if ( hasmissing(cholesky(SigmaHat)) == 0 ) {
	    val    = ln_L((ey1, ey0, et), Sigma)
		valHat = ln_L((ey1, ey0, et), SigmaHat)
        rat = valHat - val
	    alpha = min((exp(rat), 1))
	    if (runiform(1,1,0,1) < alpha) {
           av0t = av0tHat		
		    v0t =  v0tHat

		}
	    proav0t = exp(gam*(alpha - asta))*proav0t
	}
	else {
	    proav0t = exp(-asta*gam)*proav0t
	}	
			
	/* Collection */
	
	b1Hold = b1Hold \ b1
	b0Hold = b0Hold \ b0
	eHold = eHold \ e
	vHold  = vHold \ (v11, v10, v00, v1t, v0t)
	avHold = avHold \ (av11, av10, av00, av1t, av0t)
	
}	
end	

preserve
clear 
getmata (e*) =  eHold
getmata (b1*) = b1Hold
getmata (b0*) = b0Hold

sum e* b*


