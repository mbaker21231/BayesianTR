/* In this worksheet, we will describe Bayesian 0,1 endogenous treatment effects */

clear all
use http://www.stata-press.com/data/r13/union3
set more off
keep if union != .
keep if tenure != . 
describe
set seed 50

sum wage
gen highwage = wage > r(mean)

mata:
	st_view(y=., ., "highwage")
	st_view(tr=., ., "union")
	st_view(X=., ., "age grade smsa black tenure")
	st_view(W=., ., "south black tenure")
	X = X, J(rows(y), 1, 1)
	W = W, J(rows(y), 1, 1)
end

probit highwage age grade smsa black tenure if union == 1
mat binit1 = e(b)

probit highwage age grade smsa black tenure if union == 0
mat binit0 = e(b) 

probit union south black tenure
mat einit = e(b)

mata:
    real matrix invnormstab(X) {
	    XHat = editvalue(X, 0, 1e-323)
		XHat = editvalue(XHat, 1, 1-1e-16)
		return(invnormal(XHat))
    }
end

mata: 
	b1 = st_matrix("binit1"), 0, 0
	b0 = st_matrix("binit0"), 0, 0
	e = st_matrix("einit")
end

mata: 
	nb = cols(b1)
	ne = cols(e)
end


mata:
	b1Hold        = J(0, nb, .)
	b0Hold		  = J(0, nb, .)
	eHold         = J(0, ne, .)

	Tb1 = I(nb)*1000
	Tb0 = I(nb)*1000
	Te = I(ne)*1000

	ur = rnormal(rows(y), 1, 0, 1)
	uc = rnormal(rows(y), 1, 0, 1)

	/* our model seems to struggle with identification */
	/* let's see what happens if we try and pin the constant parameter */
	/* at zero...in fact, we are restricting all of the coefficients that might be 
	   trouble to a particular range */
	
	Tb1[6,6] = 1
	Tb0[6,6] = 1
	Tb1[7,7] = 1
	Tb0[7,7] = 1
	Tb1[8,8] = 1000
	Tb0[8,8] = 1000
	
	
end

mata:
draws = 100000
	for (d=1;d<=draws;d++) {

		Xu = (X, ur, uc)
		muXu0 = Xu*b0'
		muXu1 = Xu*b1'
		
		y0Hatnt = muXu0 + invnormstab(normal(-muXu0):*runiform(rows(y),1)):*(1 :- y) + 
		                invnormstab(normal(-muXu0) + (1 :- normal(-muXu0)):*runiform(rows(y),1)):*y
		y0Hatt   = muXu0 + rnormal(rows(y), 1, 0, 1)
		y0 = tr:*y0Hatt + (1 :- tr):*y0Hatnt
  
		y1Hatnt = muXu1 + rnormal(rows(y), 1, 0, 1)
		y1Hatt   = muXu1 + invnormstab(normal(-muXu1):*runiform(rows(y),1)):*(1 :- y) +
		                invnormstab(normal(-muXu1) + (1 :- normal(-muXu1)):*runiform(rows(y),1)):*y
		y1 = tr:*y1Hatt + (1 :- tr):*y1Hatnt
  
		muz = W*e' :+ uc

		et  = invnormstab( normal(-muz) :+ (1 :- normal(-muz)):*runiform(rows(muz),1) )
		ent = invnormstab( normal(-muz):*runiform(rows(muz),1))		
		z = muz :+ et:*tr :+ ent:*(1:-tr)
		
		meane = invsym(W'W + invsym(Te))*W'(z :- uc)
		vare  = invsym(W'W + invsym(Te))
		e = meane + cholesky(vare)*rnormal(cols(e),1,0,1)
		e = e'

		mb1 = invsym(Xu'Xu + invsym(Tb1))*Xu'y1
		mb0 = invsym(Xu'Xu + invsym(Tb0))*Xu'y0
		vb1 = invsym(Xu'Xu + invsym(Tb1))
		vb0 = invsym(Xu'Xu + invsym(Tb0))
		chol1b = cholesky(vb1)
		chol0b = cholesky(vb0)
		b1 = mb1 + chol1b*rnormal(nb, 1, 0, 1)
		b1 = b1'
		b0 = mb0 + chol0b*rnormal(nb, 1, 0, 1)
		b0 = b0'
		
		b1Hold  =  b1Hold \ b1
		b0Hold  = b0Hold \ b0
		eHold   =  eHold \ e
		
		/* Draw latent vars */
		
		pi1r = b1[cols(b1) - 1]
		pi0r = b0[cols(b0) - 1]
		b1z = b1
		b0z = b0
		b1z[cols(b1) - 1] = 0
		b0z[cols(b0) - 1] = 0
		A = (pi1r^2 + pi0r^2 + 1)
		B = (y1 - Xu*b1z')*pi1r + (y0 - Xu*b0z')*pi0r
		ur = B/A + rnormal(rows(y), 1 , 0, 1/A)
		
		pi1c = b1[cols(b1)]
		pi0c = b0[cols(b0)]
		b1z = b1
		b0z = b0
		b1z[cols(b1)] = 0
		b0z[cols(b0)] = 0
		A = (pi1c^2 + pi0c^2 + 2)
		B = (y1 - Xu*b1z')*pi1c + (y0 - Xu*b0z')*pi0c + (z - W*e' )
		uc = B/A + rnormal(rows(y), 1 , 0, 1/A)
	}  	

end


clear
getmata (b1*) = b1Hold
getmata (b0*) = b0Hold
getmata (e*) = eHold

gen t = _n
tsset t
