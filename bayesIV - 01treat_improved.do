/* In this worksheet, we will describe Bayesian 0,1 endogenous treatment effects */

clear all
use http://www.stata-press.com/data/r13/union3
set more off
keep if union != .
keep if tenure != . 
describe
set seed 50


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

/* New likelihood to include more stuff... */
/* We pass along six total parameters...   */
/* We also along two latent random effects */
/* Everything else is handled before this  */
/* Prior_piv should contain four terms     */
/* T is an indicator as to whether or not  */
/* the treatment was received              */
/* Note that we also add a "Flag" variable */
/* which lets us know if we should return  */
/* the sum or the individual row elements  */
/* the latter will help drawing individual */
/* level parameters                        */

mata:
	real matrix llFun(s, y1, y0, Xu, b1, b0, e, W, T, toggle) {
		ln_s1=s[1]
		ln_s0=s[2]
		Uc = Xu[,cols(Xu)]
		Ur = Xu[,cols(Xu)-1]
		
		chk = sqrt(2)*(2*T :- 1)
		P1 = normal(chk:*(W*e' :+ Uc))
		P1 = P1 :+ (P1 :== 0):*(1e-323)
		L1 = ln(P1)
		
		L2 = -(y1 - Xu*b1'):^2/(2*exp(ln_s1)^2) :- ln_s1 :+
			 -(y0 - Xu*b0'):^2/(2*exp(ln_s0)^2) :- ln_s0
		L3 = -Uc:^2
		L4 = -Ur:^2/2
		L = L1 :+ L2 :+ L3 :+ L4

		if (hasmissing(L) & toggle == 1) {
			printf("ohshit!")			
		}
		else if (hasmissing(L)) {
			printf("ohno!")
			return(.)
		}
		else if (toggle == 1) return(L)
		else return(sum(L))
	}
end

mata:
    real matrix invnormstab(X) {
	    XHat = editvalue(X, 0, 1e-323)
		XHat = editvalue(X, 1, 1-1e-16)
		return(invnormal(XHat))
    }


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
	sy1Hold       = J(0, 1, .)
	sy0Hold       = J(0, 1, .)
	uHold         = J(rows(y),0,.)

	prosy1v = 1
	prosy0v = 1

	priorpiv = 100,100,100,100
	
	Tb = I(nb)*100
	Te = I(ne)*100

	asta = .24

	b1 = b1
	b0 = b0

	e = e

	ur = rnormal(rows(y), 1, 0, 1)
	uc = rnormal(rows(y), 1, 0, 1)

	sy12 = 1
	sy02 = 1


end

mata:
draws = 10000
	for (d=1;d<=draws;d++) {

		Xu = (X, ur, uc)
		y0Hat = Xu*b0' :+ rnormal(rows(y) , 1, 0, 1):*sqrt(sy02)
		y1Hat = Xu*b1' :+ rnormal(rows(y) , 1, 0, 1):*sqrt(sy12)
		y1 = tr:*y :+ (1 :- tr):*y1Hat
		y0 = (1 :- tr):*y :+ tr:*y0Hat

		muz = W*e' :+ uc
		et  = invnormstab( normal(-muz) :+ (1 :- normal(-muz)):*runiform(rows(muz),1) )
		ent = invnormstab( normal(-muz):*runiform(rows(muz),1))		
		z = muz :+ et:*tr :+ ent:*(1:-tr)
		meane = invsym(W'W + invsym(Te))*W'(z :- uc)
		vare  = invsym(W'W + invsym(Te))
		e = meane + cholesky(vare)*rnormal(cols(e),1,0,1)
		e = e'

		Xu = X, ur, uc
		mb1 = invsym(Xu'Xu + invsym(Tb))*Xu'y1
		mb0 = invsym(Xu'Xu + invsym(Tb))*Xu'y0
		vb1 = sy12*invsym(Xu'Xu + invsym(Tb))
		vb0 = sy02*invsym(Xu'Xu + invsym(Tb))
		chol1b = cholesky(vb1)
		chol0b = cholesky(vb0)
		b1 = mb1 + chol1b*rnormal(nb, 1, 0, 1)
		b1 = b1'
		b0 = mb0 + chol0b*rnormal(nb, 1, 0, 1)
		b0 = b0'
                                  
		shp = rows(y)/2
		scl = (y-Xu*b1')'(y-Xu*b1')/2
		scl = 1/scl
	    sy12 = rgamma(1, 1, shp, scl)

		scl = (y-Xu*b0')'(y-Xu*b0')/2
		scl = 1/scl
	    sy02 = rgamma(1, 1, shp, scl)
		
		b1Hold  =  b1Hold \ b1
		b0Hold  = b0Hold \ b0
		eHold   =  eHold \ e
		sy1Hold =  sy1Hold \ sy12
		sy0Hold =  sy0Hold \ sy02
		
		/* Draw latent vars */
		
		pi1r = b1[cols(b1) - 1]
		pi0r = b0[cols(b0) - 1]
		b1z = b1
		b0z = b0
		b1z[cols(b1) - 1] = 0
		b0z[cols(b0) - 1] = 0
		A = (pi1r^2/sy12 + pi0r^2/sy02 + 1)
		B = (y1 - Xu*b1z')*pi1r/sy12 + (y0 - Xu*b0z')*pi0r/sy02
		ur = B/A + rnormal(rows(y), 1 , 0, 1/A)
		
		pi1c = b1[cols(b1)]
		pi0c = b0[cols(b0)]
		b1z = b1
		b0z = b0
		b1z[cols(b1)] = 0
		b0z[cols(b0)] = 0
		A = (pi1c^2/sy12 + pi0c^2/sy02 + 2)
		B = (y1 - Xu*b1z')*pi1c/sy12 + (y0 - Xu*b0z')*pi0c/sy02 + (z - W*e' )
		uc = B/A + rnormal(rows(y), 1 , 0, 1/A)
		
	}  	

end


clear
getmata (b1*) = b1Hold
getmata (b0*) = b0Hold
getmata (e*) = eHold
getmata sy1Hold
getmata sy0Hold

gen t = _n
tsset t
tsline sy1Hold sy0Hold e1 e2 e3 e4

/* Let's try and make all the components of the correlation matrix */

rename b17 pi1r
rename b07 pi0r
rename b18 pi1c
rename b08 pi0c

gen sigma1 = exp(sy1Hold)
gen sigma0 = exp(sy0Hold)

*Entry 1,1
gen Sigma11 = sigma1^2 + pi1r^2 + pi1c^2
hist Sigma11 if t> 8999, name(Entry11, replace)


*Entry 1,2
gen Sigma12 = pi1c*pi0c + pi1r*pi0r
hist Sigma12 if t > 8999, name(Entry12, replace)

*Entry 1,3
gen Sigma13 = pi1c
hist Sigma13 if t > 8999, name(Entry13, replace)

*Entry 2,2
gen Sigma22 = sigma0^2 + pi0r^2 + pi0c^2
hist Sigma22 if t > 8999, name(Entry22, replace)

*Entry 2,3
gen Sigma23 = pi0c
hist pi0c if t > 8999, name(Entry23, replace)

gen Ones = 1
hist Ones if t > 8999, name(Ones, replace)

graph combine Entry11 Entry12 Entry13 Entry12 Entry22 Entry23 Entry13 Entry23 Ones

/* Make an easier-to-get matrix */

mata: 




