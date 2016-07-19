clear all
use http://www.stata.com/data/jwooldridge/eacsap/card.dta
set more off
describe

regress lwage educ expersq black south smsa reg661-reg668 smsa66, vce(r)

ivregress 2sls lwage (educ=nearc4) exper expersq black south smsa reg661-reg668, first

/* Get everything into Mata */

mata:
	st_view(y=., ., "lwage")
	st_view(t=., ., "educ")
	st_view(X=., ., "educ exper expersq black south smsa reg661-reg668 smsa66")
	st_view(W=., ., "exper expersq black south smsa reg661-reg668 smsa66 nearc4")
	X = X, J(rows(y), 1, 1)
	W = W, J(rows(y), 1, 1)
end

reg lwage educ exper expersq black south smsa reg661-reg668 smsa66
mat binit = e(b)
reg educ exper expersq black south smsa reg661-reg668 smsa66 nearc4
mat einit = e(b)

mata: 
	b = st_matrix("binit"),0
	e = st_matrix("einit")
end

mata: 
	nb = cols(b)
	ne = cols(e)
end


set seed 54150

mata:
	bHold        = J(0, nb, .)
	eHold        = J(0, ne, .)
	syHold =       J(0, 1, .)
	stHold =       J(0, 1, .)
	piHold =       J(0, 1, .)
	uHold  =       J(rows(y),0,.)

	Tb = I(nb)*1000
	Te = I(ne)*1000

	draws = 4000

	ui = rnormal(rows(y),1,0,1)
	sy2 = 1
	st2 = 1
	
end

mata:
	for (d=1;d<=draws;d++) {
		bOnly =b[1..nb-1]
		pi = b[nb]
		resy = y-X*bOnly'
		rest = t-W*e'
		A = pi^2/sy2 + 1/st2 + 1
		B = pi*resy/sy2 + rest/st2
		mi = B:/A
		sd = 1/A
		ui = mi + rnormal(rows(mi),1,0,1)*sqrt(sd)
    
		Xu = (X,ui)
		mb = invsym(Xu'Xu + invsym(Tb))*Xu'y
		vb = sy2*invsym(Xu'Xu + invsym(Tb))
		cholb = cholesky(vb)
		b = mb + cholb*rnormal(nb, 1, 0, 1)
		b = b'
                                  
		tHat = t - ui
		me = invsym(W'W + invsym(Te))*W'tHat
		ve = st2*invsym(W'W + invsym(Te))
		chole = cholesky(ve)
		e = me +chole*rnormal(ne, 1, 0, 1)
		e = e'

		shp = rows(y)/2
		scl = (y-Xu*b')'(y-Xu*b')/2
		scl = 1/scl
	    sy2 = rgamma(1, 1, shp, scl)

		scl = (t - ui - W*e')'(t - ui - W*e')/2
		scl = 1/scl
	    st2 = rgamma(1, 1, shp, scl)
                             
		bHold  =  bHold \ b
		eHold  =  eHold \ e
		syHold =  syHold \ sy2
		stHold =  stHold \ st2
		piHold =  piHold \ b[nb]
	}  	  
end

clear
getmata (b*) = bHold
getmata (e*) = eHold
getmata syHold
getmata stHold
getmata piHold

gen t = _n
tsset t
sum

hist pi if t>=7000

tsline syHold stHold piHold if t>2999
