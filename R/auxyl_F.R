unimcd<-function(y,h,len){
#This code calls the UNIMCD.C routine (it's a C translation of the one found in robustbase). 
	inloc<-inmean<-incov<-0.0;
	out<-.C("R_unimcd",as.double(y),as.integer(h),as.integer(len),as.integer(length(y)),as.double(inmean),as.double(incov),as.integer(inloc),PACKAGE="DetR")
	list(initmean=as.numeric(out[[5]]),initcov=as.numeric(out[[6]]),iloc=as.numeric(out[[7]]))
}
fxOGK<-function (Data, scale_est = "Qn", intercept=1,h,doCsteps=1) {
    p <- ncol(Data)
    n <- nrow(Data)
	h0<-ceiling(n/2)
	H0<-rep(0,h0)
    lh <- length(h)
    Q <- matrix(0,n,lh)
	citer<-ovec<-rep(0,lh)
    Mtype <- match(scale_est, c("Qn", "scaleTau2"))[1]
    if (is.na(Mtype)) 
        stop("scale_est should be one of Qn or scaleTau2.")
    Mtype <- Mtype - 1
	W<-rep(0,p)
    fit2 <- .C("R_FastOGK",
		as.integer(n),			#01
		as.integer(p),			#02
		as.double(Data), 		#03
		as.integer(Q),			#04
		as.integer(Mtype),		#05
		as.integer(intercept),		#06
		as.integer(W),			#07
		as.integer(h),			#08
		as.integer(lh),			#09
		as.integer(H0),			#10
		as.integer(h0),			#11
		as.double(ovec),		#12
		as.integer(doCsteps),		#13
		as.integer(citer),		#14
	PACKAGE="DetR")
	list(bestCStep=fit2[[4]],HasZeroScale=fit2[[7]],bestRaw=fit2[[10]],Objective=fit2[[12]],citer=fit2[[14]])
}
quanf<-function(alpha,n,p) return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
#This code computes 'h' (it's an R translation of the one found in LIBRA).
DetLTS_raw<-function (x, y, intercept = 1, alpha = 0.75, h = NULL, scale_est = "scaleTau2",doCsteps=1){
    y <- as.matrix(y)
    x <- as.matrix(x)
    x <- data.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    q <- ncol(y)
    if (is.numeric(intercept) == 0) 
        stop("intercept should be set to either 1 or 0.")
    if (sum(intercept == c(0, 1)) == 0) 
        stop("intercept should be set to either 1 or 0.")
    if (q > 1) 
        stop("y is not one-dimensional.")
    na.x <- complete.cases(x)
    na.y <- complete.cases(y)
	if(min(c(sum(na.y),sum(na.x)))<n)	stop("x or y contain missing data. Remove the rows with missing data.")
    Mtype <- match(scale_est, c("Qn", "scaleTau2"))[1]
    if (is.na(Mtype)) 
        stop("scale_est should be one of Qn or scaleTau2.")
    if (sum(na.x) != sum(na.y)) 
        stop("Number of observations in x and y are not equal.")
    if (sum(na.x) < (p + 1)) 
        stop("Not enough observations with non-missing values.")
    n <- nrow(x)
    p <- ncol(x)
    hf <- ceiling((n + p + 1)/2)
    if (is.numeric(alpha) & !is.numeric(h)) {
        alpha <- sort(alpha)
        if (min(alpha) >= 0.5 & max(alpha) <= 1) {
            h <- sort(quanf(alpha, n = n, p = p + intercept))
        }
        else {
            stop("Error: invalid alpha value")
        }
    }
    if (is.numeric(h)) {
        h <- sort(h)
        if (min(h) < hf | max(h) > n) 
            stop(paste("The smallest h should be at least ", 
                hf, " and at most ", n, sep = ""))
    }
    if (is.null(alpha) & is.null(h)) {
        alpha <- 0.75
        h <- quanf(alpha, n = n, p = p + intercept)
        print("alpha was set to 0.75 [default]")
    }
    if (!is.numeric(alpha) & !is.numeric(h)) {
        alpha <- 0.75
        h <- quanf(alpha, n = n, p = p + intercept)
        print("alpha was set to 0.75 [default]")
    }
    if (min(h) < n) {
        out2 <- detlts_in9(x = x, y = y, intercept = intercept, 
            scale_est = scale_est,h=h,doCsteps=doCsteps)
    }

#    if (length(out2) == 1) 
#        out2 <- out2[[1]]
    out2
}
ordreg<-function(x,y,intercept){
	#this is classical OLS, but with the same outputs as DetLTS 
	#(so that it is compatible with the auxiliary functions for 
	#plotting the outputs of DetLTS for several values of h).
	#Mia asked for this. The variable names in ordereg are the 
	#the same as in LIBRA (if the names do not conflict with 
	#other R protected functions --otherwise I capitalized the 
	#function names). 
	n<-nrow(x)
	if(intercept)	x<-cbind(1,x)
	p<-ncol(x)
	regres<-lm(y~x-1)
	s0a<-crossprod(regres$resid)
	s0<-summary(regres)$sigma;
	if(s0>1e-7){							#KV NEW
		weights<-(abs(regres$resid/s0)<=qnorm(0.9875));
		weights<-as.numeric(weights)
		regres<-lm(y~x-1,weights=weights)
	} 
	list(raw.coefficients=coef(regres),
	objective=s0a,							#KV NEW
	wt=rep(1,n),
	fitted=fitted(regres),
	res=residuals(regres),
	scale=sqrt(sum(residuals(regres)**2)/(n-p)),			#KV NEW
	rsquared=summary(regres)$r.squared,
	h=n,
	Hsubsets=1:n,
	rd=residuals(regres)**2,
	cutoff=qnorm(0.9875),
	flag=rep(1,n),X=x,y=y)
}
detlts_in9<-function (x, y, intercept, scale_est,h,doCsteps=doCsteps) {
#x=x;y=y;intercept=1;alpha=0.5;scale_est="scaleTau2";h<-c(66,76,86)
    n <- nrow(x)
    p <- ncol(x)
    if (svd(scale(x))$d[p] < 1e-07) 
        stop("x is singular")
#    if (intercept) {
#        datamed <- c(colMedians(x), median(y))
#        x <- sweep(x, 2, datamed[1:p], check.margin = FALSE)
#        y <- y - datamed[p + 1]
#    }
#    datao <- cbind(x, y)
#    datasca <- apply(datao, 2, scale_est)
#    x <- sweep(x, 2, datasca[1:p], FUN = "/", check.margin = FALSE)
#    y <- y/datasca[p + 1]
    datao <- cbind(x, y)
    CXY <- fxOGK(Data=datao,scale_est=scale_est,intercept=intercept,h=h,doCsteps=doCsteps)
	if(sum(CXY$HasZeroScale[1:p])>0){
		stop(paste0("Column ",which(CXY$HasZeroScale>0)," of the design matrix have 0 value of ", scale_est))
	}
	if(CXY$HasZeroScale[p+1]>0){
		stop(paste0("The response vector has 0 value of ", scale_est))
	}
	lh<-length(h)
	a2<-vector("list",4)
	a2[[1]]<-vector("list",lh);
	a2[[2]]<-CXY$Objective
	if(lh>1){
		a1<-matrix(CXY$bestCStep,ncol=lh)
		for(i in 1:lh)	a2[[1]][[i]]<-a1[1:h[i],i]
	} else {
		a2[[1]][[1]]<-CXY$bestCStep[1:h]
	}
	names(a2)<-c("Subset","Objective","Raw","SubsetSize")
	a2[[4]]<-h
	names(a2[[1]])<-c(paste0("ActiveSubsetSize_",h))
	a2[[3]]<-CXY$bestRaw
    return(a2)
}
ltscheckout<-function(x,y,inbest,h,intercept,alpha,use.correction=TRUE,objfct){
	#this code does the reweighting, consistency and small sample correction given 
	#a list containing the indexes of the h observations awarded a positive weight 
	#at the end of the lts algorithm.

	#inbest:	a list containing the h indexes in \{1:n\} of those observations awarded a positive weight by lts.
	#objfct: 	so that this function returns the objective function of the subset before the reweighting/adjustments 
	#		(as in robustbase/LIBRA).

	#the rest of the code/variable names are taken from robustbase::ltsReg (see maintenair in citation("robustbase")).
	#the only differences is that I'm abiding by the R convention of putting the intercept first (instead of last as 
	#is dones in LIBRA and robustbase). 
 	raw.cnp2<-rep(1,2)
	cnp2<-rep(1,2)
	quantiel<-qnorm(0.9875)
	x<-cbind("Intercept"=1,x)
	n<-nrow(x)
	p<-ncol(x)
	coefs<-rep(NA,p)
	piv<-1:p
	ans<-list(alpha=alpha)
	ans$best<-sort(inbest)
	Fitt<-lm(y[inbest]~x[inbest,]-1)
	cf<-Fitt$coef
	fitted<-x%*%cf			#KV NEW
	resid<-y-fitted
	coefs[piv]<-cf
	ans$raw.coefficients<-coefs
	ans$quan<-h
	correct<-if(use.correction) LTScnp2(p,intercept=intercept,n,alpha) else 1
	raw.cnp2[2]<-correct
	s0<-sqrt((1/h)*sum(sort(resid^2,partial=h)[1:h]))
#	s0a<-which(resid**2<=quantile(resid**2,(h-1)/n,type=1)); ##Here for KV.NEW
#	s0<-sqrt(mean(resid[s0a]^2))
	sh0<-s0
	qn.q<-qnorm((h+n)/(2*n))
	s0<-s0/sqrt(1-(2*n)/(h/qn.q)*dnorm(qn.q))*correct
	if(abs(s0)<1e-07){
		ans$raw.weights<-weights<-as.numeric(abs(resid)<=1e-07)
		ans$scale<-ans$raw.scale<-0
		ans$coefficients<-ans$raw.coefficients
	} else {
		ans$raw.scale<-s0
		ans$raw.resid<-resid/ans$raw.scale
		ans$raw.weights<-weights<-as.numeric(abs(resid/s0)<=quantiel)
		sum.w<-sum(weights)
		## old,suboptimal: z1<-lsfit(x,y,wt=weights,intercept=FALSE)
		z1<-lm.wfit(x,y,w=weights)
		ans$coefficients<-z1$coef
		fitted<-x%*%z1$coef
		resid<-z1$residuals
		ans$scale<-sqrt(sum(weights*resid^2)/(sum.w-1))
		if (sum.w==n) {
		cdelta.rew<-1
		correct.rew<-1
		} else {
			qn.w<-qnorm((sum.w+n)/(2*n))
			cnp2[1]<-cdelta.rew<-1/sqrt(1-(2*n)/(sum.w/qn.w)*dnorm(qn.w))
			correct.rew<-if(use.correction) LTScnp2.rew(p,intercept=intercept,n,alpha) else 1
			cnp2[2]<-correct.rew
			ans$scale<-ans$scale*cdelta.rew*correct.rew
		}
		ans$resid<-resid/ans$scale
		weights<-as.numeric(abs(ans$resid)<=quantiel)
	}
	## unneeded: names(ans$coefficients)<-names(ans$raw.coefficients)
	ans$crit<-objfct
	if (intercept) {
		sh<-unimcd(y=y,h=h,len=n-h+1)$initcov
		iR2<-(sh0**2/sh)
	} else {
		#s1<-sum(sort(resid^2,partial=h)[1:h])
		s0a<-which(resid**2<=quantile(resid**2,h/n,type=1)); 	##Here for KV NEW
		s1<-sum(resid[s0a]^2)
#		sh<-sum(sort(y^2,partial=h)[1:h])
		s0b<-which(y**2<=quantile(y**2,h/n,type=1)); 		##Here for KV NEW
		sh<-sum(y[s0b]^2)
		iR2<-s1/sh
	}
	ans$rsquared<-if(is.finite(iR2)) max(0,min(1,1-iR2)) else 0
	attributes(resid)<-attributes(fitted)<-attributes(y)
	ans$method<-"Least Trimmed Squares Robust Regression."

	ans$intercept<-intercept
	if(abs(s0)<1e-07) ans$method<-paste(ans$method,"\nAn exact fit was found!")
	
	ans$lts.wt<-weights
	ans$residuals<-resid
	ans$fitted.values<-fitted
	ans$raw.cnp2<-raw.cnp2
	ans$cnp2<-cnp2
	class(ans)<-"lts"
	return(ans)
}
LTScnp2<-function(p,intercept=intercept,n,alpha){
	#This function was taken from robustbase. See citation("robustbase").
	stopifnot(0.5<=alpha,alpha<=1)
	if(intercept)	p<-p-1
	stopifnot(p==as.integer(p),p>=0)
	if(p==0){
		fp.500.n<-1-exp(0.262024211897096)/n^0.604756680630497
		fp.875.n<-1-exp(-0.351584646688712)/n^1.01646567502486
		if((0.5<=alpha) &&(alpha<=0.875)){
			fp.alpha.n<-fp.500.n+(fp.875.n-fp.500.n)/0.375*(alpha-0.5)
			fp.alpha.n<-sqrt(fp.alpha.n)
		}
		if((0.875 < alpha) &&(alpha < 1)){
			fp.alpha.n<-fp.875.n+(1-fp.875.n)/0.125*(alpha-0.875)
			fp.alpha.n<-sqrt(fp.alpha.n)
		}
	} else { ## p>=1
		if(p==1){
			if(intercept){
				fp.500.n<-1-exp(0.630869217886906 )/n^0.650789250442946
				fp.875.n<-1-exp(0.565065391014791 )/n^1.03044199012509
			}else {
				fp.500.n<-1-exp(-0.0181777452315321)/n^0.697629772271099
				fp.875.n<-1-exp(-0.310122738776431 )/n^1.06241615923172
			}
		} else { ## --- p > 1 ---
			if(intercept){
				##			 "alfaq"		"betaq"	   "qwaarden"
				coefgqpkwad875<-matrix(c(-0.458580153984614,1.12236071104403,3,-0.267178168108996,1.1022478781154,5),ncol=2)
				coefeqpkwad500<-matrix(c(-0.746945886714663,0.56264937192689,3,-0.535478048924724,0.543323462033445,5),ncol=2)
			} else {
				##			 "alfaq"		"betaq"	   "qwaarden"
				coefgqpkwad875<-matrix(c(-0.251778730491252,0.883966931611758,3,-0.146660023184295,0.86292940340761,5),ncol=2)
				coefeqpkwad500<-matrix(c(-0.487338281979106,0.405511279418594,3,-0.340762058011,0.37972360544988,5),ncol=2)
			}
			y.500<-log(- coefeqpkwad500[1,]/p^coefeqpkwad500[2,])
			y.875<-log(- coefgqpkwad875[1,]/p^coefgqpkwad875[2,])

			A.500<-cbind(1,-log(coefeqpkwad500[3,]*p^2))
			coeffic.500<-solve(A.500,y.500)
			A.875<-cbind(1,-log(coefgqpkwad875[3,]*p^2))

			coeffic.875<-solve(A.875,y.875)
			fp.500.n<-1-exp(coeffic.500[1])/n^coeffic.500[2]
			fp.875.n<-1-exp(coeffic.875[1])/n^coeffic.875[2]
		}
		if(alpha<=0.875)
			fp.alpha.n<-fp.500.n+(fp.875.n-fp.500.n)/0.375*(alpha-0.5)
		else ##	 0.875 < alpha<=1
			fp.alpha.n<-fp.875.n+(1-fp.875.n)/0.125*(alpha-0.875)
	}## else(p>=1)
	return(1/fp.alpha.n)
} ## LTScnp2
LTScnp2.rew<-function(p,intercept=intercept,n,alpha){
	#This function was taken from robustbase. See citation("robustbase").
	stopifnot(0.5<=alpha,alpha<=1)
	if(intercept) p<-p-1
	stopifnot(p==as.integer(p),p>=0)
	if(p==0){
		fp.500.n<-1-exp(1.11098143415027)/n^1.5182890270453
		fp.875.n<-1-exp(-0.66046776772861)/n^0.88939595831888

		if(alpha<=0.875)
			fp.alpha.n<-fp.500.n+(fp.875.n-fp.500.n)/0.375*(alpha-0.5)
		else ##	 0.875 < alpha<=1
			fp.alpha.n<-fp.875.n+(1-fp.875.n)/0.125*(alpha-0.875)
		## MM: sqrt(){below} is ''different logic'' than below..(??)
		fp.alpha.n<-sqrt(fp.alpha.n)
	} else {
		if(p==1){
			if(intercept){
				fp.500.n<-1-exp(1.58609654199605 )/n^1.46340162526468
				fp.875.n<-1-exp(0.391653958727332)/n^1.03167487483316
			} else {
				fp.500.n<-1-exp(0.6329852387657)/n^1.40361879788014
				fp.875.n<-1-exp(-0.642240988645469)/n^0.926325452943084
			}
		} else { ##  --- p > 1 ---
			if(intercept){
			##			 "alfaq"		"betaq"	   "qwaarden"
				coefqpkwad875<-matrix(c(-0.474174840843602,1.39681715704956,3,-0.276640353112907,1.42543242287677,5),ncol=2)
				coefqpkwad500<-matrix(c(-0.773365715932083,2.02013996406346,3,-0.337571678986723,2.02037467454833,5),ncol=2)
			} else {
				##			 "alfaq"		"betaq"	   "qwaarden"
				coefqpkwad875<-matrix(c(-0.267522855927958,1.17559984533974,3,-0.161200683014406,1.21675019853961,5),ncol=2)
				coefqpkwad500<-matrix(c(-0.417574780492848,1.83958876341367,3,-0.175753709374146,1.8313809497999,5),ncol=2)
			}
			y.500<-log(-coefqpkwad500[1,]/p^coefqpkwad500[2,])
			y.875<-log(-coefqpkwad875[1,]/p^coefqpkwad875[2,])
			A.500<-cbind(1,-log(coefqpkwad500[3,]*p^2))
			coeffic.500<-solve(A.500,y.500)
			A.875<-cbind(1,-log(coefqpkwad875[3,]*p^2))
			coeffic.875<-solve(A.875,y.875)
			fp.500.n<-1-exp(coeffic.500[1])/n^coeffic.500[2]
			fp.875.n<-1-exp(coeffic.875[1])/n^coeffic.875[2]
		}
		if(alpha<=0.875)
			fp.alpha.n<-fp.500.n+(fp.875.n-fp.500.n)/0.375*(alpha-0.5)
		else ##	 0.875 < alpha<=1
			fp.alpha.n<-fp.875.n+(1-fp.875.n)/0.125*(alpha-0.875)
	}## else(p>=1)
	return(1/fp.alpha.n)
} ## LTScnp2.rew
#lightly modified code, originally 
#from http://www.stat.ubc.ca/~matias/fasts.txt
our.solve <- function(a,b) {
    a <- qr(a)
    da <- dim(a$qr)
    if(a$rank < (p <- da[2]))
        return(NA)
    else qr.coef(a, b)
}
rho_fastS <- function(u, cc=cc) {
    w <- abs(u)<=cc
    v <- (u^2/(2)*(1-(u^2/(cc^2))+(u^4/(3*cc^4))))*w +(1-w)*(cc^2/6)
    v <- v*6/cc^2
    return(v)
}
norm_fastS <- function(x) sqrt( sum( x^2 ) )
loss.S <- function(u,s,cc) mean(rho_fastS(u/s,cc) )
f.w <- function(u, cc){
        # weight function = psi(u)/u
        tmp <- (1 - (u/cc)^2)^2
        tmp[ abs(u/cc) > 1 ] <- 0
        return(tmp)
}
re.s <- function(x,y,initial.beta,initial.scale,k,conv,b,cc) {
    # does "k" IRWLS refining steps from "initial.beta"
    #
    # if "initial.scale" is present, it's used, o/w the MAD is used
    # k = number of refining steps
    # conv = 0 means "do k steps and don't check for convergence"
    # conv = 1 means "stop when convergence is detected, or the
    #                 maximum number of iterations is achieved"
    # b and cc = tuning constants of the equation
    # 

  

    n <- dim(x)[1]
    p <- dim(x)[2]
    res <- y - x %*% initial.beta
    if( missing( initial.scale ) ) 
        initial.scale <- scale <- median(abs(res))/.6745
    else
        scale <- initial.scale
    
    if( conv == 1) k <- 50
    #
    # if conv == 1 then set the max no. of iterations to 50
    # magic number alert!!!
    
    beta <- initial.beta

    lower.bound <- median(abs(res))/cc

    
    for(i in 1:k) {
        # do one step of the iterations to solve for the scale
        scale.super.old <- scale
        #lower.bound <- median(abs(res))/1.56
        scale <- sqrt( scale^2 * mean( rho_fastS( res / scale, cc ) ) / b     )
        # now do one step of IRWLS with the "improved scale"
        weights <- f.w( res/scale, cc )
        W <- matrix(weights, n, p)
        xw <- x * sqrt(W)
        yw <- y *   sqrt(weights)
        beta.1 <- our.solve( t(xw) %*% xw ,t(xw) %*% yw )
        if(any(is.na(beta.1))) { beta.1 <- initial.beta
                                    scale <- initial.scale
                                     break
        }
        if( (conv==1) )
        {
            # check for convergence
            if( norm_fastS( beta - beta.1 ) / norm_fastS(beta) < 1e-20 ) break
            # magic number alert!!!
        }
        res <- y - x %*% beta.1
        beta <- beta.1
    }
    
    res <- y - x %*% beta
    # get the residuals from the last beta
    return(list(beta.rw = beta.1, scale.rw = scale))
    
}
scale1 <- function(u, b, cc, initial.sc=median(abs(u))/.6745) {
        # find the scale, full iterations
        max.it <- 200
        # magic number alert
        #sc <- median(abs(u))/.6745
        sc <- initial.sc
        i <- 0
        eps <- 1e-20
        # magic number alert
        err <- 1
        while( ( (i <- i+1) < max.it ) && (err > eps) ) {
            sc2 <- sqrt( sc^2 * mean( rho_fastS( u / sc, cc ) ) / b   )
            err <- abs(sc2/sc - 1)
            sc <- sc2
        }
        return(sc)
}
fast.s <- function(x, y, int=1, H1,k=2,best.r=5, b=.5, cc=1.56, seed,conv=1) {
	N<-1
    #
    # int = 1 -> add a column of ones to x
    # N = cant de sub-samples
    # k = number of refining iterations in ea. 
    #     subsample
    # k = 0 means "raw-subsampling"
    # b = right hand side of the equation
    # cc = corresponding tuning constant
    # best.r = number of "best betas" to remember
    #          from the subsamples. These will be later
    #          iterated until convergence


    # conv = 0 means "do k steps and don't check for convergence"
    # conv = 1 means "stop when convergence is detected, or the
    #                 maximum number of iterations is achieved"

    # the objective function, we solve loss.S(u,s,cc)=b for "s"
 




    
    n <- dim(x)[1]
    p <- dim(x)[2]
    
    if( int == 1) {
            x <- cbind(rep(1,n), x)
            p <- p + 1
    }

    if(!missing(seed)) set.seed(seed)

    best.betas <- matrix(0, best.r, p)
    best.scales <- rep(1e20, best.r)
    s.worst <- 1e20
    n.ref <- 1
    

        #
        # get a subsample
        #

	xs <- x[H1,]
	ys <- y[H1]
	beta <- our.solve(xs,ys)

        if(k>0) { 
            # do the refining
            tmp <- re.s(x=x,y=y,initial.beta=beta,k=k,conv=conv,b=b,cc=cc)
            beta.rw <- tmp$beta.rw
            scale.rw <- tmp$scale.rw
            res.rw <- y - x %*% beta.rw
        } else { #k = 0 means "no refining"
            beta.rw <- beta
            res.rw <- y - x %*% beta.rw
            scale.rw <- median(abs(res.rw))/.6745
        }


                best.scales[best.r]  <- scale1(res.rw,b,cc,scale.rw)
                best.betas[best.r,] <- beta.rw

    # do the complete refining step until convergence (conv=1) starting
    # from the best subsampling candidate (possibly refined)
    super.best.scale <- Inf
    # magic number alert
    for(i in best.r:1) {
        tmp <- re.s(x=x,y=y,initial.beta=best.betas[i,],
                initial.scale=best.scales[i],k=k,conv=conv,b=b,cc=cc)
        if(tmp$scale.rw < super.best.scale) {
            super.best.scale <- tmp$scale.rw
            super.best.beta <- tmp$beta.rw
        }   
    }
    return(list(as.vector(super.best.beta), super.best.scale))

}
