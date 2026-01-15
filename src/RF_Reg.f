C       Last change:  LB   13 Mar 2002   12:01 pm
	program RRforest
c	copyright 1999 by leo Breiman
c	this is free software and can be used for any purpose. 
c	It comes with no guarantee.  
	


	parameter(mdim=12,nsample=506,ntest=1,nthsize=5,
     1	nrnodes=2*(nsample/nthsize)+1,jprint=1,jbt=500,mtry=mdim/3,
     1	imp=0,nimp=imax0(imp*nsample,1),mimp=imax0(imp*mdim,1),
     1  th=.5,impprint=0)
     
c	mdim=number of variables in data set

c	nsample=number of cases

c	number of samples in the test set.  ntest=1 implies no test set

c	nthsize=number of cases in a node below which the tree will not split,
c	setting nthsize=5 generally gives good results.

c	jprint=1 turns on output as trees are being added giving tree number,
c	test set mean-sum-of-squared prediction errors,
c	bagged mean-sum-of-squared prediction errors, and the bagged error 
c	divided by var(y).The bagged error estimates are known to be are 
c 	accurate estimates of the true error values.

c	jbt=number of trees in run.  200-500 gives pretty good results

c	mtry=number of variables to pick to split on at each node.  mdim/3
c		seems to give genrally good performance, but it can be 
c		altered uo or down

c	imp=1 turns on variable importance.  This is computed for the
c	mth variable as the percent rise in the test set mean sum-of-
c	squared errors when the mth variable is randomly permuted.

c	if impprint =1 then th=.5 says print out all variables with 
c	importance greater than .5 as the run is progressing.  
c	To keep from having a lot of output, set th high. 

c	at the end of the run all variable importances >0 will be printed out.

c	for each data set, the number of categoricals in each variable has to be
c	specified.  This is done by filling in the cat(mdim) vector See below.
C	cat(m)=1 says that variable m is numerical.  The program expects to see
c	categoricals specified by integers consecutively upward from 1.  


c	MAIN PROGRAM

	real x(mdim,nsample),y(nsample),yb(nsample),rsnodecost(nrnodes),
     1	rstreecost(0:nrnodes),bestcrit(nrnodes),sd(mdim),wts(nsample),
     1	v(nsample),ut(nsample),xt(nsample),xb(mdim,nsample),
     1	errimp(mimp), ytr(nsample),yptr(nsample),yl(nsample),xa(3*mdim),
     1	avnode(nrnodes),utr(nsample),predimp(nimp,mimp),za(mdim),
     1  tgini(mdim),upper(nrnodes),
     1	ypred(ntest),ytree(ntest),xts(mdim,ntest),yts(ntest)
	
	integer jdex(nsample),treemap(2,nrnodes), nodestatus(nrnodes),
     1	nodepop(nrnodes),npert(nsample),ip(mdim),nterm(0:nrnodes),
     1	nperm(nsample),parent(nrnodes),cat(mdim),nout(nsample),
     1	jin(nsample),isort(nsample),nodestart(nrnodes),ncase(nsample),
     1	nbrterm(nrnodes),jperm(jbt),mbest(nrnodes),incl(mdim)
	
c	initialize rn generator

	do n=1,717
	zz=rand(1)
	end do
	
c	declare number of values for each catgegorical

	do m=1,mdim
	cat(m)=1
	end do
	
	

c	BEGIN MAIN

	
c	read in data==example follows
	
c	open(7,file='data5',status='old')  
c	do k=1,500
c	read(7,*) (x(j,k),j=1,mdim),y(k)
c	end do
c	close(7)
c	if there is a test set
c	open(7,file='data5test',status='old')  
c	do k=1,300
c	read(7,*) (xts(j,k),j=1,mdim),yts(k)
c	end do
c	close(7)


 	open(7,file='boshouse',status='old')
	do k=1,nsample
	read(7,*) y(k),(x(j,k),j=1,mdim)
	end do
	close(7)
	
	
	qverrts=0
	averrb=0
	
		avy=0
		vary=0
		do n=1,nsample
		ntrue=n-1
		vary=vary+ntrue*(y(n)-avy)**2/(ntrue+1)
		avy=(ntrue*avy+y(n))/(ntrue+1)
		end do
		vary=vary/nsample
	
	
	call zervr(yptr,nsample)
	call zerv(nout,nsample)
	astr=0
	asd=0

	
	do jb=1,jbt
        call zerv(jin,nsample)
	do n=1,nsample
	k=int(rand1(1)*nsample)+1
	jin(k)=1
	yb(n)=y(k)
	do m=1,mdim
	xb(m,n)=x(m,k)
	end do
	end do
	
	
	
	nls=nsample
	
	call buildtree(xb,yb,yl,mdim,nls,nsample,treemap,jdex,
     1  upper,avnode,bestcrit, nodestatus,nodepop,nodestart,
     1  nrnodes,nthsize,rsnodecost,ncase,parent,ut,v,iprint,
     1  xt,mtry,ip,nlim,mbest,cat,tgini)

	ndbigtree=nrnodes
	do k=nrnodes,1,-1
	if (nodestatus(k).eq.0) ndbigtree=ndbigtree-1
	if (nodestatus(k).eq.2) nodestatus(k)=-1
	end do

	do n=1,50
	!write(*,*) n,avnode(n)
	end do
	!stop
	
	call zervr(ytr,nsample)
	
	call testreebag(x,nsample,mdim,treemap,nodestatus,
     1	nrnodes,ndbigtree,ytr,upper,avnode,mbest,cat)
     
     	
     
    	errb=0
	jout=0
	do n=1,nsample
	if(jin(n).eq.0) then
	yptr(n)=(nout(n)*yptr(n)+ytr(n))/(nout(n)+1)
	nout(n)=nout(n)+1
	end if
	
	
	if(nout(n).gt.0) jout=jout+1
	errb=errb+(y(n)-yptr(n))**2
	end do
	errb=errb/nsample
	
	call zervr(ytree,ntest)
	
	call testreebag(xts,ntest,mdim,treemap,nodestatus,
     1	nrnodes,ndbigtree,ytree,upper,avnode,mbest,cat)
    
	
	errts=0
	do n=1,ntest
	ypred(n)=((jb-1)*ypred(n)+ytree(n))/jb
	errts=errts+(yts(n)-ypred(n))**2
	end do
	errts=errts/ntest
	
	if(jprint.eq.1) then
	if(ntest.gt.1) then
	write(*,'(i5,3f11.3)')jb,errts,errb,errb/vary
	else
	write(*,'(i5,2f11.3)')jb,errb,errb/vary
	end if
	end if
	
	if(imp.eq.1) then
	do mr=1,mdim
	call permobmr(mr,x,utr,xt,jin,nsample,mdim)
	call testreebag(x,nsample,mdim,treemap,nodestatus,
     1	nrnodes,ndbigtree,ytr,upper,avnode,mbest,cat)
	do n=1,nsample
	x(mr,n)=xt(n)
	end do
	em=0
	do n=1,nsample
	if(jin(n).eq.0) then
	predimp(n,mr)=(nout(n)*predimp(n,mr)+ytr(n))
     1  /(nout(n)+1)
	end if
	em=em+(y(n)-predimp(n,mr))**2
	end do
	errimp(mr)=em/nsample
	if(impprint.eq.1) then 
	if(100*(errimp(mr)-errb)/errb.gt.th) 
     1	write(*,*) jb,mr,100*(errimp(mr)-errb)/errb
     	end if
     	end do 
	end if
	
	end do
c	end of tree iterations=======================================

	if(ntest.gt.1) then
	write(*,'(i5,3f11.3)')jbt,errts,errb,errb/vary
	else
	write(*,'(i5,2f11.3)')jbt,errb,errb/vary
	end if
	
	
	if(imp.eq.1) then
	write(*,*) "importances"
	do m=1,mdim
	errimp(m)=100*((errimp(m)/errb)-1)
	if(errimp(m).le.0) errimp(m)=0
	write(*,*) m,errimp(m)
	end do
	end if
	
	do m=1,mdim
c	write(*,*) m,tgini(m)
	end do
	
	
	end
	
	
c	END MAIN
	

c	SUBROUTINE BUILDTREE
	
	subroutine buildtree(x,y,yl,mdim,nls,nsample,treemap,
     1  jdex,upper,avnode,bestcrit, nodestatus,
     1  nodepop,nodestart,nrnodes,nthsize,rsnodecost,
     1  ncase,parent,ut,v,iprint,xt,mtry,ip,nlim,
     1  mbest,cat,tgini)
	
	integer treemap(2,nrnodes),parent(nrnodes),
     1  nodestatus(nrnodes),ip(mdim),nodepop(nrnodes),
     1  nodestart(nrnodes),jdex(nsample),ncase(nsample),
     1  mbest(nrnodes),cat(mdim)
	
	real  	y(nsample),bestcrit(nrnodes),x(mdim,nsample),
     1  avnode(nrnodes),xt(nsample),upper(nrnodes),
     1  v(nsample),ut(nsample),rsnodecost(nrnodes),
     1  yl(nsample),tgini(mdim)
	
	
	call zerv(nodestatus,nrnodes)
	call zerv(nodestart,nrnodes)
	call zerv(nodepop,nrnodes)
	call zervr(avnode,nrnodes)
	
	
	
	do n=1,nsample
	ut(n)=0
	jdex(n)=n
	end do

	ncur=1
	nodestart(1)=1
	nodepop(1)=nls
	nodestatus(1)=2
	
	av=0
	ss=0
	do n=1,nls
	d=y(jdex(n))
	ss=ss+(n-1)*(av-d)*(av-d)/n
	av=((n-1)*av+d)/n
	end do !n
	avnode(1)=av
	rsnodecost(1)=ss/nls
	
c	start main loop

	do 30 kbuild=1,nrnodes

	if (kbuild.gt.ncur) goto 50
	if (nodestatus(kbuild).ne.2) goto 30
	
c		initialize for next call to findbestsplit

	ndstart=nodestart(kbuild)
	ndend=ndstart+nodepop(kbuild)-1
	nodecnt=nodepop(kbuild)
	sumnode=nodecnt*avnode(kbuild)
	jstat=0
	
	
	call findbestsplit(x,xt,ut,jdex,y,mdim,nsample,
     1  ndstart,ndend,msplit,decsplit,ubest,ncase,ndendl,
     1  jstat,v,mtry,ip,nlim,sumnode,nodecnt,yl,cat)
	
	
	iprint=0
	if(iprint.eq.1) then
	write(*,*) kbuild,decsplit,ubest
	write(*,*) ndstart,ndendl,ndend
	write(*,*) ""
	end if
	
		
	if (jstat.eq.1) then
	nodestatus(kbuild)=-1
	go to 30
	else
		mbest(kbuild)=msplit
		upper(kbuild)=ubest
		bestcrit(kbuild)=decsplit
	end if
	
	tgini(msplit)=tgini(msplit)+decsplit
	
	
	
	
c	leftnode no.= ncur+1, rightnode no. = ncur+2.

	nodepop(ncur+1)=ndendl-ndstart+1
	nodepop(ncur+2)=ndend-ndendl
	nodestart(ncur+1)=ndstart
	nodestart(ncur+2)=ndendl+1
	
	av=0
	ss=0
	do n=ndstart,ndendl
	d=y(jdex(n))
	k=n-ndstart
	ss=ss+k*(av-d)*(av-d)/(k+1)
	av=(k*av+d)/(k+1)
	end do !n
	avnode(ncur+1)=av
	rsnodecost(ncur+1)=ss/nls
	
	av=0
	ss=0
	do n=ndendl+1,ndend
	d=y(jdex(n))
	k=n-ndendl-1
	ss=ss+k*(av-d)*(av-d)/(k+1)
	av=(k*av+d)/(k+1)
	end do !n
	avnode(ncur+2)=av
	rsnodecost(ncur+2)=ss/nls


	

c	check on nodestatus

	nodestatus(ncur+1)=2
	nodestatus(ncur+2)=2
	if (nodepop(ncur+1).le.nthsize)
     1  nodestatus(ncur+1)=-1
	if (nodepop(ncur+2).le.nthsize)
     1  nodestatus(ncur+2)=-1

	treemap(1,kbuild)=ncur+1
	treemap(2,kbuild)=ncur+2
	parent(ncur+1)=kbuild
	parent(ncur+2)=kbuild
	nodestatus(kbuild)=1
	ncur=ncur+2
	
	if (ncur.ge.nrnodes) goto 50
	
30	continue
50	continue

	end


c	SUBROUTINE FINDBESTSPLIT

	subroutine findbestsplit(x,xt,ut,jdex,y,mdim,
     1  nsample,ndstart,ndend,msplit,decsplit,ubest,
     1  ncase,ndendl,jstat,v,mtry,ip,nlim,
     1	sumnode,nodecnt,yl,cat)
	
	integer ncase(nsample),jdex(nsample),ip(mdim),
     1  ncat(32),icat(32),cat(mdim)
		
	real x(mdim,nsample),ut(nsample),xt(nsample),
     1  v(nsample),y(nsample),yl(nsample),
     1  sumcat(32),avcat(32),tavcat(32)
	
	
	
	
c		START BIG LOOP
		
		critmax=0
200		call zerv(ip,mdim)
		non=0
		do mt=1,mtry
		critvar=0
100		kv=int(rand(1)*mdim)+1
		if(ip(kv).eq.1) goto 100
		ip(kv)=100
		lc=cat(kv)
		if(lc.eq.1) then
		do n=ndstart,ndend
		xt(n)=x(kv,jdex(n))
		yl(n)=y(jdex(n))
		end do
		
		else
		
		call zervr(sumcat,32)
		call zerv(ncat,32)
		do n=ndstart,ndend
		l=nint(x(kv,jdex(n)))
		d=y(jdex(n))
		sumcat(l)=sumcat(l)+d
		ncat(l)=ncat(l)+1
		end do
		do  j=1,lc
		if(ncat(j).gt.0) then
		avcat(j)=sumcat(j)/ncat(j)
		else
		avcat(j)=0
		end if
		end do
		do n=1,nsample
		xt(n)=avcat(nint(x(kv,jdex(n))))
		yl(n)=y(jdex(n))
		end do
		end if
		
		
		
		 
		do n=ndstart,ndend
		v(n)=xt(n)
		end do
		do n=1,nsample
		ncase(n)=n
		end do
		call quicksort (v,ncase,ndstart,ndend,nsample)
		if(v(ndstart).ge.v(ndend))then
		non=non+1
		if(non.ge.3*mdim) then
		jstat=1
		return
		end if
		goto 100
		end if
		
		 
c		ncase(n)=case number of v nth from bottom

		suml=0
		sumr=sumnode
		npopl=0
		npopr=nodecnt
		
		
		do nsp=ndstart,ndend-1
		d=yl(ncase(nsp))
		suml=suml+d
		sumr=sumr-d
		npopl=npopl+1
		npopr=npopr-1
		if (v(nsp).lt.v(nsp+1)) then
		crit=(suml*suml/npopl)+(sumr*sumr/npopr)
		if (crit.gt.critvar) then
			ubestt=(v(nsp)+v(nsp+1))/2.0
			critvar=crit
			nbestt=nsp
			endif
		end if
		end do

                if(critvar.gt.critmax) then
                ubest=ubestt
                nbest=nbestt
                msplit=kv
		critmax=critvar
		do n=ndstart,ndend
		ut(n)=xt(n)
		end do
		if (cat(kv).gt.1) then
                ic=cat(kv)
                do j=1,ic
                tavcat(j)=avcat(j)
                end do
                end if
                end if
                end do

                nl=ndstart-1
		do nsp=ndstart,ndend
		if(ut(nsp).le.ubest) then
		nl=nl+1
		ncase(nl)=jdex(nsp)
		end if
		end do
		
                ndendl=max0(nl,ndstart+1)
		nr=ndendl
		do nsp=ndstart,ndend
		if(ut(nsp).gt.ubest) then
		nr=nr+1
		if(nr.gt.nsample) goto 765
		ncase(nr)=jdex(nsp)
		end if
		end do
765		continue

                if(ndendl.ge.ndend) ndendl=ndend-1


                do n=ndstart,ndend
		jdex(n)=ncase(n)
		end do
                lc=cat(msplit)
		if(lc.gt.1) then
		do j=1,lc
		if(tavcat(j).lt.ubest) then
		icat(j)=1
		else
		icat(j)=0
		end if
		end do
                call packlb(lc,icat,nubest)
		ubest=real(nubest)
		end if
		
		decsplit=critmax-(sumnode*sumnode/nodecnt)

 		end
		
	

	subroutine testreebag(x,nsample,mdim,treemap,nodestatus,
     1	nrnodes,ndbigtree,ytree,upper,avnode,mbest,cat)
	
	
	real x(mdim,nsample),
     1	upper(nrnodes),avnode(nrnodes),ytree(nsample)
	
	integer treemap(2,nrnodes),nodestatus(nrnodes),
     1	mbest(nrnodes),cat(mdim),icat(32)
     	
	
	do n=1,nsample
	kt=1
	do k=1,ndbigtree
	if(nodestatus(kt).eq.-1) then
	ytree(n)=avnode(kt)
	goto 100
	end if
	m=mbest(kt)
	lc=cat(m)
	if(lc.eq.1) then
		if (x(m,n).le.upper(kt)) then 
			kt=treemap(1,kt)
		else
			kt=treemap(2,kt)
		endif 
	else
	mm=nint(upper(kt))
	call unpacklb(lc,mm,icat)
	j=nint(x(m,n))
		if(icat(j).eq.1) then 
			kt=treemap(1,kt)
		else
			kt=treemap(2,kt)
		endif
	end if
	end do
100	continue
	end do

	end
	
	subroutine permobmr(mr,x,tp,tx,jin,nsample,mdim)
	dimension x(mdim,nsample), tp(nsample),tx(nsample),
     1  jin(nsample)
	kout=0
	call zervr(tp,nsample)
	do n=1,nsample
	if(jin(n).eq.0) then
	kout=kout+1
	tp(kout)=x(mr,n)
	end if
	end do !n
	call perm1(kout,nsample,tp)
	iout=0
	do n=1,nsample
	tx(n)=x(mr,n)
	if(jin(n).eq.0) then
	iout=iout+1
	x(mr,n)=tp(iout)
	end if
	end do
	end
	
	
	
	
	
	




c	MISCELLANOUS SMALL SUBROUTINES	

	
	subroutine zerv(ix,m1)
	integer ix(m1)
	do 10 n=1,m1
	ix(n)=0
10	continue
	end
	
	subroutine zervr(rx,m1)
	real rx(m1)
	do 10 n=1,m1
	rx(n)=0
10	continue
	end
	
	subroutine zerm(mx,m1,m2)
	integer mx(m1,m2)
	do 10 i=1,m1
	do 20 j=1,m2
	mx(i,j)=0
20	continue
10	continue
	end
	
	subroutine packlb(l,icat,npack)
	integer icat(32)
	
c	icat is a binary integer with ones for categories going left
c	and zeroes for those going right.  The sub returns npack- the integer 
c	corresponding whose binary expansion is icat.
	
	npack=0
	do 10,k=1,l
	npack=npack+icat(k)*(2**(k-1))
10	continue
	end
	
	subroutine unpacklb(l,npack,icat)
	
c	npack is a long integer.  The sub. returns icat, an integer of zeroes and
c	ones corresponding to the coefficients in the binary expansion of npack.
	
	integer icat(32)
	call zerv(icat,32)
	n=npack
	icat(1)=mod(n,2)
	do 10 k=2,l
	n=(n-icat(k-1))/2
	icat(k)=mod(n,2)
10	continue
	end
	
	
	
	subroutine quicksort (v,iperm,ii,jj,kk)
c
c     puts into iperm the permutation vector which sorts v into
c     increasing order.  only elementest from ii to jj are considered.
c     array iu(k) and array il(k) permit sorting up to 2**(k+1)-1 elements
c
c     this is a modification of acm algorithm #347 by r. c. singleton,
c     which is a modified hoare quicksort.
c
      dimension iperm(kk),v(kk),iu(32),il(32)
      integer t,tt
     
      m=1
      i=ii
      j=jj
 10   if (i.ge.j) go to 80
 20   k=i
      ij=(j+i)/2
      t=iperm(ij)
      vt=v(ij)
      if (v(i).le.vt) go to 30
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
 30   l=j
      if (v(j).ge.vt) go to 50
      iperm(ij)=iperm(j)
      iperm(j)=t
      t=iperm(ij)
      v(ij)=v(j)
      v(j)=vt
      vt=v(ij)
      if (v(i).le.vt) go to 50
      iperm(ij)=iperm(i)
      iperm(i)=t
      t=iperm(ij)
      v(ij)=v(i)
      v(i)=vt
      vt=v(ij)
      go to 50
 40   iperm(l)=iperm(k)
      iperm(k)=tt
      v(l)=v(k)
      v(k)=vtt
 50   l=l-1
      if (v(l).gt.vt) go to 50
      tt=iperm(l)
      vtt=v(l)
 60   k=k+1
      if (v(k).lt.vt) go to 60
      if (k.le.l) go to 40
      if (l-i.le.j-k) go to 70
      il(m)=i
      iu(m)=l
      i=k
      m=m+1
      go to 90
 70   il(m)=k
      iu(m)=j
      j=l
      m=m+1
      go to 90
 80   m=m-1
      if (m.eq.0) return
      i=il(m)
      j=iu(m)
 90   if (j-i.gt.10) go to 20
      if (i.eq.ii) go to 10
      i=i-1
 100  i=i+1
      if (i.eq.j) go to 80
      t=iperm(i+1)
      vt=v(i+1)
      if (v(i).le.vt) go to 100
      k=i
 110  iperm(k+1)=iperm(k)
      v(k+1)=v(k)
      k=k-1
      if (vt.lt.v(k)) go to 110
      iperm(k+1)=t
      v(k+1)=vt
      go to 100
      end
      
      	subroutine perm1(np,ns,tp)
	real tp(ns)
	j=np
11	rnd = rand(1)
	k=int(j*rnd)+1
	tx=tp(j)
	tp(j)=tp(k)
	tp(k)=tx
	j=j-1
	if(j.gt.1) go to 11
	end

	
	subroutine perm(ns,ntp)
	integer ntp(ns)
	do 1 n= 1,ns
	ntp(n)=n
1	continue        
	j=ns
11	rnd = rand1(1)
	k=int(j*rnd)+1
	jx=ntp(j)
	ntp(j)=ntp(k)
	ntp(k)=jx
	j=j-1
	if(j.gt.1) go to 11
	end
	
	
	real function rand(j)
	double precision dseed
	real u
	save dseed
	data dseed /17395/
	call lrnd(dseed,u)
	rand=u
	end

	subroutine lrnd(dseed,u)
	real u
	double precision dseed, d31m1
	data d31m1 /2147483647/
	dseed=dmod(16087*dseed,d31m1)
	u=dseed/d31m1
	return
	end
	
	real function rand1(j)
	double precision dseed1
	real u
	save dseed1
	data dseed1 /17395/
	call lrnd1(dseed1,u)
	rand1=u
	end

	subroutine lrnd1(dseed1,u)
	real u
	double precision dseed1, d31m1
	data d31m1 /2147483647/
	dseed1=dmod(16087*dseed1,d31m1)
	u=dseed1/d31m1
	return
	end
	
	
	
	function rnorm(j)
	double precision u,v
	u=dble(rand(1))
	v=dble(rand(1))
	rnorm=real(dsqrt(-2*dlog(u))*dcos(6.28318531*v))
	end 
	
	
	
	

	
	

 
	
	

	

	
	 
	
	


	

	

	

