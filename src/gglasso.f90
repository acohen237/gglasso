
! --------------------------------------------------
SUBROUTINE ls_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
  ! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: mnl
  INTEGER:: bn
  INTEGER::bs(bn)
  INTEGER::ix(bn)
  INTEGER::iy(bn)
  INTEGER:: nobs
  INTEGER::nvars
  INTEGER::dfmax
  INTEGER::pmax
  INTEGER::nlam
  INTEGER::nalam
  INTEGER::npass
  INTEGER::jerr
  INTEGER::maxit
  INTEGER::intr
  INTEGER:: idx(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION:: b0(nlam)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  DOUBLE PRECISION::d
  DOUBLE PRECISION::t
  DOUBLE PRECISION::dif
  DOUBLE PRECISION::unorm
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  INTEGER:: jx
  INTEGER:: jxx(bn)
  DOUBLE PRECISION:: ga(bn)
  DOUBLE PRECISION:: vl(nvars)
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars))
  ALLOCATE(oldbeta(0:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(oidx(1:bn))
  ! - - - checking pf - - -
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  jxx = 0
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  idx = 0
  oidx = 0
  npass = 0
  ni = npass
  alf = 0.0D0
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! flmin <1 means NOT user-supplied, so this is default
     flmin = Max (mfl, flmin)
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  vl = matmul(r, x)/nobs
  DO g = 1,bn
     ALLOCATE(u(bs(g)))
     u = vl(ix(g):iy(g))
     ga(g) = sqrt(dot_product(u,u))
     DEALLOCATE(u)
  ENDDO
  DO l=1,nlam
     al0 = al
     IF(flmin>=1.0D0) THEN
        al=ulam(l)
     ELSE
        IF(l > 2) THEN
           al=al*alf
        ELSE IF(l==1) THEN
           al=big
        ELSE IF(l==2) THEN
           al0 = 0.0D0
           DO g = 1,bn
              IF(pf(g)>0.0D0) THEN
                 al0 = max(al0, ga(g) / pf(g))
              ENDIF
           ENDDO
           al = al0 * alf
        ENDIF
     ENDIF
     tlam = (2.0*al-al0)
     DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
     ENDDO
     ! --------- outer loop ----------------------------
     DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
           DO j=1,ni
              g=idx(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           npass=npass+1
           dif=0.0D0
           DO g=1,bn
              IF(jxx(g) == 0) CYCLE
              startix=ix(g)
              endix=iy(g)
              ALLOCATE(u(bs(g)))
              ALLOCATE(dd(bs(g)))
              ALLOCATE(oldb(bs(g)))
              oldb=b(startix:endix)
              u=matmul(r,x(:,startix:endix))/nobs
              u=gam(g)*b(startix:endix)+u
              unorm=sqrt(dot_product(u,u))
              t=unorm-pf(g)*al
              IF(t>0.0D0) THEN
                 b(startix:endix)=u*t/(gam(g)*unorm)
              ELSE
                 b(startix:endix)=0.0D0
              ENDIF
              dd=b(startix:endix)-oldb
              IF(any(dd/=0.0D0)) THEN
                 dif=max(dif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(oidx(g)==0) THEN
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    oidx(g)=ni
                    idx(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(u,dd,oldb)
           ENDDO
           IF(intr /= 0) THEN
              d=sum(r)/nobs
              IF(d/=0.0D0) THEN
                 b(0)=b(0)+d
                 r=r-d
                 dif=max(dif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (dif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              npass=npass+1
              dif=0.0D0
              DO j=1,ni
                 g=idx(j)
                 startix=ix(g)
                 endix=iy(g)
                 ALLOCATE(u(bs(g)))
                 ALLOCATE(dd(bs(g)))
                 ALLOCATE(oldb(bs(g)))
                 oldb=b(startix:endix)
                 u=matmul(r,x(:,startix:endix))/nobs
                 u=gam(g)*b(startix:endix)+u
                 unorm=sqrt(dot_product(u,u))
                 t=unorm-pf(g)*al
                 IF(t>0.0D0) THEN
                    b(startix:endix)=u*t/(gam(g)*unorm)
                 ELSE
                    b(startix:endix)=0.0D0
                 ENDIF
                 dd=b(startix:endix)-oldb
                 IF(any(dd/=0.0D0)) THEN
                    dif=max(dif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(u,dd,oldb)
              ENDDO
              IF(intr /= 0) THEN
                 d=sum(r)/nobs
                 IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    dif=max(dif,d**2)
                 ENDIF
              ENDIF
              IF(dif<eps) EXIT
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! end inner loop
        ENDDO ! end middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(jxx(g) == 1) CYCLE
           ALLOCATE(u(bs(g)))
           u = vl(ix(g):iy(g))
           ga(g) = sqrt(dot_product(u,u))
           IF(ga(g) > al*pf(g))THEN
              jxx(g) = 1
              jx = 1
           ENDIF
           DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
     ENDDO ! End outer loop
     !---------- final update variable and save results------------
     IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
     ENDIF
     IF(ni>0) THEN
        DO j=1,ni
           g=idx(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     b0(l)=b(0)
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,oidx)
  RETURN
END SUBROUTINE ls_f

! --------------------------------------------------
SUBROUTINE ls_f_sparse (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
     eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr,alsparse)
  ! --------------------------------------------------
  
  IMPLICIT NONE
  ! - - - arg types - - -
  DOUBLE PRECISION, PARAMETER :: big=9.9E30
  DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
  INTEGER, PARAMETER :: mnlam = 6
  INTEGER:: mnl  
  INTEGER:: bn
  INTEGER::bs(bn)
  INTEGER::ix(bn)
  INTEGER::iy(bn)
  INTEGER:: nobs 
  INTEGER::nvars
  INTEGER::dfmax
  INTEGER::pmax
  INTEGER::nlam
  INTEGER::nalam
  INTEGER::npass
  INTEGER::jerr
  INTEGER::maxit
  INTEGER::intr
  INTEGER:: idx(pmax)
  INTEGER::nbeta(nlam)
  DOUBLE PRECISION:: flmin
  DOUBLE PRECISION::eps
  DOUBLE PRECISION:: x(nobs,nvars)
  DOUBLE PRECISION::y(nobs)
  DOUBLE PRECISION::pf(bn)
  DOUBLE PRECISION::ulam(nlam)
  DOUBLE PRECISION::gam(bn)
  DOUBLE PRECISION:: b0(nlam)
  DOUBLE PRECISION::beta(nvars,nlam)
  DOUBLE PRECISION::alam(nlam)
  DOUBLE PRECISION::alsparse
  ! - - - local declarations - - -
  DOUBLE PRECISION:: max_gam
  DOUBLE PRECISION::d
  ! DOUBLE PRECISION::t ! No longer using this
  DOUBLE PRECISION::dif
  ! DOUBLE PRECISION::unorm ! No longer using this for ls_new
  DOUBLE PRECISION::al
  DOUBLE PRECISION::alf
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r ! Residue y-beta_k*x etc
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u ! No longer using this for update step, but still need for other parts 
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
  INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
  INTEGER:: g
  INTEGER::j
  INTEGER::l
  INTEGER::ni
  INTEGER::me
  INTEGER::startix
  INTEGER::endix
  ! - - - Aaron's declarations
  DOUBLE PRECISION::snorm
  DOUBLE PRECISION::t_for_s(bn) ! this is for now just 1/gamma
  DOUBLE PRECISION::tea ! this takes the place of 't' in the update step for ls
  DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: s ! takes the place of 'u' in update for ls
  INTEGER::soft_g ! this is an iterating variable 'vectorizing' the soft thresholding operator
  !  INTEGER::al_t ! just for the lambda loop
  DOUBLE PRECISION::al2 ! just for the lambda loop
  INTEGER::vl_iter ! for iterating over columns(?) of x*r
  INTEGER::jk ! to break out of the while loop
  INTEGER::kill_count
  ! - - - begin - - -
  ! - - - local declarations - - -
  DOUBLE PRECISION:: tlam
  INTEGER:: jx
  INTEGER:: jxx(bn)
  DOUBLE PRECISION:: ga(bn) ! What is this for??
  DOUBLE PRECISION:: vl(nvars) ! What is this for?
  DOUBLE PRECISION:: al0
  ! - - - allocate variables - - -
  ALLOCATE(b(0:nvars))
  ALLOCATE(oldbeta(0:nvars))
  ALLOCATE(r(1:nobs))
  ALLOCATE(oidx(1:bn))
  !    ALLOCATE(al_sparse)
  ! - - - checking pf - ! pf is the relative penalties for each group
  IF(maxval(pf) <= 0.0D0) THEN
     jerr=10000
     RETURN
  ENDIF
  pf=max(0.0D0,pf)
  ! - - - some initial setup - - -
  jxx = 0 
  al = 0.0D0
  mnl = Min (mnlam, nlam)
  r = y
  b = 0.0D0
  oldbeta = 0.0D0
  idx = 0
  kill_count = 0
  oidx = 0
  npass = 0 ! This is a count, correct?
  ni = npass ! This controls so-called "outer loop"
  alf = 0.0D0
  ! al_sparse = 0.05 ! This is alpha for sparsity, controls sparse vs group; eventually should be an input parameter
  t_for_s = 1/gam ! might need to use a loop if no vectorization.........
  ! --------- lambda loop ----------------------------
  IF(flmin < 1.0D0) THEN ! THIS is the default...
     flmin = Max (mfl, flmin) ! just sets a threshold above zero 
     alf=flmin**(1.0D0/(nlam-1.0D0))
  ENDIF
  vl = matmul(r, x)/nobs ! Note r gets updated in middle and inner loop 
  DO g = 1,bn ! For each group...
     ALLOCATE(u(bs(g))) 
     u = vl(ix(g):iy(g))
     ga(g) = sqrt(dot_product(u,u)) 
     DEALLOCATE(u) 
  ENDDO
  al0 = 0.0D0
  DO vl_iter = 1,nvars
     al0 = max(al0, vl(vl_iter)) ! Infty norm of X'y, big overkill for lam_max
  ENDDO
  PRINT *, alsparse
  al2 = al0 ! / (1-alsparse) ! this value ensures all betas are 0 , divide by 1-a?
  jk = 1
  l = 0
  al = al2 ! start al at this big value
  DO WHILE (l < nlam) !! This is the start of the loop over all lambda values...
     ! print *, "l = ", l
     ! print *, "al = ", al
     IF(kill_count > 1000) RETURN
     kill_count = kill_count + 1
     al0 = al ! store old al value on subsequent loops, first set to al
     IF(flmin>=1.0D0) THEN ! user supplied lambda value, break out of everything
        l = l+1
        al=ulam(l)
        ! print *, "This is at the flmin step of the while loop"
     ELSE
        al=al*alf
        IF(l > 1) THEN ! have some active groups
           l = l+1
           ! print *, "This is the l>1 step of while loop"
        ELSE IF(l==0) THEN
           ! print *, "This is at l=0 step of while loop"
           ! Trying to find an active group
           jk = jk+1
           IF (jk > 50) RETURN ! here just to kill the while loop if we're stuck
        ENDIF
     ENDIF
     ! This is the start of the algorithm, for a given lambda...
     IF(l==0) THEN
        tlam = al
     ELSE
        tlam = max((2.0*al-al0), 0.0) ! Here is the strong rule...
     ENDIF
     print *, "Here is tlam = ", tlam
     DO g = 1, bn 
        IF(jxx(g) == 1)  CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1 ! Implementing the strong rule
     ENDDO
     ! --------- outer loop ---------------------------- ! 
     DO
        oldbeta(0)=b(0) 
        ! print *, "Here is the outer loop, and here's oldbeta:", oldbeta
        IF(ni>0) THEN
           DO j=1,ni
              g=idx(j)
              oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
           ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
           ! print *, "This is where we enter the middle loop"
           npass=npass+1 
           dif=0.0D0
           DO g=1,bn
              IF(jxx(g) == 0) THEN
                 ! print *, "in middle loop, jxx(g), for group (g) = ", g, ", is 0, so we CYCLE"
                 CYCLE
              ENDIF
              startix=ix(g)
              endix=iy(g)
              ALLOCATE(dd(bs(g)))
              ALLOCATE(oldb(bs(g)))
              oldb=b(startix:endix)
              ALLOCATE(s(bs(g)))
              s = matmul(r,x(:,startix:endix))/nobs
              s = s*t_for_s(g) + b(startix:endix)
              DO soft_g = 1, bs(g)
                 s(soft_g) = sign(abs(s(soft_g))-alsparse*t_for_s(g)*al, s(soft_g))
              ENDDO
              snorm = sqrt(dot_product(s,s))
              tea = snorm - t_for_s(g)*(1-alsparse)*al
              IF(tea>0.0D0) THEN
                 b(startix:endix) = s*tea/snorm
              ELSE
                 b(startix:endix) = 0.0D0
              ENDIF
              dd=b(startix:endix)-oldb
              IF(any(dd/=0.0D0)) THEN
                 dif=max(dif,gam(g)**2*dot_product(dd,dd))
                 r=r-matmul(x(:,startix:endix),dd)
                 IF(oidx(g)==0) THEN ! Here is where middle loop is different; if group g was not in oidx (active), and the
                    ! difference was nonzero, put it in active (ni)
                    ni=ni+1
                    IF(ni>pmax) EXIT
                    oidx(g)=ni
                    idx(ni)=g
                 ENDIF
              ENDIF
              DEALLOCATE(s,dd,oldb)
              ! DEALLOCATE(u,dd,oldb)
           ENDDO ! End middle loop
           IF(intr /= 0) THEN
              d=sum(r)/nobs
              IF(d/=0.0D0) THEN
                 b(0)=b(0)+d
                 r=r-d
                 dif=max(dif,d**2)
              ENDIF
           ENDIF
           IF (ni > pmax) EXIT
           IF (dif < eps) EXIT
           IF(npass > maxit) THEN
              jerr=-l
              RETURN
           ENDIF
           ! --inner loop----------------------
           DO
              ! PRINT *, "Here is where the inner loop starts"
              npass=npass+1
              dif=0.0D0
              DO j=1,ni
                 g=idx(j)
                 startix=ix(g)
                 endix=iy(g)
                 ALLOCATE(s(bs(g)))
                 ALLOCATE(dd(bs(g)))
                 ALLOCATE(oldb(bs(g)))
                 oldb=b(startix:endix)
                 s = matmul(r,x(:,startix:endix))/nobs
                 s = s*t_for_s(g) + b(startix:endix)
                 DO soft_g = 1, bs(g)
                    s(soft_g) = sign(abs(s(soft_g))-alsparse*t_for_s(g)*al, s(soft_g))
                 ENDDO
                 snorm = sqrt(dot_product(s,s))
                 tea = snorm - t_for_s(g)*(1-alsparse)*al
                 IF(tea>0.0D0) THEN
                    b(startix:endix) = s*tea/snorm
                 ELSE
                    b(startix:endix) = 0.0D0
                 ENDIF
                 dd=b(startix:endix)-oldb
                 IF(any(dd/=0.0D0)) THEN
                    dif=max(dif,gam(g)**2*dot_product(dd,dd))
                    r=r-matmul(x(:,startix:endix),dd)
                 ENDIF
                 DEALLOCATE(s,dd,oldb)
                 ! DEALLOCATE(u,dd,oldb)
              ENDDO ! END INNER LOOP
              IF(intr /= 0) THEN ! intr is whether to include intercept
                 d=sum(r)/nobs
                 IF(d/=0.0D0) THEN
                    b(0)=b(0)+d
                    r=r-d
                    dif=max(dif,d**2)
                 ENDIF
              ENDIF
              IF(dif<eps) EXIT ! Exit nearest loop. This is till convergence.
              IF(npass > maxit) THEN
                 jerr=-l
                 RETURN
              ENDIF
           ENDDO ! End Inner loop
        ENDDO ! End middle loop
        IF(ni>pmax) EXIT
        !--- final check ------------------------ ! This checks which violate KKT condition
        ! PRINT *, "Here is where the final check starts"
        jx = 0
        max_gam = maxval(gam)
        IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        vl = matmul(r, x)/nobs
        DO g = 1, bn
           IF(jxx(g) == 1) CYCLE
           ALLOCATE(u(bs(g)))
           u = vl(ix(g):iy(g))
           ga(g) = sqrt(dot_product(u,u))
           IF(ga(g) > al*pf(g))THEN
              jxx(g) = 1
              jx = 1
           ENDIF
           DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE ! This goes all the way back to outer loop
        EXIT
     ENDDO ! Ends outer loop
     !---------- final update variable and save results------------
     IF(l==0) THEN
        IF(maxval(jxx)==0) THEN
           CYCLE ! don't save anything, we're still decrementing lambda
        ELSE
           l=2
           alam(1) = al / alf ! store previous, larger value
        ENDIF
     ENDIF
     ! PRINT *, "Here is where the final update starts"
     IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
     ENDIF
     IF(ni>0) THEN
        DO j=1,ni
           g=idx(j)
           beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
     ENDIF
     nbeta(l)=ni
     b0(l)=b(0)
     alam(l)=al
     nalam=l
     IF (l < mnl) CYCLE
     me=0
     DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
     ENDDO
     IF(me>dfmax) EXIT
  ENDDO ! end lambda loop
  DEALLOCATE(b,oldbeta,r,oidx)
  RETURN
END SUBROUTINE ls_f_sparse

! --------------------------------------------------
SUBROUTINE log_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER:: mnl
    INTEGER:: bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER:: nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::maxit
    INTEGER::intr
    INTEGER:: idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::gam(bn)
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION:: max_gam
    DOUBLE PRECISION::d
    DOUBLE PRECISION::t
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER:: g
    INTEGER::j
    INTEGER::l
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - begin - - -
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars))
    ALLOCATE(oldbeta(0:nvars))
    ALLOCATE(r(1:nobs))
    ALLOCATE(oidx(1:bn))
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b=0.0D0
    oldbeta=0.0D0
    idx=0
    oidx=0
    npass=0
    ni=npass
    alf = 0.0D0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    vl = matmul(y/(1.0D0+exp(r)), x) / nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)))
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1, nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al0 = 0.0D0
                DO g = 1,bn
                    IF(pf(g)>0.0D0) THEN
                        al0 = max(al0, ga(g) / pf(g))
                    ENDIF
                END DO
                al = al0 * alf
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)))
                    ALLOCATE(dd(bs(g)))
                    ALLOCATE(oldb(bs(g)))
                    oldb=b(start:end)
                    u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(abs(dd)>0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r+y*matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                IF(intr /= 0) THEN
                    d = sum(y/(1.0D0+exp(r)))
                    d = 4.0D0*d/nobs
                    IF(d /= 0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif,d**2)
                    ENDIF
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)))
                        ALLOCATE(dd(bs(g)))
                        ALLOCATE(oldb(bs(g)))
                        oldb=b(start:end)
                        u=matmul(y/(1.0D0+exp(r)),x(:,start:end))/nobs
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(abs(dd)>0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r+y*matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    IF(intr /= 0) THEN
                        d = sum(y/(1.0D0+exp(r)))
                        d = 4.0D0*d/nobs
                        IF(d/=0.0D0) THEN
                            b(0)=b(0)+d
                            r=r+y*d
                            dif=max(dif, d**2)
                        ENDIF
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            vl = matmul(y/(1.0D0+exp(r)), x) / nobs
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)))
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g))THEN
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    DEALLOCATE(b,oldbeta,r,oidx)
    RETURN
END SUBROUTINE log_f

! --------------------------------------------------
SUBROUTINE hsvm_f (delta,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER:: mnl
    INTEGER:: bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER:: nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::maxit
    INTEGER::intr
    INTEGER::idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::delta
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::gam(bn)
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION:: max_gam
    DOUBLE PRECISION::d
    DOUBLE PRECISION::t
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::dl(nobs)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER::i
    INTEGER::g
    INTEGER::j
    INTEGER::l
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - begin - - -
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars))
    ALLOCATE(oldbeta(0:nvars))
    ALLOCATE(r(1:nobs))
    ALLOCATE(oidx(1:bn))
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf = max(0.0D0, pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b = 0.0D0
    oldbeta = 0.0D0
    idx = 0
    oidx = 0
    npass = 0
    ni = npass
    alf = 0.0D0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
    ENDIF
    vl = 0.0
    CALL hsvmdrv(delta,nobs,nvars,x,y,r,vl)
    DO g = 1,bn
            ALLOCATE(u(bs(g)))
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al0 = 0.0D0
                DO g = 1,bn
                    IF(pf(g)>0.0D0) THEN
                        al0 = max(al0, ga(g) / pf(g))
                    ENDIF
                END DO
                al = al0 * alf
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)))
                    ALLOCATE(dd(bs(g)))
                    ALLOCATE(oldb(bs(g)))
                    oldb=b(start:end)
                    u = 0.0D0
                    DO i = 1,nobs
                        IF (r(i) > 1.0D0) THEN
                            dl(i) = 0.0D0
                        ELSEIF (r(i) <= (1-delta)) THEN
                            dl(i) = 1.0D0
                        ELSE
                            dl(i) = (1.0D0 - r(i)) / delta
                        ENDIF
                        u = u + dl(i)*y(i)*x(i,start:end)/nobs
                    ENDDO
                    u=gam(g)*b(start:end) + u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(dd/=0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r+y*matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                IF(intr /= 0) THEN
                    d = 0.0D0
                    DO i = 1,nobs
                        IF (r(i) > 1.0D0) THEN
                            dl(i) = 0.0D0
                        ELSEIF (r(i) <= (1-delta)) THEN
                            dl(i) = 1.0D0
                        ELSE
                            dl(i) = (1.0D0 - r(i)) / delta
                        ENDIF
                        d = d + dl(i)*y(i)
                    ENDDO
                    d = 0.5 * delta * d / nobs
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif, d**2)
                    ENDIF
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)))
                        ALLOCATE(dd(bs(g)))
                        ALLOCATE(oldb(bs(g)))
                        oldb=b(start:end)
                        u = 0.0D0
                        DO i = 1,nobs
                            IF (r(i) > 1.0D0) THEN
                                dl(i) = 0.0D0
                            ELSEIF (r(i) <= (1-delta)) THEN
                                dl(i) = 1.0D0
                            ELSE
                                dl(i) = (1.0D0 - r(i)) / delta
                            ENDIF
                            u = u + dl(i)*y(i)*x(i,start:end)/nobs
                        ENDDO
                        u=gam(g)*b(start:end) + u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(dd/=0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r+y*matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    IF(intr /= 0) THEN
                        d = 0.0D0
                        DO i = 1,nobs
                            IF (r(i) > 1.0D0) THEN
                                dl(i) = 0.0D0
                            ELSEIF (r(i) <= (1-delta)) THEN
                                dl(i) = 1.0D0
                            ELSE
                                dl(i) = (1.0D0 - r(i)) / delta
                            ENDIF
                            d = d + dl(i)*y(i)
                        ENDDO
                        d = 0.5 * delta * d / nobs
                        IF(d/=0.0D0) THEN
                            b(0)=b(0)+d
                            r=r+y*d
                            dif=max(dif, d**2)
                        ENDIF
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                     jerr=-l
                     RETURN
                  ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            CALL hsvmdrv(delta,nobs,nvars,x,y,r,vl)
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)))
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g))THEN
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    DEALLOCATE(b,oldbeta,r,oidx)
    RETURN
END SUBROUTINE hsvm_f





SUBROUTINE hsvmdrv(delta,nobs,nvars,x,y,r,vl)
IMPLICIT NONE
INTEGER:: nobs
INTEGER:: nvars
INTEGER:: i
DOUBLE PRECISION:: delta
DOUBLE PRECISION:: dl(nobs)
DOUBLE PRECISION:: y(nobs)
DOUBLE PRECISION:: r(nobs)
DOUBLE PRECISION:: x(nobs,nvars)
DOUBLE PRECISION:: vl(nvars)
 vl = 0.0
 DO i = 1, nobs
     IF (r(i) > 1.0D0) THEN
        dl (i) = 0.0D0
     ELSEIF (r(i) <= (1-delta)) THEN
        dl (i) = 1.0D0
     ELSE
        dl (i) = (1.0D0 - r(i)) / delta
     ENDIF
 ENDDO
 vl = matmul(dl*y, x) / nobs
END SUBROUTINE hsvmdrv

! --------------------------------------------------
SUBROUTINE sqsvm_f (bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER:: mnl
    INTEGER:: bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER:: nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::maxit
    INTEGER::intr
    INTEGER:: idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::gam(bn)
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION:: max_gam
    DOUBLE PRECISION::d
    DOUBLE PRECISION::t
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::dl(nobs)
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER:: g
    INTEGER::j
    INTEGER::l
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - begin - - -
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars))
    ALLOCATE(oldbeta(0:nvars))
    ALLOCATE(r(1:nobs))
    ALLOCATE(oidx(1:bn))
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = 0.0D0
    b=0.0D0
    oldbeta=0.0D0
    idx=0
    oidx=0
    npass=0
    ni=npass
    alf = 0.0D0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    dl = 2.0D0 * dim (1.0D0, r)
    vl = matmul(dl*y, x) / nobs
    DO g = 1,bn
            ALLOCATE(u(bs(g)))
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al0 = 0.0D0
                DO g = 1,bn
                    IF(pf(g)>0.0D0) THEN
                        al0 = max(al0, ga(g) / pf(g))
                    ENDIF
                END DO
                al = al0 * alf
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)))
                    ALLOCATE(dd(bs(g)))
                    ALLOCATE(oldb(bs(g)))
                    oldb=b(start:end)
                    dl = 2.0D0 * dim(1.0D0, r)
                    u=matmul(y*dl,x(:,start:end))/nobs
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(abs(dd)>0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r+y*matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                IF(intr /= 0) THEN
                    dl = 2.0D0 * dim(1.0D0, r)
                    d = dot_product(y,dl)
                    d = 0.25*d/nobs
                    IF(d /= 0.0D0) THEN
                        b(0)=b(0)+d
                        r=r+y*d
                        dif=max(dif,d**2)
                    ENDIF
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
                ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)))
                        ALLOCATE(dd(bs(g)))
                        ALLOCATE(oldb(bs(g)))
                        oldb=b(start:end)
                        dl = 2.0D0 * dim(1.0D0, r)
                        u=matmul(y*dl,x(:,start:end))/nobs
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(abs(dd)>0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r+y*matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    IF(intr /= 0) THEN
                        dl = 2.0D0 * dim(1.0D0, r)
                        d = dot_product(y,dl)
                        d = 0.25*d/nobs
                        IF(d/=0.0D0) THEN
                            b(0)=b(0)+d
                            r=r+y*d
                            dif=max(dif,d**2)
                        ENDIF
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                     jerr=-l
                     RETURN
                  ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            dl = 2.0D0 * dim (1.0D0, r)
            vl = matmul(dl*y, x) / nobs
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)))
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g))THEN
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    DEALLOCATE(b,oldbeta,r,oidx)
    RETURN
END SUBROUTINE sqsvm_f


! --------------------------------------------------
SUBROUTINE wls_f (bn,bs,ix,iy,wrs,wx,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9E30
    DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
    INTEGER, PARAMETER :: mnlam = 6
    INTEGER::mnl
    INTEGER::bn
    INTEGER::bs(bn)
    INTEGER::ix(bn)
    INTEGER::iy(bn)
    INTEGER::nobs
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::maxit
    INTEGER::intr
    INTEGER::idx(pmax)
    INTEGER::nbeta(nlam)
    DOUBLE PRECISION:: flmin
    DOUBLE PRECISION:: eps
    DOUBLE PRECISION:: wrs(nobs)
    DOUBLE PRECISION:: wx(nobs,nvars)
    DOUBLE PRECISION:: x(nobs,nvars)
    DOUBLE PRECISION::y(nobs)
    DOUBLE PRECISION::pf(bn)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::gam(bn)
    DOUBLE PRECISION:: b0(nlam)
    DOUBLE PRECISION::beta(nvars,nlam)
    DOUBLE PRECISION::alam(nlam)
    ! - - - local declarations - - -
    DOUBLE PRECISION:: max_gam
    DOUBLE PRECISION::d
    DOUBLE PRECISION::t
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldbeta
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: r
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: oldb
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: u
    DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: dd
    INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
    INTEGER:: g
    INTEGER::j
    INTEGER::l
    INTEGER::ni
    INTEGER::me
    INTEGER::start
    INTEGER::end
    ! - - - begin - - -
    ! - - - local declarations - - -
    DOUBLE PRECISION:: tlam
    INTEGER:: jx
    INTEGER:: jxx(bn)
    DOUBLE PRECISION:: ga(bn)
    DOUBLE PRECISION:: vl(nvars)
    DOUBLE PRECISION:: al0
! - - - allocate variables - - -
    ALLOCATE(b(0:nvars))
    ALLOCATE(oldbeta(0:nvars))
    ALLOCATE(r(1:nobs))
    ALLOCATE(oidx(1:bn))
! - - - checking pf - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    jxx = 0
    al = 0.0D0
    mnl = Min (mnlam, nlam)
    r = y
    b = 0.0D0
    oldbeta = 0.0D0
    idx = 0
    oidx = 0
    npass = 0
    ni = npass
    alf = 0.0D0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    vl = matmul(r, wx)
    DO g = 1,bn
            ALLOCATE(u(bs(g)))
            u = vl(ix(g):iy(g))
            ga(g) = sqrt(dot_product(u,u))
            DEALLOCATE(u)
    END DO
    DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al0 = 0.0D0
                DO g = 1,bn
                    IF(pf(g)>0.0D0) THEN
                        al0 = max(al0, ga(g) / pf(g))
                    ENDIF
                END DO
                al = al0 * alf
            ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
            IF(jxx(g) == 1) CYCLE
            IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
! --------- outer loop ----------------------------
        DO
            oldbeta(0)=b(0)
            IF(ni>0) THEN
                DO j=1,ni
                    g=idx(j)
                    oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
                ENDDO
            ENDIF
! --middle loop-------------------------------------
            DO
                npass=npass+1
                dif=0.0D0
                DO g=1,bn
                    IF(jxx(g) == 0) CYCLE
                    start=ix(g)
                    end=iy(g)
                    ALLOCATE(u(bs(g)))
                    ALLOCATE(dd(bs(g)))
                    ALLOCATE(oldb(bs(g)))
                    oldb=b(start:end)
                    u=matmul(r,wx(:,start:end))
                    u=gam(g)*b(start:end)+u
                    unorm=sqrt(dot_product(u,u))
                    t=unorm-pf(g)*al
                    IF(t>0.0D0) THEN
                        b(start:end)=u*t/(gam(g)*unorm)
                    ELSE
                        b(start:end)=0.0D0
                    ENDIF
                    dd=b(start:end)-oldb
                    IF(any(dd/=0.0D0)) THEN
                        dif=max(dif,gam(g)**2*dot_product(dd,dd))
                        r=r-matmul(x(:,start:end),dd)
                        IF(oidx(g)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            oidx(g)=ni
                            idx(ni)=g
                        ENDIF
                    ENDIF
                    DEALLOCATE(u,dd,oldb)
                ENDDO
                IF(intr /= 0) THEN
                    d=dot_product(r,wrs)
                    IF(d/=0.0D0) THEN
                        b(0)=b(0)+d
                        r=r-d
                        dif=max(dif,d**2)
                    ENDIF
                ENDIF
                IF (ni > pmax) EXIT
                IF (dif < eps) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
               ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    dif=0.0D0
                    DO j=1,ni
                        g=idx(j)
                        start=ix(g)
                        end=iy(g)
                        ALLOCATE(u(bs(g)))
                        ALLOCATE(dd(bs(g)))
                        ALLOCATE(oldb(bs(g)))
                        oldb=b(start:end)
                        u=matmul(r,wx(:,start:end))
                        u=gam(g)*b(start:end)+u
                        unorm=sqrt(dot_product(u,u))
                        t=unorm-pf(g)*al
                        IF(t>0.0D0) THEN
                            b(start:end)=u*t/(gam(g)*unorm)
                        ELSE
                            b(start:end)=0.0D0
                        ENDIF
                        dd=b(start:end)-oldb
                        IF(any(dd/=0.0D0)) THEN
                            dif=max(dif,gam(g)**2*dot_product(dd,dd))
                            r=r-matmul(x(:,start:end),dd)
                        ENDIF
                        DEALLOCATE(u,dd,oldb)
                    ENDDO
                    IF(intr /= 0) THEN
                        d=dot_product(r,wrs)
                        IF(d/=0.0D0) THEN
                            b(0)=b(0)+d
                            r=r-d
                            dif=max(dif,d**2)
                        ENDIF
                    ENDIF
                    IF(dif<eps) EXIT
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- final check ------------------------
            jx = 0
            max_gam = maxval(gam)
            IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
            IF (jx /= 0) CYCLE
            vl = matmul(r, wx)
            DO g = 1, bn
                IF(jxx(g) == 1) CYCLE
                ALLOCATE(u(bs(g)))
                u = vl(ix(g):iy(g))
                ga(g) = sqrt(dot_product(u,u))
                IF(ga(g) > al*pf(g))THEN
                    jxx(g) = 1
                    jx = 1
                ENDIF
                DEALLOCATE(u)
            ENDDO
            IF(jx == 1) CYCLE
            EXIT
        ENDDO
!---------- final update variable and save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) THEN
            DO j=1,ni
                g=idx(j)
                beta(ix(g):iy(g),l)=b(ix(g):iy(g))
            ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
            g=idx(j)
            IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
    ENDDO
    DEALLOCATE(b,oldbeta,r,oidx)
    RETURN
END SUBROUTINE wls_f





