        program rkformr
!  This program integrates the full Maxey-Riley equation with the
!  standard BBH using Runge-Kutta (stage 2 & 4) algorithm using Eq. 4.1
!  in the paper.
!  Written by: Divya Jagannathan

	implicit real*8 (a-h,k-l,o-z)
        integer, parameter :: nk = 51
        integer, parameter :: nstage = 4
        integer, parameter :: n2pow = 3
        real, parameter :: pi=4.d0*datan(1.d0)
        real*8 :: startTime, endTime
        real*8 :: hx(nk+1), hy(nk+1), k(nk+1), jacobk(nk+1), w(nk+1),
     c        ckdummy(nk+1), qxi(nstage), frcxi(nstage),
     c        a(nstage), b(nstage), axy(nstage,nstage), d(nk+1,nstage)
        character*15 filename, fstatefile

!       user input: step-length (dt), alpha (alp), gamma (gmma), final
!       time (tfinal), initial condition (qx0)
        
        dt = 1.d0/2.d0**n2pow; tfinal = 5.d0
        alp = 33.d-2; gmma = 1.d0; freq = 5.d0; 
        qx0 = 1.d0

!       derived quantities
        gsqt = gmma*dsqrt(dt)
        nsteps = floor(tfinal/dt)+1

!       initialise everything to 0
        call initialise (hx,hy,qx,qy,a,b,axy,d,k,w,jacobk,nk,nstage)

!       get the quadrature weights (w) and nodes (k) to compute history
!       integral using clencurt. get the runge-kutta coefficients for 
!       time integration (f stands for fractional)

        call clencurt (gsqt,k,ckdummy,w,jacobk,nk)
       
        if (nstage.eq.2) then
          if (n2pow.le.9) then
            write (filename,'(A,I1.1,A)') 'traj2f_h', n2pow,'.dat' 
          else
            write (filename,'(A,I2.1,A)') 'traj2f_h', n2pow,'.dat'
          endif
          call rkstage2f (k,gsqt,a,axy,b,d,nk,nstage)
          open(unit=1, file=filename)
        else
          if (n2pow.le.9) then
             write (filename,'(A,I1.1,A)') 'traj4f_h', n2pow,'.dat'
          else
             write (filename,'(A,I2.1,A)') 'traj4f_h', n2pow,'.dat'
          endif
          call rkstage4f (k,gsqt,a,axy,b,d,nk,nstage)
          open(unit=1, file=filename)
        endif

!       initialise the state variables [w,h]
        qx = qx0
        hx = (2.d0/pi)*(gsqt/(gsqt**2 + k**2))*qx
        time = 0.d0

        write(1,*) time, qx
        call CPU_TIME (startTime)

!       time-integration begins here

        do istep=2,nsteps
         call stage_initialise (qxi, frcxi, nstage)
         qxi(1) = qx
         frcxi(1) = sin(freq*(time+a(1)*dt)) - alp*qxi(1)
!        construct the intermediate stages (stage 2 onwards)     
         if (nstage.ne.1) then
          do istage = 2, nstage
           call compute_hintegral(k,hx,a(istage),w,jacobk,hxintegral,nk)
           qxi(istage) = hxintegral + dt*sum(axy(istage,:)*frcxi)
           frcxi(istage) = sin(freq*(time+a(istage)*dt))-alp*qxi(istage)
          enddo
         endif
!        construct the final stage
         call compute_hintegral (k,hx,1.d0,w,jacobk,hintegral,nk)
         qx = hintegral + dt*(sum(b*frcxi))
         hx = (2.d0/pi)*dt*(gsqt/(gsqt**2 + k**2))*matmul(d,frcxi)+
     c        dexp(-k**2)*hx
         time = time + dt
         write(1,*) time, qx
        enddo        

!       time-integration ends here

        call CPU_TIME (endTime)
        close(1)
        stop
        end

!       -------------main program ends here ------------

!--------AHOY! SR1:Clenshaw-Curtis
!        (from Trefethen's book)
        subroutine clencurt (gmod,k,chebnodes,weights,jacobk,nk)
        implicit real*8 (a-h,k-l,o-z)
        real*8 k(nk+1), chebnodes(nk+1), weights(nk+1), jacobk(nk+1),
     c         v(nk-1), thta(nk+1)

        pi = 4.d0*datan(1.d0)
        do niter=1,nk+1        
                thta(niter) = (niter-1)*pi/nk
                chebnodes(niter) = cos(thta(niter))
        enddo
        weights = 0.d0
        v = 1.d0
        if (mod(nk,2).eq.0.d0) then
               weights(1)=1.d0/(nk**2-1)
               weights(nk+1) = weights(1)
               do m=1,nk/2-1
                  v = v - 2.d0*cos(2*m*thta(2:nk))/(4*m**2-1)
               enddo
               v = v - cos(nk*thta(2:nk))/(nk**2-1)
        else 
               weights(1) = 1.d0/nk**2
               weights(nk+1) = weights(1)
               do m=1,(nk-1)/2
                 v = v - 2.d0*cos(2*m*thta(2:nk))/(4*m**2-1)
               enddo
        endif        
        weights(2:nk) = 2.d0*v/nk
        L = dsqrt(gmod)
        k = L*(1.d0+chebnodes)/(1.d0-chebnodes)
        jacobk = L*(2.d0/(1.d0-chebnodes)**2)
        jacobk(1) = 0.0
        k(1) = k(2)*10.d0  
        return
        end

!-------AHOY! SR2a: RK2f constants
        subroutine rkstage2f (k, gmod, a, axy, b, d,nk,nstage)
        implicit real*8 (a-h, k-l, o-z) 
        integer, parameter :: nx = 52
        real*8 a(nstage),b(nstage),axy(nstage,nstage), d(nk+1,nstage), 
     c         k(nk+1), psi0(nk+1), psih(nk+1),t(nx+1),
     c         cx(nx+1), wx(nx+1),jdummy(nx+1), kdummy(nx+1)
   
        a(1) = 0.d0
        a(2) = 1.d0
        call clencurt (gmod,kdummy,cx,wx,jdummy,nx)
        t = (1.d0+cx)*5.d-1
        phi0 = dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2),wx)
        phih = dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2)*
     c        dsqrt(1.d0-t), wx)
        psi0 = 1.d0
        psih = 2.d0/3.d0
        do ik = 1,nk
          psi0(ik)=dot_product(5.d-1*dexp(-t*(k(ik))**2),wx)
          psih(ik)=dot_product(5.d-1*dexp(-t*(k(ik))**2)*
     c           (1.d0-t)**0.5,wx)
        enddo
        b = 0.d0        
        b(2) = phih/dsqrt(a(2))
        b(1) = phi0-b(2)
        axy = 0.d0
        axy(2,1) = a(2)*phi0
        d = 0.d0
        d(:,2) = psih/dsqrt(a(2))
        d(:,1) = psi0 - d(:,2)
        return
        end

!-------AHOY! S2b: RK4f constants
        subroutine rkstage4f (k,gmod,a,axy,b,d,nk,ns)
        implicit real*8 (a-h, k-l, o-z)
        integer, parameter :: nx = 52       
        real*8 a(ns),b(ns),axy(ns,ns), d(nk+1,ns), 
     c      k(nk+1), d1(nk+1),d2(nk+1),d3(nk+1),d4(nk+1),t(nx+1),
     c      cx(nx+1),wx(nx+1),dummyjx(nx+1),kdummy(nx+1),df(nk+1,ns),
     c      bmatdummy(ns,ns), bmat(ns,ns),bmatinv(ns,ns),bf(ns)

        call clencurt (gmod, kdummy, cx, wx, dummyjx,nx)
        t = (1.d0+cx)*5.d-1 
        a = (/0.d0, 25.d-2, 9.d-1, 1.d0/)
        bmat = 0.d0
        bmat(1,:) = 1.d0
        bmat(2,:) = dsqrt(a)
        bmat(3,:) = a
        bmat(4,:) = a**1.5 
        bmatdummy = bmat
        call matinverse (bmatdummy, bmatinv, ns)
        bf(1)=dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2),wx)
        bf(2)=dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2)*
     c        (1.d0-t)**0.5, wx)
        bf(3)=dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2)*
     c        (1.d0-t), wx)
        bf(4)=dot_product(5.d-1*erfc(gmod*dsqrt(t))*dexp(t*gmod**2)*
     c        (1.d0-t)**1.5, wx)
        b = matmul(bmatinv, bf) 
        d1 = 1.d0
        d2 = 2.d0/3.d0
        d3 = 5.d-1
        d4 = 4.d-1
        do n=1,nk
           d1(n) = dot_product(5.d-1*dexp(-k(n)**2*t), wx)
           d2(n) = dot_product(5.d-1*dexp(-k(n)**2*t)*(1.d0-t)**0.5, wx)
           d3(n) = dot_product(5.d-1*dexp(-k(n)**2*t)*(1.d0-t), wx)
           d4(n) = dot_product(5.d-1*dexp(-k(n)**2*t)*(1.d0-t)**1.5, wx)
        enddo
        df(:,1) = d1
        df(:,2) = d2
        df(:,3) = d3
        df(:,4) = d4
        d = transpose(matmul(bmatinv, transpose(df)))
        axy = 0.d0
        phi02=dot_product(0.5*erfc(gmod*dsqrt(a(2)*t))*
     c        dexp(gmod**2*a(2)*t), wx)
        phi03 = dot_product(0.5*erfc(gmod*dsqrt(a(3)*t))*
     c        dexp(gmod**2*a(3)*t), wx)
        phi04 = dot_product(0.5*erfc(gmod*dsqrt(a(4)*t))*
     c        dexp(gmod**2*a(4)*t), wx)
        phih2 = dot_product(0.5*erfc(gmod*dsqrt(a(2)*t))*
     c        dexp(gmod**2*a(2)*t)*(1.d0-t)**0.5, wx)
        phih3 = dot_product(0.5*erfc(gmod*dsqrt(a(3)*t))*
     c        dexp(gmod**2*a(3)*t)*(1.d0-t)**0.5, wx)
        phih4 = dot_product(0.5*erfc(gmod*dsqrt(a(4)*t))*
     c        dexp(gmod**2*t)*(1.d0-t)**0.5, wx)
        phi14 = dot_product(0.5*erfc(gmod*dsqrt(a(4)*t))*
     c        dexp(gmod**2*a(4)*t)*(1.d0-t), wx)
        axy(2,1) = a(2)*phi02
        axy(3,2) =(d(nk+1,3)*a(3)**1.5*phih3+d(nk+1,2)*a(2)**1.5*phih2)/
     c             (d(nk+1,3)*a(2)**0.5)
        axy(3,1) = a(3)*phi03 - axy(3,2)
        axy(4,3) = (a(2)**0.5*a(4)**1.5*phih4-a(4)**2*phi14)/
     c             ((a(2)*a(3))**0.5-a(3))
        axy(4,2) = (a(4)**2*phi14/a(2))-(a(3)/a(2))*axy(4,3)
        axy(4,1) = a(4)*phi04 - axy(4,2) - axy(4,3)
        return
        end

!-------AHOY! S3:initialise data
        subroutine initialise (hx,hy,qx,qy,a,b,axy,d,k,w,jacobk,N,ns)
        real*8 qx, qy, a(ns), b(ns), d(N+1,ns), axy(ns,ns), 
     c         hx(N+1), hy(N+1), k(N+1), w(N+1), jacobk(N+1)
        qx = 0.d0; qy = 0.d0
        hx = 0.d0; hy = 0.d0
        a = 0.d0; b=0.d0; axy = 0.d0; d=0.d0
        k = 0.d0; w=0.d0; jacobk=0.d0
        return
        end

!-------AHOY! S4: stage initialise
        subroutine stage_initialise (q,frc,ns)
        implicit real*8 (a-h,k-l,o-z)
        real*8 q(ns), frc(ns)
        q = 0.d0
        frc = 0.d0
        return 
        end

!-------AHOY! S5: compute history integral
        subroutine compute_hintegral (k,his,ai,w,jacobk,hisintegral,nk)
        implicit real*8 (a-h, k-l, o-z)
        real*8 k(nk+1), w(nk+1), jacobk(nk+1), his(nk+1)
        hisintegral=dot_product((dexp(-ai*k**2)*jacobk)*his, w)
        return
        end

!-------AHOY! S6: invert a matrix-LU
!       (from wikipedia/online source)  
        subroutine matinverse (a,c,n)
        integer n, i, j, k
        real*8 L(n,n), U(n,n), a(n,n), c(n,n),
     c         b(n), d(n), x(n), coeff
        L = 0.d0
        U = 0.d0
        b = 0.d0
!       S1:
        do k=1,n-1
           do i=k+1,n
                coeff=a(i,k)/a(k,k)
                L(i,k) = coeff
                do j=k+1,n
                   a(i,j) = a(i,j)-coeff*a(k,j)
                enddo
           enddo
        enddo
!       S2:        
        do i=1,n
           L(i,i) = 1.d0
        enddo
        do j=1,n
         do i=1,j
            U(i,j) = a(i,j)
         enddo
        enddo
!       S3:
        do k=1,n
           b(k) = 1.d0
           d(1) = b(1)
           do i=2,n
              d(i) = b(i)
              do j=1,i-1
                 d(i)=d(i)-L(i,j)*d(j)
              enddo
            enddo
            x(n)=d(n)/U(n,n)
            do i=n-1,1,-1
               x(i)=d(i)
               do j=n,i+1,-1
                  x(i) = x(i)-U(i,j)*x(j)
               enddo
               x(i)=x(i)/U(i,i)
            enddo
            do i=1,n
               c(i,k)=x(i)
            enddo
            b(k) = 0.d0
        enddo
        end
