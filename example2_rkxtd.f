        program rk4mr
!  This program integrates the full Maxey-Riley equation with the
!  standard BBH using Runge-Kutta algorithm for a particle in 2D
!  stationary Lamb-Oseen vortex. 
!  Written by: Divya Jagannathan

	implicit real*8 (a-h,k-l,o-z)
        integer, parameter :: nk = 51
        integer, parameter :: nstage = 4
        integer, parameter :: n2pow = 3
        real, parameter :: pi = 4.d0*datan(1.d0)
        real*8 :: startTime, endTime
        real*8 :: hx(nk+1), hy(nk+1), k(nk+1), jacobk(nk+1), w(nk+1),
     c        ckdummy(nk+1), qxi(nstage), qyi(nstage), frcxi(nstage),
     c        frcyi(nstage), xpi(nstage), ypi(nstage), uxi(nstage),
     c        uyi(nstage),a(nstage), b(nstage), axy(nstage,nstage), 
     c        d(nk+1,nstage),xp,yp,
     c        pb(nstage), paxy(nstage,nstage), flowatxy(8)
        character*15 filename

!       user input: step-length (dt), R (denR), Stokes# (St),
!       initial conditions (pos-x0,y0; slip vel-qx0,qy0, final time
!       (tfine)

        dt = 1.d0/2.d0**n2pow; tfinal = 2.d2
        denR = 3.d0; St = 1.d0
        qx0 = 1.d0; qy0 = 0.d0; x0 = 1.d0; y0 = 0.d0

!       derived quantities for the integrator
        alp = 1.d0/St
        gmma = dsqrt(3.d0*alp/denR)
        gsqt = gmma*dsqrt(dt)
        nsteps = floor(tfinal/dt)+1

!       initialise everything to 0
        call initialise (hx,hy,a,b,axy,d,pb,paxy,
     c       k,w,jacobk,nk,nstage)

!       get the quadrature weights (w) and nodes (k) to compute history
!       integral using clencurt. get the runge-kutta coefficients using 
!       rkstage2/4f (f stands for fractional)
        
        call clencurt (gsqt,k,ckdummy,w,jacobk,nk)

        if (nstage.eq.2) then
          if (n2pow.le.9) then
            write (filename,'(A,I1.1,A)') 'traj2f_h', n2pow,'.dat' 
          else
            write (filename,'(A,I2.1,A)') 'traj2f_h', n2pow,'.dat'
          endif
          call rkstage2f (k,gsqt,a,axy,b,d,paxy,pb,nk,nstage)
        else
          if (n2pow.le.9) then
             write (filename,'(A,I1.1,A)') 'traj4f_h', n2pow,'.dat'
          else
             write (filename,'(A,I2.1,A)') 'traj4f_h', n2pow,'.dat'
          endif
          call rkstage4f (k,gsqt,a,axy,b,d,paxy,pb,nk,nstage)
        endif 

!       initialise to non-zero values, if need be
        qx = qx0; qy = qy0
        xp = x0; yp = y0
        hx = (2.d0/pi)*(gsqt/(gsqt**2+k**2))*qx
        hy = (2.d0/pi)*(gsqt/(gsqt**2+k**2))*qy
        time = 0.d0

        open(unit=1,file=filename)

        call CPU_TIME (startTime)

!       time-integration begins here
        do istep=2,nsteps
          call stage_initialise (qxi,qyi,frcxi,frcyi,xpi,ypi,
     c         uxi,uyi,flowatxy,nstage)
          qxi(1) = qx; qyi(1) = qy
          xpi(1) = xp; ypi(1) = yp
          call lambOseen (flowatxy, xpi(1), ypi(1), 1.d0)
          uxi(1) = flowatxy(1); uyi(1) = flowatxy(2)
          call qFrcng (qxi(1),qyi(1),time+a(1)*dt,flowatxy,
     c         frcxi(1),frcyi(1),alp,denR)

!         File writes: 1. time, 2. x-slip vel, 3. y-slip vel, 4. x-particle
!         position, 5. y-particle position
          if (mod(istep,1).eq.0.d0) then
             write(1,*) time, qx, qy, xp, yp
          endif

!        construct the intermediate stages (stage 2 onwards)     
         do istage = 2, nstage
          call compute_hintegral (k,hx,hy,a(istage),w,jacobk,
     c         hxintegral,hyintegral,nk)
          qxi(istage)=hxintegral+dt*sum(axy(istage,:)*frcxi)
          qyi(istage)=hyintegral+dt*sum(axy(istage,:)*frcyi)
          xpi(istage)=xp+dt*sum(paxy(istage,:)*(qxi+uxi))
          ypi(istage)=yp+dt*sum(paxy(istage,:)*(qyi+uyi))
          call lambOseen (flowatxy, xpi(istage),ypi(istage),1.d0)
          uxi(istage) = flowatxy(1)
          uyi(istage) = flowatxy(2)
          call qFrcng (qxi(istage),qyi(istage),time+a(istage)*dt,
     c         flowatxy,frcxi(istage),frcyi(istage),alp,denR)
         enddo
         
         call compute_hintegral (k,hx,hy,1.d0,w,jacobk,
     c        hxintegral,hyintegral,nk)
         qx = hxintegral + dt*(sum(b*frcxi))
         qy = hyintegral + dt*(sum(b*frcyi))
         xp = xp + dt*(sum(pb*(qxi+uxi)))
         yp = yp + dt*(sum(pb*(qyi+uyi)))
         hx = (2.d0/pi)*dt*(gsqt/(gsqt**2+k**2))*matmul(d,frcxi)+
     c        dexp(-k**2)*hx
         hy = (2.d0/pi)*dt*(gsqt/(gsqt**2+k**2))*matmul(d,frcyi)+
     c        dexp(-k**2)*hy
         time = time + dt
        enddo        

        call CPU_TIME (endTime)
        write(*,*) "iterations, time: ", nsteps, endTime-startTime
        close(1)
        stop
        end


!--------AHOY! SR1:Clenshaw-Curtis
        subroutine clencurt (gmod,k, cnodes, weights,jacobk,nk)
        implicit real*8 (a-h,k-l,o-z)
        real*8 k(nk+1), cnodes(nk+1), weights(nk+1), jacobk(nk+1),
     c          v(nk-1), thta(nk+1)
        pi = 4.d0*datan(1.d0)
        do niter=1,nk+1        
                thta(niter) = (niter-1)*pi/nk
                cnodes(niter) = cos(thta(niter))
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
        k = L*(1.d0+cnodes)/(1.d0-cnodes)
        jacobk = L*(2.d0/(1.d0-cnodes)**2)
        jacobk(1) = 0.0
        k(1) = k(2)*10.d0  
        return
        end

!-------AHOY! SR2: RK2f constants
        subroutine rkstage2f (k,gmod,a,axy,b,d,paxy,pb,nk,nstage)
        implicit real*8 (a-h, k-l, o-z)
        integer, parameter :: nx = 57
        real*8 a(nstage),b(nstage),axy(nstage,nstage), d(nk+1,nstage), 
     c         pa(nstage),pb(nstage),paxy(nstage,nstage),
     c         k(nk+1), psi0(nk+1), psih(nk+1), t(nx+1),
     c         cx(nx+1), wx(nx+1), jdummy(nx+1), kdummy(nx+1)
   
        a = (/0.d0, 1.d0/)
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
     c             (1.d0-t)**0.5,wx)
        enddo
        b = 0.d0        
        b(2) = phih/dsqrt(a(2))
        b(1) = phi0-b(2)
        axy = 0.d0
        axy(2,1) = a(2)*phi0
        d(:,2) = psih/dsqrt(a(2))
        d(:,1) = psi0 - d(:,2)
!       for position
        paxy = 0.d0
        paxy(2,1)= 1.d0
        pb = (/5.d-1,5.d-1/)
        return
        end

!-------AHOY! S3:initialise data
        subroutine initialise (hx,hy,a,b,axy,d,pb,paxy,
     c         k,w,jacobk,N,ns)
        real*8 a(ns), b(ns), d(N+1,ns), axy(ns,ns), 
     c         pa(ns), pb(ns), paxy(ns,ns),
     c         hx(N+1), hy(N+1), k(N+1), w(N+1), jacobk(N+1)
        hx = 0.d0; hy = 0.d0
        a = 0.d0; b = 0.d0; axy = 0.d0; d = 0.d0
        pb = 0.d0; paxy = 0.d0
        k = 0.d0; w = 0.d0; jacobk = 0.d0
        return
        end

!-------AHOY! S4: stage initialise
        subroutine stage_initialise (qx,qy,frcx,frcy,xp,yp,u,v,fxy,ns)
        implicit real*8 (a-h,k-l,o-z)
        real*8 qx(ns),qy(ns),frcx(ns),frcy(ns),xp(ns),yp(ns),
     c         u(ns), v(ns), fxy(8)
        qx = 0.d0; qy=0.d0
        frcx = 0.d0; frcy = 0.d0
        xp = 0.d0; yp = 0.d0
        u = 0.d0; v = 0.d0
        fxy = 0.d0
        return 
        end

!-------AHOY! S5: compute history integral
        subroutine compute_hintegral (k,hx,hy,ai,w,jacobk,hxi,hyi,nk)
        implicit real*8 (a-h, k-l, o-z)
        real*8 k(nk+1), w(nk+1), jacobk(nk+1), hx(nk+1),hy(nk+1)
        hxi=dot_product((dexp(-ai*k**2)*jacobk)*hx, w)
        hyi=dot_product((dexp(-ai*k**2)*jacobk)*hy, w)
        return
        end

!-------AHOY! S6: invert a matrix-LU
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

!-------AHOY! S7: rk stage 4
        subroutine rkstage4f (k,gmod,a,axy,b,d,paxy,pb,nk,ns)
        implicit real*8 (a-h, k-l, o-z)
        integer, parameter :: nx = 57       
        real*8 a(ns),b(ns),axy(ns,ns), d(nk+1,ns), 
     c      pb(ns), paxy(ns,ns), pbf(ns),
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
        pbf(1)=dot_product(5.d-1*(1.d0-t)**0.d0,wx)
        pbf(2)=dot_product(5.d-1*(1.d0-t)**0.5,wx)
        pbf(3)=dot_product(5.d-1*(1.d0-t),wx)
        pbf(4)=dot_product(5.d-1*(1.d0-t)**1.5,wx)
        b = matmul(bmatinv, bf) 
        pb = matmul(bmatinv, pbf)
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
        axy(3,2) = (b(3)*a(3)**1.5*phih3+b(2)*a(2)**1.5*phih2)/
     c             (b(3)*a(2)**0.5)
        axy(3,1) = a(3)*phi03 - axy(3,2)
        axy(4,3) = (a(2)**0.5*a(4)**1.5*phih4-a(4)**2*phi14)/
     c             ((a(2)*a(3))**0.5-a(3))
        axy(4,2) = (a(4)**2*phi14/a(2))-(a(3)/a(2))*axy(4,3)
        axy(4,1) = a(4)*phi04 - axy(4,2) - axy(4,3)
!       for position
        paxy = 0.d0
        vone02=dot_product(5.d-1*(1.d0-t)**0.d0, wx)
        vone03 = dot_product(5.d-1*(1.d0-t)**0.d0, wx)
        vone04 = dot_product(5.d-1*(1.d0-t)**0.d0, wx)
        voneh2 = dot_product(5.d-1*(1.d0-t)**0.5, wx)
        voneh3 = dot_product(5.d-1*(1.d0-t)**0.5, wx)
        voneh4 = dot_product(5.d-1*(1.d0-t)**0.5, wx)
        vone14 = dot_product(5.d-1*(1.d0-t), wx)
        paxy(2,1) = a(2)*vone02
        paxy(3,2) = (pb(3)*a(3)**1.5*voneh3+pb(2)*a(2)**1.5*voneh2)/
     c             (pb(3)*a(2)**0.5)
        paxy(3,1) = a(3)*vone03 - paxy(3,2)
        paxy(4,3) = (a(2)**0.5*a(4)**1.5*voneh4-a(4)**2*vone14)/
     c             ((a(2)*a(3))**0.5-a(3))
        paxy(4,2) = (a(4)**2*vone14/a(2))-(a(3)/a(2))*paxy(4,3)
        paxy(4,1) = a(4)*vone04 - paxy(4,2) - paxy(4,3)
        return
        end

!       AHOY! S8: forcing for the slip vel eqn
        subroutine qFrcng (qx,qy,time,fatxy,frcx,frcy,alp,denR)
	implicit real*8 (a-h,k-l,o-z)
        real*8 frcx,frcy,time,qx,qy, fatxy(8)
        frcx=(1.d0/denR-1.d0)*(fatxy(3)+fatxy(1)*fatxy(5)+
     c       fatxy(2)*fatxy(6))-(alp*qx)-(qx*fatxy(5)+qy*fatxy(6))
        frcy=(1.d0/denR-1.d0)*(fatxy(4)+fatxy(1)*fatxy(7)+
     c       fatxy(2)*fatxy(8))-(alp*qy)-(qx*fatxy(7)+qy*fatxy(8))
        return
        end

!       AHOY! S9: Lamb Oseen flow field
        subroutine lambOseen (flowatxy,x,y,St)
        implicit real*8 (a-h,k-l,o-z)
        real*8 u, v, dudt, dvdt, dudx, dudy, dvdx, dvdy, flowatxy(8)
        rp = dsqrt(x**2+y**2)
        u = -(1.d0/St)*(y/rp**2)*(1.d0-dexp(-St*rp**2))
        v = (1.d0/St)*(x/rp**2)*(1.d0-dexp(-St*rp**2))
        dudt = 0.d0
        dvdt = 0.d0
        dudx = ((2.d0*x*y)/(St*rp**4))*(1.d0-dexp(-St*rp**2))-
     c         ((2.d0*x*y)/(rp**2))*dexp(-St*rp**2)
        dudy = -(1.d0/St)*(1.d0/rp**2 - (2.d0*y**2)/rp**4)*
     c       (1.d0-dexp(-St*rp**2))-((2.d0*y**2)/rp**2)*dexp(-St*rp**2)
        dvdx = (1.d0/St)*(1.d0/rp**2 - (2.d0*x**2)/rp**4)*
     c       (1.d0-dexp(-St*rp**2))+((2.d0*x**2)/rp**2)*dexp(-St*rp**2)
        dvdy = -1.d0*dudx
        flowatxy = (/u,v,dudt,dvdt,dudx,dudy,dvdx,dvdy/)
        return
        end

