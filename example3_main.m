% Solve H-equation using ETD2RK 
% Written by Divya Jagannathan

clear all, clc, close all

% Physical parameters
alp = .33; gmma= 1; omg = 5; Tfinal = 5; y0 = 1;

npow = 10; ipow=1;
err=linspace(0,0,npow); h = linspace(0,0,npow);

while ipow<=npow
% iterative parameters
dt = 1/2^ipow; nsteps = ceil(Tfinal/dt)+1;

% Scheme related quantities
N = 101;
L = Lmat(N); [expLh,M1,M2] = expMat(-L,dt); f = fvec(N,gmma);

% this one's only for the position y...
somen = 0:1:(N-1); alter = (-1).^somen; basis4y = sqrt(2*pi)*alter'.*hermiteFunction(N,0);

% initial condition
c = y0*f;
y = linspace(0,0,nsteps);  y(1)=y0; t = linspace(0,0,nsteps);
ye = linspace(0,0,nsteps); ye(1)=y0;

% k - vector, t - scalar
dc11 = @(k,x) (2/pi)*exp(-k.^2*x).*((gmma*k.^2)./((k.^2-alp).^2+gmma^2*k.^2));
f1 = @(k,x) (2/pi)*((gmma*k.^2)./((k.^2-alp).^2+gmma^2*k.^2)).*((omg*exp(-x*(k.^2))-omg*cos(omg*x)+(k.^2).*sin(omg*x))./(omg^2+k.^4));

istep = 1; 
while istep<=nsteps
    Nforce0 = forcing(y(istep),t(istep),alp,omg);
    avec = expLh*c + M1*f*Nforce0;
    yhalf = dot(avec,basis4y);
    Nforce1 = forcing(yhalf,t(istep)+dt,alp,omg);
    c = avec + (M2/dt)*f*(Nforce1-Nforce0);
    y(istep+1) = dot(c,basis4y);
    t(istep+1) = t(istep)+dt;
    % Exact solution from Section 4.13 in Prasath et al (2019)
    ye(istep+1) = integral(@(k) dc11(k,t(istep+1)),0,Inf)*y0 + integral(@(k) f1(k,t(istep+1)),0,Inf);
    istep = istep+1;
end
err(ipow)  = sqrt(dt)*norm(y-ye,2);
h(ipow) = dt;
ipow=ipow+1;
end

figure(1)
hop = 10;
p1=plot(t,ye,'b-','LineWidth',2); hold on
p2=plot(t(1:hop:end),y(1:hop:end),'ko','MarkerSize',3.5,'MarkerFaceColor','k'); 
xlabel('time'); ylabel('solution');
legend([p1,p2],{'Exact','Numerical'});
title('dt=', h(end));
pbaspect([6,4,1])

figure(2)
loglog(h,err,'bo','MarkerFaceColor','b','MarkerSize',3); hold on;
grid on; xlabel('$\Delta t$','Interpreter','latex'); ylabel('$l_2$ error','Interpreter','latex'); axis square;


