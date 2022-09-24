% This program solves the 2D Schnakenberg reaction-diffusion system using
% the numerical method intriduced by Kristiansen, 2008, where we use
% Chebysehv interpolation to approximate the Laplacian terms, which
% converts the PDEs into a system of ODEs. The ODEs system is then solved
% by MATLAB's built-in solver `ode15s`. 
% Mathematical Biology modeling project, 2018 Michaelmas Term, Oxford
% University. 


clear all; 
global gamma ab diag10 u0u0 v0v0 u0v0
global N2 Nx Ny
global d Dx Dy
global a b

%% parameters 
% parameter set 1: spot patterns 
% Lx = pi/2; Ly = pi/2; a = 0.05; b = 1; d = 80; t = 200; gamma = 100; Nx = 15;
% parameter set 2: spot patterns (more spots on a larger domain)
% Lx = pi; Ly = pi; a = 0.1; b = 0.9; d = 40; t = 150; gamma = 100; Nx = 20;
% parameter set 3: stripe patterns 
 Lx = pi/2; Ly = pi/10; a = 0.05; b = 1; d = 80; t = 200; gamma = 80; Nx = 15;

tic 
tspan = [0 t];
Ny = floor(Ly/Lx*Nx);
if mod(Ny,2) ~= 0;
    Ny = Ny+1;
end
N2 = (Nx+1)*(Ny+1);
plotN = 100;

x = linspace(0,Lx,Nx+1);
y = linspace(0,Ly,Ny+1);

L = ((d*(b-a)-(a+b)^3)-sqrt((d*(b-a)-(a+b)^3)^2-...
    4*d*(a+b)^4))/(2*d*(a+b));
M = ((d*(b-a)-(a+b)^3)+sqrt((d*(b-a)-(a+b)^3)^2-...
    4*d*(a+b)^4))/(2*d*(a+b));


kvec = zeros(Nx+1,Ny+1);
for n= 0:Nx;
    for m = 0:Ny;
        kvec(n+1,m+1) = sqrt((n*pi)^2/Lx^2+(m*pi)^2/Ly^2);
    end
end
kvec = kvec(:);
kvec = sort(kvec(:));
[X,Y] = meshgrid(x,y);

%% cheb. approximation
[dx,vecx] = cheb(Nx);
dx = -(2/Lx)*dx;
[dy,vecy] = cheb(Ny); dy = -(2/Ly)*dy;
x = -Lx*vecx'/2+Lx/2; y = -Ly*vecy'/2+Ly/2;
I = eye(Ny+1);
dx = kron(dx,I);
I1 = [1 0;0 0]; I2 = [0 0;0 1];
Dx = sparse(kron(I1,dx)+kron(I2,dx)); I = eye(Nx+1);
dy = kron(I,dy);
Dy = sparse(kron(I1,dy)+kron(I2,dy));

% initial conditions
u0 = a+b; v0 = b/(a+b)^2; uold = zeros(Ny+1,Nx+1);
for i=1:Nx
    for j = 1:Ny
        m = max(i-1,j-1);
        c = rand(1)*(-1)^(round(rand(1)))*m/max(m,1)*...
            1/(j*i)^(2^(round(rand(1))));
        uold = uold+c*cos((j-1)*pi*Y/Ly).*cos((i-1)*pi*X/Lx);
    end
end

uold = reshape(uold,N2,1);
vold = zeros(size(uold));

%% contour plots 
% maxu = max(abs(uold));
% uold = u0/100*uold/maxu; %vold = v0/100*(-uold)/maxu/d;
% figure(1);
% axes('fontsize',12)
% contourf(X,Y,reshape(uold,Ny+1,Nx+1))
% xlabel('x')
% ylabel('y')
% colorbar
% title('(u-u^*)(0)')

% gamma = 197; 
u00 = ones(N2,1)*u0; v00 = ones(N2,1)*v0;
u0u0 = [u00;u00]; v0v0 = [v00;-v00]; u0v0 = [u00;v00];
diag10 = diag([ones(N2,1);zeros(N2,1)]);
ab = zeros(size(u0u0)); ab(1:N2) = a; ab(N2+1:2*N2) = b;

%% ode solver
% Tolerance for built-in MATLAB ODE solver 
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
sol = [uold;vold];
[T,sol] = ode15s(@func2d,tspan,sol,options); 
sol = sol';
uend = sol(1:N2,end) ;
vend = sol(N2+1:2*N2,end);
toc 

%% plotting 
uplot = reshape(uend,[Ny+1,Nx+1]);
vplot = reshape(vend,[Ny+1,Nx+1]);
f1 = figure(1); 
subplot(1,2,1)
surf(X,Y,uplot); xlabel('x'); ylabel('y'); title('Activator u'); colorbar
axis equal 
set(gca,'fontsize',18)
subplot(1,2,2)
surf(X,Y,vplot); xlabel('x'); ylabel('y'); title('Inhibitor v'); colorbar
axis equal 
set(gca,'fontsize',18)
f1.Position(3:4) = [1000 400]

f2 = figure (2)
subplot(1,2,1)
contourf(X,Y,uplot);xlabel('x'); ylabel('y'); title('Activator u'); colorbar
axis equal 
set(gca,'fontsize',18)
subplot(1,2,2)
contourf(X,Y,vplot);xlabel('x'); ylabel('y'); title('Inhibitor v'); colorbar
axis equal 
set(gca,'fontsize',18)
f2.Position(3:4) = [1000 400]

function [D,x] = cheb(N)
if N==0, D=0;
    x=1;
    return,
end
x = cos(pi*(0:N)/N)';
c = [2; ones(N-1,1); 2].*(-1).^(0:N)'; X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1)));
D = D - diag(sum(D'));
% off-diagonal entries
end

function w=func2d(t,u)
global gamma diag10 u0u0 v0v0 
global N2 Nx Ny
%global LL
global d Dx Dy
global u0v0 ab
%uu
uu = [u(1:N2);u(1:N2)];
% v -v
vv = [u(N2+1:2*N2);-u(N2+1:2*N2)];

ux = Dx*u;
ux = reshape(ux,Ny+1,2*Nx+2); 
ux(:,1) = 0;
ux(:,Nx+1) = 0;
ux(:,Nx+2) = 0;
ux(:,2*Nx+2) = 0;
ux = reshape(ux,2*N2,1);

uy = Dy*u;
uy = reshape(uy,Ny+1,2*Nx+2); 
uy(1,:) = 0;
uy(Ny+1,:) = 0;
uy = reshape(uy,2*N2,1);

lapl = Dx*ux+Dy*uy;
lapl(N2+1:2*N2) = d*lapl(N2+1:2*N2);

w = gamma*(ab-diag10*(u+u0v0)+(uu+u0u0).^2.*(vv+v0v0))+lapl;
end