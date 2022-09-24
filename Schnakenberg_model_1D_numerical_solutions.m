% This program numerically solves the 1D Schnakenberg reaction-diffusion
% system using the finite element method formulated by ?A Numerical 
% Approach to the Study of Spatial Pattern Formation? by Madzvamuse, 2000/ 
% Mathematical Biology modeling project, 2018 Michaelmas Term, Oxford
% University. 

clear all; close all; 
%PARAMETERS
gamma = 197; a = 0.1; b = 0.9; Du = 1; Dv = 40;

% SET UP GRIDPOINTS 
% Spacial grid 
a1 = 0; b1 = 2;
% J = (b-a)/0.5;
J = 200;
h = (b1-a1)/(J+1);
x = @(i) a1+ i*h;
X = x([0:J+1]');
% Time grid 
T1 = 1;
% N = max(32*J,T1/0.1);
% N = T1/0.1;
N = 32*J;
dt = T1/N;
t = @(n) dt*n;
T = t([0:N]');

% IMPLEMENT BOUNDARY CONDITIONS
U = zeros(J+2,N+1);
V = zeros(J+2,N+1);

% COMPUTE U0 
eps = 0.01*randn(J+2,1);
U(:,1) = 1+eps; V(:,1) = 0.9+eps;

% CONSTRUCT STIFFNESS AND MASS MATRICES
S = ( 1/h * full(gallery('tridiag',J+2,-1,2,-1))); 
S(1,1) = 1/h; S(end,end) = 1/h;
M = ( h * full(gallery('tridiag',J+2,1/6,2/3,1/6)));
M(1,1) = h/3; M(end,end) = h/3;
 
% COMPUTE LOCAL FORCE VECTOR
% localF = zeros(2*(J+1),1);
% for l = 1:J+1
%     localF(2*(l-1)+1:2*(l-1)+2) = h/2*[1;1]; 
% end

% CONSTRUCT GLOBAL FORCE VECTOR
F = ones(J+2,1); F(2:J+1) = 2; F = h/2*F;
% F = zeros(J+2,1);
% for l = 1:J+1
%     for i = 1:2
%             F(l-1+i) = F(l-1+i) + localF((l-1)*2+i);
%     end
% end

% SLOVE FOR U & V
for n = 2:N+1
    % CONSTRUCT NONLINEAR LOCAL MATRICES K1 & K2
    for l = 1:J+1
        %         N1 = @(x) (x-X(l))./h; N2 = @(x) (X(l+1) - x)./h;
        k11 = h/5*U(l,n-1)*V(l,n-1) + h/20*(U(l,n-1)*V(l+1,n-1) + U(l+1,n-1)*V(l,n-1)) + h/30*U(l+1,n-1)*V(l+1,n-1);
        k12 = h/20*U(l,n-1)*V(l,n-1) + h/30*(U(l,n-1)*V(l+1,n-1) + U(l+1,n-1)*V(l,n-1)) + h/20*U(l+1,n-1)*V(l+1,n-1);
        k21 = k12;
        k22 = h/30*U(l,n-1)*V(l,n-1) + h/20*(U(l,n-1)*V(l+1,n-1) + U(l+1,n-1)*V(l,n-1)) + h/5*U(l+1,n-1)*V(l+1,n-1);
        q11 = h/5*U(l,n-1)^2  +h/30*U(l+1,n-1)^2 + 2*h/20*U(l,n-1)*U(l+1,n-1);
        q12 = h/20*U(l,n-1)^2  +h/20*U(l+1,n-1)^2 + 2*h/30*U(l,n-1)*U(l+1,n-1);
        q21 = q12;
        q22 = h/30*U(l,n-1)^2  +h/5*U(l+1,n-1)^2 + 2*h/20*U(l,n-1)*U(l+1,n-1);
        localK1(2*(l-1)+1:2*(l-1)+2,:) = [k11, k12; k21,k22];
        localK2(2*(l-1)+1:2*(l-1)+2,:) = [q11,q12;q21,q22];
    end
    
    % CONSTRUCT GLOBAL NONLINEAR MATRICES K1 & K2
    K1 = zeros(J+2,J+2); K2 = zeros(J+2,J+2);
    for l = 1:J+1
        for i = 1:2
            for j = 1:2
                K1(l-1+i,l-1+j) = K1(l-1+i,l-1+j) + localK1((l-1)*2+i,j);
                K2(l-1+i,l-1+j) = K2(l-1+i,l-1+j) + localK2((l-1)*2+i,j);
            end
        end
    end
    
    A1 = sparse(M + dt*gamma*M + dt*Du*S);
    RHS1 = dt*gamma*a*F + M*U(:,n-1) + dt*gamma*K1*U(:,n-1);
    U(:,n) = A1\RHS1;
    A2 = sparse(M + dt*Dv*S);
    RHS2 = dt*gamma*b*F + M*V(:,n-1) - dt*gamma*K2*V(:,n-1);
    V(:,n) = A2\RHS2;
    
end

% PLOTTING
figure(1)
subplot(2,1,1); plot(X,U(:,end),'linewidth',2); title('activator u'); 
xlabel('x'); ylabel('u'); 
% xlim([0,2]); ylim([0,2]);
set(gca,'fontsize',13)
subplot(2,1,2); plot(X,V(:,end),'linewidth',2); title('inhibitor v')
xlabel('x'); ylabel('v'); 
set(gca,'fontsize',13)
% xlim([0,2]); ylim([0,2]);
figure(3)
subplot(2,1,1); mesh(X,T,U'); xlabel('x - distance'); ylabel('t');title('activator u')
zlabel('u'); colorbar
view([0,90])
set(gca,'fontsize',13)
subplot(2,1,2); mesh(X,T,V'); xlabel('x - distance'); ylabel('t');title('inhibitor v')
zlabel('v'); colorbar
view([0,90])
set(gca,'fontsize',13)
