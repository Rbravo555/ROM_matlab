% This program solves a static one-dimensional structural mechanics
%        K u = f
% with Dirichlet boundary conditions using the finite element method.

clear all; close all; clc
svd_tolerance=cinput('Tolerance for SVD truncation', 1e-4);
 
% Reference element: numerical quadrature and shape functions
p = 1;
referenceElement = SetRefereceElement(p);

% Spatial discretization
dom = [0,1];
example.dom = dom; % Computational domain
disp(' ')
nElem = 3; %Number of elements
nStep = 30; %Number of time steps
nPt = nElem + 1;
if p ~= 1
	warning('This program is only working for linear elements');
end
h = (dom(2) - dom(1))/nElem;
X = (dom(1):h:dom(2))';
T = [1:nPt-1; 2:nPt]';

[M,K,C] = FEM_matrices(X,T,referenceElement);

% Reduced system to impose du(0) = 0:
ind_unk = 2:nPt;
K = K(ind_unk,ind_unk);
F = zeros(length(ind_unk), 1);


u = zeros(nPt, nStep);

for i=1:nStep
    F(end)=(i/nStep);
    u(ind_unk, i) = K\F;
end


%imposing u2=u3
T_red = [1,0;1,0;0,1];

for i=1:nStep
    F(end)=(i/nStep);
    K_red = T_red'* K * T_red;
    F_red = T_red' * F;
    u_hat = K_red\F_red;
    u(ind_unk, i) = T_red * u_hat;
end


Basis = CreateBasis(u, svd_tolerance);


%% Running POD simulation
[M,K,C] = FEM_matricesPOD(X,T,referenceElement,Basis);

u_pod = zeros(nPt,nStep);
F = zeros(nPt, 1);


% Solution at each time step
for i=1:nStep
    F(end)=(i/nStep);
    F_red = Basis' * F;
    q = K\F_red;
    u_pod(:, i) = Basis*q;
end


%% ErrorMeasure
%L2 error
UP = sum((u_pod - u).^2);
DOWN = sum(u.^2);
L2 = sqrt(UP/DOWN);
fprintf('\n\n\nL2 error: %e \n', L2);




%% Plot simulation results 
figure(55);clf;
for n = [1:round(nStep/30):nStep,nStep]
    tt = (n-1)*.1;
    plot(X,u(:,n),'k:', X,u_pod(:,n),'ro', 'LineWidth',2)
    axis([dom(1), dom(2),-0.25, max(u(:))+.1 ])
    set(gcf,'color',[1,1,1])
    l = legend('FOM', ['ROM']); 
    set(l,'FontSize',24,'Location','NorthWest')
    title(['t=',num2str(tt)])
    pause(0.1)
end