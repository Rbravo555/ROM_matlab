% This program solves a one-dimensional transient convection equation
%        u_t + a u_x = 0
% with Dirichlet boundary conditions using the finite element method.

clear; close all; clc
addpath('./Isomap')


disp(' ')
disp('This program solves a transient convection equation')
disp('         u_t + a u_x = 0')

% Reference element: numerical quadrature and shape functions
p = 1; 
referenceElement = SetRefereceElement(p); 

% PDE coefficients
disp(' ')
a  = 0.5; %Convection coefficient a

% Spatial discretization
dom = [0,1]; 
example.dom = dom;  % Computational domain
disp(' ')
nElem = 200; %Number of elements
nPt = nElem + 1;
if p ~= 1
    warning('This program is only working for linear elements'); 
end
h = (dom(2) - dom(1))/nElem;
X = (dom(1):h:dom(2))'; 
T = [1:nPt-1; 2:nPt]';   

% Initial condition
problem = 1; 
[u0, du0] = InitialCondition(X); 

% Time discretization
tEnd = 1.2;%cinput('End time', 0.6);
nStep = 600;%cinput('Number of time-steps', 40);
dt = tEnd / nStep; 
Courant = a*dt/h;
disp(['Courant number: ', num2str(Courant)]); 

[M,K,C] = FEM_matrices(X,T,referenceElement); 

%Get system contributions (Crank-Nicolson + Galerkin)
A = M + 1/2*a*dt*C;
B = -a*dt*C; 

% Reduced system to impose du(0) = 0: 
ind_unk = 2:nPt; 
A = A(ind_unk,ind_unk); 
if problem == 4
    f = [B(2,1); zeros(nPt-2,1)];
else
    f = zeros(nPt-1,1); 
end
B = B(ind_unk,ind_unk); 

% Solution at each time step
u = zeros(nPt,nStep+1); 
u(:,1) = u0;

for n = 1:nStep
    Du = A\(B*u(ind_unk,n) + f); 
    u(ind_unk,n+1) = u(ind_unk,n) + Du; 
end


%% Solving using Isomap 
number_of_neighbors = 7;
dimensions_of_embedding = 1;
embedding = RunIsomap( u, dimensions_of_embedding, number_of_neighbors);

%Linearliy mapping from Parameters to Embedding
mapping_try1 = linspace(embedding(1), embedding(end),601);

u_isomap = zeros(nPt,nStep+1); 
u_isomap(:,1) = u0;

%Reconstructing
for n = 2:nStep+1
  y_star = mapping_try1(n);
  Neighborhood = FindNeighborhood(y_star,embedding,number_of_neighbors);
  reconstructed = IsomapBackMapping(y_star, embedding(:,Neighborhood), u(:,Neighborhood));
  u_isomap(:,n) = reconstructed;  
end

%% plotting isomap
figure(1); clf;
plot(X,u(:,1),'k--', X,u_isomap(:,nStep+1),'r','LineWidth',1.5);
l = legend('Initial condition', ['Isomap at t=',num2str(tEnd)]); 
set(l,'FontSize',24,'Location','NorthWest')
%axis([dom(1), dom(2),-0.25,1.25])
axis([dom(1), dom(2),-0.25,max(u(:,1))+.2])
set(gcf,'color',[1,1,1])


figure(2); clf;
xx = linspace(X(1),X(end),2*nElem+1); 
[uu, duu] = InitialCondition(xx - a*tEnd);
plot(xx,uu,'k--', X,u_isomap(:,nStep+1),'r','LineWidth',1.5);
l = legend('Exact solution', [' Isomap']); 
set(l,'FontSize',24,'Location','NorthWest')
%axis([dom(1), dom(2),-0.25,1.25])
axis([dom(1), dom(2),-0.25,max(u(:,1))+.2])
v = axis; 
text(v(1)+0.15,v(3)+0.15,['Cn=',num2str(Courant)],'FontSize',32)
set(gcf,'color',[1,1,1])


figure(3); clf;
fprintf('\n Press any key to see simulation results \n')
pause;
for n = [1:round(nStep/30):nStep+1,nStep+1]
    tt = (n-1)*dt;
    plot(X,u(:,n),'k:', X,u_isomap(:,n),'r', 'LineWidth',2)
    %axis([dom(1), dom(2),-0.25,1.25])
    axis([dom(1), dom(2),-0.25, max(u(:,1))+.2 ])
    set(gcf,'color',[1,1,1])
    l = legend('FOM', ['Isomap']); 
    set(l,'FontSize',24,'Location','NorthWest')
    title(['t=',num2str(tt)])
    pause(0.01)
end
