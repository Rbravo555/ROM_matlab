% This program solves a one-dimensional transient convection equation
%        u_t + a u_x = 0
% with Dirichlet boundary conditions using the finite element method.



function L2 = Launch_simulation_SINGLE_SET_REDUCED_ELEMENTS(NumberOfClusters)

if nargin<1
  clear all; close all; clc
  NumberOfClusters=cinput('Number of clusters for the local POD', 3);
  plotting = true;
else
  plotting = false;
end


addpath('./Isomap')

disp(' ')
disp('This program solves a transient convection equation')
disp('         u_t + a u_x = 0')

% Reference element: numerical quadrature and shape functions
p = 1;
referenceElement = SetRefereceElement(p);

% PDE coefficients
disp(' ')
a = 0.5; %Convection coefficient a

% Spatial discretization
dom = [0,1];
example.dom = dom; % Computational domain
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
B = M - 1/2*a*dt*C;

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
	u(ind_unk,n+1) = Du;
end

[Bases, Centers] = CreateLocalBases( u ,NumberOfClusters,plotting);
ResidualProjected = cell(size(Bases));

[d,w] = Pre_compute_distances(u0, Centers, Bases);

w_pod = w;
d_pod = d;
z_pod = d_pod;
idx_stateOld = find(all(z_pod'<0));
idx_state = idx_stateOld;



%% Running POD simulation
u_pod = zeros(nPt,nStep+1);
u_pod(:,1) = u0;
[~,idx_stateOld] = pdist2(Centers,u0','squaredeuclidean','Smallest',1);

current_basis = Bases{idx_stateOld};

% Solution at each time step
q = zeros(size(current_basis,2),nStep+1);
q_old = current_basis' * u0;
u_pod(:,1)= current_basis*q_old;

for n = 2:nStep+1
    
    [M,K,C] = FEM_matricesPOD(X,T,referenceElement,current_basis);    
    A = M + 1/2*a*dt*C;
    B = -a*dt*C;

        
	Dq = A\(B*q_old);
    u_pod(:,n) = u_pod(:,n-1) + (current_basis * Dq);    
    %q_old = q_old + Dq;    
    
    z_pod = UpdateZMatrix(z_pod,w_pod, Dq, idx_state);
    
    idx_state = find(all(z_pod'<0));
    if idx_stateOld ~= idx_state
        %change from one cluster to another
        current_basis = Bases{idx_state} ;
        q_old = current_basis'*u_pod(:,n); %project the last step's low dim representation onto the new basis  
        idx_stateOld=idx_state;
    else
        q_old = q_old + Dq;
    end
    
    %hyper-reduction training   
    system_info.u_pod = u_pod;
    system_info.Bases = Bases;
    system_info.a = a;
    system_info.dt = dt;
    system_info.n = n;
    system_info.idx_state = idx_state;
    if isempty(ResidualProjected{idx_state})
        ResidualProjected{idx_state} = GetResidual(X,T,referenceElement, system_info);
    else
        ResidualProjected{idx_state} = [ResidualProjected{idx_state}, GetResidual(X,T,referenceElement, system_info)];
    end
end



%% Hyper-reduction
[elements, weights] = HyperReduce_SINGLE_SET_ELEMENTS(ResidualProjected);


%  %plotting minimization results
%  figure
%  for i=1:NumberOfClusters
%      plot(ResidualProjected{i}' * ones( size(ResidualProjected{i},1), 1) - (ResidualProjected{i}(elements,:))' * weights)
%      pause
%  end

 
 %visualizing the selected elements in the domain
 
 
 

w_hrom = w;
d_hrom = d;
z_hrom = d_hrom;
idx_stateOld = find(all(z_hrom'<0));
idx_state = idx_stateOld;
 
u_hrom = zeros(nPt,nStep+1);
u_hrom(:,1) = u0;

current_basis = Bases{idx_stateOld};

% Solution at each time step
q_old = current_basis' * u0;
u_hrom(:,1)= current_basis*q_old;


for n = 2:nStep+1
    
    idx_state = find(all(z_hrom'<0));
    if idx_stateOld ~= idx_state
        current_basis = Bases{idx_state};
        q_old = current_basis'*u_hrom(:,n-1); %project the last step's low dim representation onto the new basis  
        idx_stateOld=idx_state;        
    end
    
    
    [M,K,C] = FEM_matricesPOD_hrom(X,T,referenceElement,current_basis, weights, elements);
    A = M + 1/2*a*dt*C;
    B = -a*dt*C;
        
	Dq = A\(B*q_old);
    u_hrom(:,n) = u_hrom(:,n-1) + (current_basis * Dq);   
    z_hrom = UpdateZMatrix(z_hrom,w_hrom, Dq, idx_state);

    
    q_old = q_old +Dq;

end



%% ErrorMeasure
%L2 error ROM
UP = sum((u_pod - u).^2);
DOWN = sum(u.^2);
L2 = sqrt(UP/DOWN);
fprintf('\n\n\nL2 error rom: %e \n', L2);


%L2 error H-ROM
UP = sum((u_hrom - u).^2);
DOWN = sum(u.^2);
L2 = sqrt(UP/DOWN);
fprintf('\n\n\nL2 error hrom: %e \n', L2);


% if plotting
%     %% Plot simulation results
%     fprintf('\n Press any key to see simulation results \n')
%     pause;
%     figure(55);clf;
%     for n = [1:round(nStep/30):nStep+1,nStep+1]
%         tt = (n-1)*dt;
%         plot(X,u(:,n),'k:', X,u_pod(:,n),'r', 'LineWidth',2)
%         axis([dom(1), dom(2),-0.25, max(u(:,1))+.2 ])
%         set(gcf,'color',[1,1,1])
%         l = legend('FOM', ['ROM']); 
%         set(l,'FontSize',24,'Location','NorthWest')
%         title(['t=',num2str(tt)])
%         pause(0.1)
%     end
% end
% 
% 
if plotting
    %% Plot simulation results
    fprintf('\n Press any key to see simulation results \n')
    pause;
    figure(55);clf;
    for n = [1:round(nStep/30):nStep+1,nStep+1]
        tt = (n-1)*dt;
        plot(X,u(:,n),'k:', X,u_hrom(:,n),'r', 'LineWidth',2)
        axis([dom(1), dom(2),-0.25, max(u(:,1))+.2 ])
        set(gcf,'color',[1,1,1])
        l = legend('FOM', ['H-ROM']); 
        set(l,'FontSize',24,'Location','NorthWest')
        title(['t=',num2str(tt)])
        pause(0.1)
    end
end