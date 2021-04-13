function [ResidualsProjected] = GetResidual(X,T,referenceElement, system_info) 


% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeights; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 

% Number of nodes and elements
nPt = length(X); 
nElem = size(T,1); 
a = system_info.a;
dt = system_info.dt;
n = system_info.n;
idx_state = system_info.idx_state;
Phi = system_info.Bases{idx_state};

ResidualsProjected = zeros(nElem, size(Phi,2)*2);

% New system is of compact dimensions
%M = zeros(size(Phi,2),size(Phi,2) );
%K = zeros(size(Phi,2),size(Phi,2) ); 
%C = zeros(size(Phi,2),size(Phi,2) );

% Loop on elements
for ielem=1:nElem
    Te = T(ielem,:); 
    Xe = X(Te);    
    h = Xe(end) - Xe(1);
    Me = zeros(nen); 
    Ke = zeros(nen); 
    Ce = zeros(nen); 
    % Loop on Gauss points
    for ig = 1:ngaus
        N_ig = N(ig,:);
        Nx_ig = Nxi(ig,:)*2/h;
        w_ig = wgp(ig)*h/2;        
        
        Me = Me + w_ig*(N_ig'*N_ig);
        Ke = Ke + w_ig*(Nx_ig'*Nx_ig); 
        Ce = Ce + w_ig*(N_ig'*Nx_ig);
    end
    
    A = Me + 1/2*a*dt*Ce;
    B = Me - 1/2*a*dt*Ce;
    
    
    Phie = Phi(Te,:); %The elemental basis
    ResidualsProjected(ielem,:) = [Phie'* ( A * system_info.u_pod(Te, n) - B * system_info.u_pod(Te, n-1));  Phie'* ( A * system_info.u_pod(Te, n-1) - B * system_info.u_pod(Te, n-1)) ];    
    
    %[Phie'* ( A * system_info.u_pod(Te, n) - B * system_info.u_pod(Te, n-1) ) * weights{kkk}(mmm) ; Phie'* ( A * system_info.u_pod(Te, n-1) - B * system_info.u_pod(Te, n-1) ) * weights{kkk}(mmm)]
    
end



