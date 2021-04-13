function [M,K,C] = FEM_matrices(X,T,referenceElement) 


% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeights; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 

% Number of nodes and elements
nPt = length(X); 
nElem = size(T,1); 

M = zeros(nPt);
K = zeros(nPt); 
C = zeros(nPt);
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
    % Assembly
    M(Te,Te) = M(Te,Te) + Me; 
    K(Te,Te) = K(Te,Te) + Ke; 
    C(Te,Te) = C(Te,Te) + Ce;
end




