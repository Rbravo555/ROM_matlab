function [M,K,C,Ml] = FEM_matricesPOD_hrom(X,T,referenceElement, Phi, w, z) 


% reference element information
nen = referenceElement.nen; 
ngaus = referenceElement.ngaus; 
wgp = referenceElement.GaussWeights; 
N = referenceElement.N; 
Nxi = referenceElement.Nxi; 

% Number of nodes and elements
nPt = length(X); 
nElem = size(T,1); 

% New system is of compact dimensions
M = zeros(size(Phi,2),size(Phi,2) );
K = zeros(size(Phi,2),size(Phi,2) ); 
C = zeros(size(Phi,2),size(Phi,2) );

% Loop on elements
for i=1:size(z)
    ielem = z(i);
    Te = T(ielem,:); 
    Xe = X(Te);
    Phie = Phi(Te,:); %The elemental basis
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
    M = M + Phie' * Me * Phie *w(i); 
    K = K + Phie' * Ke * Phie *w(i); 
    C = C + Phie' * Ce * Phie *w(i);
end




