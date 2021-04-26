function [elements, weights] = HyperReduce(ResidualProjected)

[U_svd,S_svd,V_svd] = svd(ResidualProjected);

Sigma = diag(S_svd);

%selecting modes
FrobeniusTolerance = 1e-4;%cinput('Tolerance for modes truncation', 0.01);
DOWN =sum(Sigma.^2);
UP=DOWN;

for i=1:length(Sigma)
    UP = UP - Sigma(i)^2;
    if sqrt(UP/DOWN)<FrobeniusTolerance
    %pruning Sigma
        Sigma = Sigma(1:i);
        break
    end
end

%pruning U
modes_taken = length(Sigma);
U_svd=U_svd(:,1:modes_taken);
S_svd=Sigma;%S_svd(1:modes_taken);

W_svd = ones(size(ResidualProjected,1),1);
% Empirical Cubature Method
DATA = [] ; 
DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm
DATA.TOLFilterCandidatePoints = 1e-10;
[elements,weights]= EmpiricalCubatureMethod(U_svd,S_svd,W_svd,DATA);


%  %plotting minimization results
%  figure
%  plot(U_svd' * ones( size(U_svd,1), 1) - (U_svd(elements,:))' * weights)
%  pause
 

