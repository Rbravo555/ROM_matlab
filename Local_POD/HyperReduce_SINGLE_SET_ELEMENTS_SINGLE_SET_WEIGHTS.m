function [elements, weights] = HyperReduce_SINGLE_SET_ELEMENTS_SINGLE_SET_WEIGHTS(ResidualProjected)

for j=1:size(ResidualProjected,2)
    [U_svd,S_svd,~] = svd(ResidualProjected{j}, 'econ');

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
   
%     % single set of selected elements and positive weights
%     W = ones(size(U_svd,1),1);
%     DATA = [] ;
%     DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm           
%     DATA.TOLFilterCandidatePoints = 1e-10;
%     if j>1
%         DATA.IND_POINTS_CANDIDATES = elements;
%     end
%     [elements,weights] = EmpiricalCubatureMethod_CANDcompl(U_svd,S_svd,W,DATA);   
    
    if j==1
        U_combined = U_svd;
        S_combined = S_svd;
        
    else
        U_combined = [U_combined,U_svd];
        S_combined = [S_combined;S_svd];
    end
    
    %last cluster contains the set of elements that integrate the last
    %basis(including as much as possible the elements that integrate other
    %bases...) If true, this approach is not good enough!!
    
end


%taking a second svd

[U_svd,S_svd,~] = svd(U_combined, 'econ');

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



% single set of selected elements and positive weights
W = ones(size(U_svd,1),1);
DATA = [] ;
DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm           
DATA.TOLFilterCandidatePoints = 1e-10;
[elements,weights] = EmpiricalCubatureMethod(U_svd,S_svd,W,DATA);   




%  %  % naive way, many sets of elements and corresponding weights

% W_svd = ones(size(ResidualProjected,1),1);
% % Empirical Cubature Method
% DATA = [] ; 
% DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm
% DATA.TOLFilterCandidatePoints = 1e-10;
% [elements,weights]= EmpiricalCubatureMethod(U_svd,S_svd,W_svd,DATA);



