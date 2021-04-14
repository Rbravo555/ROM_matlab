function [elements, weights, global_elements] = HyperReduce_SINGLE_SET_ELEMENTS_MANY_SETS_WEIGHTS(ResidualProjected)

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
   
    % single set of selected elements and positive weights
    W = ones(size(U_svd,1),1);
    DATA = [] ;
    DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm           
    DATA.TOLFilterCandidatePoints = 1e-10;
    if j>1
        DATA.IND_POINTS_CANDIDATES = global_elements;
    end
    [e,w] = EmpiricalCubatureMethod_CANDcompl(U_svd,S_svd,W,DATA);%EmpiricalCubatureMethod_CANDcompl or  EmpiricalCubatureMethod
    elements{j} = e;
    weights{j} = w;
    if j==1
       global_elements = e;
    else
       global_elements = unique([global_elements;e]);
    end
    
    
    
    
end

% W_svd = ones(size(ResidualProjected,1),1);
% % Empirical Cubature Method
% DATA = [] ; 
% DATA.IncludeSingularValuesF  = 0 ; % Singular Values are not included in the minimization norm
% DATA.TOLFilterCandidatePoints = 1e-10;
% [elements,weights]= EmpiricalCubatureMethod(U_svd,S_svd,W_svd,DATA);



