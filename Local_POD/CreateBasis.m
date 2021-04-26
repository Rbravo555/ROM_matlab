function [U] = CreateBasis( snapshots, svd_tolerance )

if nargin<2
    svd_tolerance = 1e-4;
end

%taking the SVD
[U,Sigma,Vt] = svd(snapshots);
Sigma = diag(Sigma);

%selecting modes
FrobeniusTolerance = svd_tolerance;
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

modes_taken = length(Sigma);

disp(['Number of modes considered : ', num2str(modes_taken)]);

%pruning U
modes_taken = length(Sigma);
U=U(:,1:modes_taken);

% % see the modes
% figure
% for i=1:modes_taken
%     plot(linspace(1,size(U,1),size(U,1)), U(:,i)*Sigma(i), 'LineWidth',1.5)
%     hold on
% end
% hold off
% title('Modes Visualization (\Phi_i * \sigma_i)')
end   