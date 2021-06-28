function [Bases, Centers] = CreateLocalBases( snapshots ,number_of_clusters, plotting,add_overlapping, r )

if nargin<4
  add_overlapping = true; %threshold for overlaping 10%
end


if nargin<5
  r = 0.1; %threshold for overlaping 10%
end

% k-means
opts = statset('Display','final');

%explore 10 different initalization of the centroids, pick the best
[idx,Centers] = kmeans(snapshots',number_of_clusters,'Replicates',10,'Options',opts);  

if plotting
% plot distribution of snapshots in the clusters
histogram(idx)
title('Number of snapshots per cluster')
end
% Split snapshots into subsets
for i=1:number_of_clusters
    MyData.States{i} = snapshots(:,idx==i);
    MyData.Neighbors{i} = [];    
end

id_neighbors = zeros(size(idx));

if add_overlapping
    if number_of_clusters~=1
        % Add overlapping
        for i=1:size(snapshots,2)
            %identify the two nearest clusters centroids to state i and mark these clusters as neighbors
            [~,idx_state] = pdist2(Centers,snapshots(:,i)','euclidean','Smallest',2);
            id_neighbors(i) = idx_state(2);
            if ~sum(ismember(MyData.Neighbors{idx(i)},id_neighbors(i)))
              MyData.Neighbors{idx(i)} = [MyData.Neighbors{idx(i)},id_neighbors(i)];
            end
        end


        for j=1:size(MyData.States,2)
            N_snaps = size(MyData.States{j},2);
            N_neighbors = size(MyData.Neighbors{j},2);

            N_add = ceil(N_snaps*r / N_neighbors); %number of snapshots to add to subset j

            for i=1:N_neighbors
                [~,idx_state] = pdist2(MyData.States{MyData.Neighbors{j}(i)}',Centers(j,:),'euclidean','Smallest',N_add);
                if i==1
                    snapshots_to_add{j}= MyData.States{MyData.Neighbors{j}(i)}(:,idx_state);
                else
                    snapshots_to_add{j} = [snapshots_to_add{j}, MyData.States{MyData.Neighbors{j}(i)}(:,idx_state)];
                end
                %MyData.States{j} = [MyData.States{j},MyData.States{MyData.Neighbors{j}(i)}(:,idx_state)];
            end
        end
        
        for j=1:size(MyData.States,2)
            MyData.States{j} = [MyData.States{j}, snapshots_to_add{j}];
        end
        
    end
end





fprintf('\n\n\n')

for j=1:length(MyData.States)
    %taking the SVD
    [U,Sigma,Vt] = svd(MyData.States{j});
    Sigma = diag(Sigma);

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

    modes_taken = length(Sigma);
    
    disp(['Number of modes considered : ', num2str(modes_taken)]);

    %pruning U
    modes_taken = length(Sigma);
    U=U(:,1:modes_taken);
    
    if plotting
    % see the modes
    figure
    for i=1:modes_taken
        plot(linspace(1,size(U,1),size(U,1)), U(:,i)*Sigma(i), 'LineWidth',1.5)
        hold on
    end
    hold off
    title('Modes Visualization (\Phi_i * \sigma_i)')
    end
    
% %     for m=1:size(U,1)
% %         for k=1:size(U,2)
% %             if abs(U(m,k))<1e-3
% %                 U(m,k)=0;
% %             end
% %         end
% %     end
    
    Bases{j} = U;
    
end