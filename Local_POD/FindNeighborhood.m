function Neighborhood = FindNeighborhood(y_star,embedding,number_of_neighbors);
%calculate distance to other points
d_star = L2_distance(y_star,embedding,0);

[d_star_sorted, index] = sort(d_star);
Neighborhood = index(1:number_of_neighbors); %index of the neighbors