function z = BuildZMatrix(current_u, centers);

if nargin<1
    centers = [1,2,3,4,5,6];
    current_u = 1;
end
number_of_bases = length(centers);
z = cell(number_of_bases,number_of_bases);

for i=1:length(centers)
    k=i+1;
    for j=k:length(centers)
        fprintf('%d,%d\n',i,j)
        z{i,j} = 'hi'
    end
end