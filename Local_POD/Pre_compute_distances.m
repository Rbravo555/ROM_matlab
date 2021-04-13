function [d,w] = Pre_compute_distances(u0, centers, bases)

number_of_bases = size(bases,2);
d = -1.*ones(number_of_bases,number_of_bases);
w = cell(number_of_bases,number_of_bases, number_of_bases);

% %regularized centroids
centers = (centers'-u0)';

%d - loop
for m=1:number_of_bases
    k=m+1;
    for p=k:number_of_bases
        %fprintf('%d,%d\n',i,j)
        %d(m,p) = m*p;
        d(m,p) = ((centers(m,:))*(centers(m,:)')) - ((centers(p,:))*(centers(p,:)'));
        %d(p,m) = -m*p;
        d(p,m) = - d(m,p);
        for l=1:number_of_bases
            w{l,m,p} = 2*(bases{l})'*((centers(p,:)') - (centers(m,:)') );
        end
    end
end
k=66;
%d_m_p : distance from current time step to center m and center p. If value is negative, current state is closer to center m
