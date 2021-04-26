function C = ConstrainedMatrix(y_star, vector_y, epsilon);

if nargin<3
  epsilon = 0.01;
end

n = 4;
c = zeros(size(vector_y,2), 1);
x  = -vector_y + y_star;
NormColums = vecnorm(x);
MaxNorm = max(NormColums);

for i=1:length(c)
c(i)= epsilon*(norm(x(:,i)) / MaxNorm)^n;
end

C = diag(c);
end