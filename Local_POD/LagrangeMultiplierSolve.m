function results = LagrangeMultiplierSolve(Gtilde)
%build rhs
n = size(Gtilde,1);
A = ones(n+1,n+1);
A(1:n,1:n) = 2*Gtilde;
A(end,end) = 0;

b = zeros(n+1,1);
b(end) = 1;

results = linsolve(A,b);
results = results(1:end-1);
end