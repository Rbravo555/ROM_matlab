function reconstructed = IsomapBackMapping(y_star, vector_y, vector_x)

G = GrammMatrix(y_star,vector_y);
C = ConstrainedMatrix(y_star,vector_y);
Gtilde = G+C;
results = LagrangeMultiplierSolve(Gtilde);

%reconstructing
reconstructed = vector_x * results;
end