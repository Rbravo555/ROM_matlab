function G = GrammMatrix(y_star, vector_y);

x  = -vector_y + y_star;
G = x'*x;


end