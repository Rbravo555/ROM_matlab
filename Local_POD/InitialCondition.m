function [u, du] = InitialCondition(X)
%     x0 = 0.2; sigma = 0.12; 
%     ind = find( abs(X-x0)<=sigma ); 
%     u = zeros(size(X)); 
%     u(ind) = 0.5*(1 + cos( (X(ind)-x0)*pi/sigma ));
a = 0.5;


sigma = 0.05;
mu = 0.2;
u = (1 / sqrt(2 * pi * sigma^2)) * exp(-(X - mu).^2/(2*sigma^2));


du = u .* (a/(sigma^2)).*(X - mu);

end