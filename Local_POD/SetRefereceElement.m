function Element = SetRefereceElement(p)

if p == 1
	Xe_ref = [-1,1];
    zgp = [-1/sqrt(3); 1/sqrt(3)]; 
    wgp = [1   1];
    ngaus = length(wgp); 
    N =  [(1-zgp)/2     (1+zgp)/2];  
    Nxi =  [-1/2*ones(ngaus,1)  1/2*ones(ngaus,1)];
    N2xi = zeros(ngaus,2);
elseif p == 2
    Xe_ref = [-1,0,1];
    zgp = [-sqrt(15)/5; 0; sqrt(15)/5]; 
    wgp = [5/9   8/9   5/9];
    ngaus = length(wgp);
    N   = [(zgp-1).*zgp/2   1-zgp.^2     (zgp+1).*zgp/2];  
    Nxi = [zgp-1/2   -2*zgp   zgp+1/2];
    N2xi = [ones(ngaus,1)    -2*ones(ngaus,1)   ones(ngaus,1)]; 
else
    error('unavailable element')
end


Element.degree = p;
Element.Xe_ref = Xe_ref; 
Element.nen = length(Xe_ref); 
Element.ngaus = ngaus; 
Element.GaussPoints  = zgp;
Element.GaussWeights = wgp;
Element.N = N; 
Element.Nxi = Nxi; 
Element.N2xi = N2xi; 
