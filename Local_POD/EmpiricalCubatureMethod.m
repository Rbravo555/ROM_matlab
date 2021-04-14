function [z,w]= EmpiricalCubatureMethod(BasisF,SingVal_F,W,DATA)
% Empirical Cubature Method
% Joaquín A. Hernández, UPC/CIMNE. 2018 

if DATA.IncludeSingularValuesF == 1
    G = bsxfun(@times,BasisF',SingVal_F);
    %G = [G; ones(1,length(G))];  %adding a row of ones 
    b = G*sqrt(W) ;  % b Vector (exact integral)
    bEXACT = b ;
else
    G = BasisF' ;
    %G = [G; ones(1,length(G))];  %adding a row of ones
    b = G*sqrt(W) ;  % b Vector (exact integral)
    bEXACT = b.*SingVal_F ; 
end
nbEXACT = norm(bEXACT) ; 


Gnorm =sqrt(sum(G.*G,1)) ; % Norm of Columns of G
M = size(G,2) ;  % Number of FE points
DATA = DefaultField(DATA,'TOL',0) ; % Default tolerance for convergence
TOL = DATA.TOL ;
% INITIALIZATIONS
% ------------------------
z = [] ; % Set of integration points
% Set of candidate points (those whose associated column has low norm are removed)
y=1:M ;
DATA = DefaultField(DATA,'TOLFilterCandidatePoints',1e-6) ;
if DATA.TOLFilterCandidatePoints >0
    TOL_REMOVE = DATA.TOLFilterCandidatePoints*norm(b) ;
    rmvpin = find(Gnorm(y)<TOL_REMOVE) ;
    y(rmvpin) = [] ;
end
alpha = [] ; % Vector of weights
mPOS = 0 ; % Number of nonzero weights
r = b ; % Residual vector
k = 1;  % Number of iterations
errorGLO = [] ; %(for storing error)
% Default number of points
DATA = DefaultField(DATA,'npoints',length(b)) ;
m = min(DATA.npoints,length(b)) ;
% END INITIALIZATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normB = norm(b) ;
nerror = norm(r)/normB  ;
H = [] ; % Inverse of (Gz'*Gz)
ERROR_GLO = [] ;
NPOINTS =[] ; 
nerrorACTUAL = nerror; 
while  nerrorACTUAL >TOL && mPOS <m   && ~isempty(y)
    % TIMELOC_k = tic ;
    % STEP 1. Compute new point
    ObjFun = G(:,y)'*r ;
    ObjFun = ObjFun./Gnorm(y)';
    [maxLOC, indSORT] = max(ObjFun)  ;
    i = y(indSORT(1));
    % STEP 2.  Update alpha and H  (unrestricted least-squares)
    if k==1
        alpha =  G(:,i)\b ;
        H = 1/(G(:,i)'*G(:,i)) ;
    else
        [H alpha] = UpdateWeightsInverse(G(:,z),H,G(:,i),alpha,r) ;
    end
    % STEP 3. Move i from set y to set z
    z = [z;i] ;     y(indSORT(1)) = [] ;
    % STEP 4. Find possible negative weights
    n = find(alpha<=0) ;
    if  ~isempty(n)
        % STEP 5        
        y = [y z(n)'];  z(n)=[] ;
        H = MultiUpdateInverseHermitian(H,n) ;
        % Recomputing alpha
        alpha = H*(G(:,z)'*b );
        disp('Negative Weight Found');
    end
    
    % STEP 6
    r = b-G(:,z)*alpha ;
    nerror = norm(r)/norm(b) ; % Relative error (using r and b)
    if DATA.IncludeSingularValuesF == 0
        nerrorACTUAL = SingVal_F.*r ;
        nerrorACTUAL = norm(nerrorACTUAL/nbEXACT); 
    else
        nerrorACTUAL = nerror ;
    end
    % STEP 7
    disp(['k = ',num2str(k),', m=',num2str(length(z)),' ,','; error n(res)/n(b) (%) = ',...
        num2str(nerror*100),';  Actual error % =',num2str(nerrorACTUAL*100)]) ; 
    ERROR_GLO(k) = nerrorACTUAL ;
    NPOINTS(k) =  length(z) ; 

    mPOS = length(z) ; 
    k = k + 1 ;
%     if k==13
%         break;
%        %fprintf('Stop here')
%     end
    
%     if length(z) == m
%         dbstop('88')
%         disp('')
%     end
   
end



w = alpha.*sqrt(W(z)) ;

disp(['Total number of iterations =',num2str(k)])

figure(500)
hold on
xlabel('Number of points')
ylabel('Error (%)')
plot(NPOINTS,ERROR_GLO*100,'k')

% figure(501)
% hold on
% xlabel('Number of iterations')
% ylabel('Number of points')
% plot(NPOINTS,'k')

end

function  DATAGEN = DefaultField(DATAGEN,FIELDVAR,dEFval) ; 
% DefaultField adds the field FIELDVAR to structure array DATAGEN and
% assigns the specified value dEFval. 
if ~isfield(DATAGEN,FIELDVAR)
    DATAGEN.(FIELDVAR) = dEFval ; 
end

end

function [Bast x] = UpdateWeightsInverse(A,Aast,a,xold,r)

if nargin == 0
    format long g
    P = 100 ; 
    k = 49 ; 
    A = randn(P,k) ; 
    a = randn(P,1) ;
    Aast = inv(A'*A) ; 
    B = [A a] ; 
    BastReal = inv(B'*B) ;
    b = randn(P,1) ; 
    xold = Aast*(A'*b) ; 
    r = b - A*xold ; 
    xREAL = BastReal*(B'*b) ;
    
 
end
c = A'*a ; 
d = Aast*c ; 
s = a'*a-c'*d ; 
Bast = [Aast + (d*d')/s  -d/s; -d'/s  1/s] ; 
v = a'*r/s ; 
x = [(xold -d*v ); v] ;

 
end


function Ahinv = MultiUpdateInverseHermitian(Binv,jrowMAT)
% Recursive application of UpdateInverseHermitian
%  J.A. Hdez, 24 March 2017
if nargin == 0
    m = 10; n =7; 
    jrowMAT =[ 8] ; 
    A = randn(m,n) ; a = randn(m,length(jrowMAT)) ;
    Bor = zeros(m,n+length(jrowMAT)) ; 
    indOR = 1:size(Bor,2) ; 
    indOR(jrowMAT) = [] ; 
    Bor(:,indOR) = A ; 
    Bor(:,jrowMAT) = a;     
    B = [Bor'*Bor] ;  
    Binv = inv(B);     
    AhinvREAL = inv(A'*A)      ;
    
end
 
jrowMAT = sort(jrowMAT) ; 
BinvOLD = Binv ; 
for i = 1:length(jrowMAT)
    jrow = jrowMAT(i)-i+1 ; 
    Ahinv = UpdateInverseHermitian(BinvOLD,jrow) ; 
    BinvOLD = Ahinv ; 
end


end


function Ahinv = UpdateInverseHermitian(Binv,jrow)
% Let B = ([A a]^T*[A a]), where A is a m-by-n full rank  matrix, and "a" a vector of
% m entries. Suppose that we are given the inverse of B (Binv). Function
% UpdateInverseHermitian(Binv,j) returns the value of inv(A^T A) based on the block
% decomposition of % matrix Binv
%
% If the second argument jrow is passed, then it is assumed that B is given
% by
%
%  B = [A(:,1:jrow-1) a A(:,jrow:end)] ;
%
%  J.A. Hdez, 24 March 2017
if nargin == 0
    m = 10; n =4;
    jrow =1 ;
    A = randn(m,n) ; a = randn(m,1) ;
    Bor = [A(:,1:jrow-1) a A(:,jrow:end)] ;
    B = [Bor'*Bor] ;
    Binv = inv(B);
    AhinvREAL = inv(A'*A) 
    
    
end

if nargin ==1
    jrow = size(Binv,2) ;
end

if jrow == size(Binv,2)
    Ahinv = Binv(1:end-1,1:end-1) -  (Binv(1:end-1,end)*Binv(end,1:end-1))/Binv(end,end) ;
else
    AUX = [Binv(:,1:jrow-1)   Binv(:,jrow+1:end)   Binv(:,jrow)] ;
    AUX = [AUX(1:jrow-1,:) ;  AUX(jrow+1:end,:)  ; AUX(jrow,:)] ;
    
    Ahinv = AUX(1:end-1,1:end-1) -  (AUX(1:end-1,end)*AUX(end,1:end-1))/AUX(end,end) ;
    
end

end



