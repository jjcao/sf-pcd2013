function [Z,E] = admreglrr(X,lambda,K,beta,rho,DEBUG)

if (~exist('DEBUG','var'))
    DEBUG = 0;
end
if nargin < 5
    rho = 1.1;
end

normfX = norm(X,'fro');
tol =  1e-3;%threshold for the error in constraint
[d n] = size(X);
maxIter = 1000;
max_mu = 1e10;
mu = 1e-1;
xtx = X'*X;
inv_x = inv(xtx+2*eye(n));

%% Initializing optimization variables
% intialize
E = zeros(d,n);
YA = zeros(d,n);
YB = zeros(n,n);
YC = zeros(n,n);
Z = zeros(n,n);
W = zeros(n,n);
L = zeros(n,n);

%% Start main loop
convergenced = 0;
iter = 0;

if DEBUG
    disp(['initial,rank(Z)=' num2str(rank(Z))]);
end

while iter<maxIter
    iter = iter + 1;
    
   temp = Z+YC/mu;
    [U,sigma,V] = svd(temp,'econ');
    sigma = diag(sigma);
    svp = length(find(sigma>1/mu));
    if svp>=1
        sigma = sigma(1:svp)-1/mu;
    else
        svp = 1;
        sigma = 0;
    end
    J = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
    
    tempe=Z+YB/mu;
    W = max(0,tempe - beta/mu*K)+min(0,tempe + beta/mu*K);
    
    Z =inv_x*(xtx-X'*E+W+J+(X'*YA-YB-YC)/mu);
    
    E = solve_l1l2(X - X*Z + YA/mu,lambda/mu);

    dYA = X - X*Z - E;
    dYB = Z - W;
    dYC = Z - J;
    stopC = max([max(max(abs(dYA))),max(max(abs(dYB))),max(max(abs(dYC)))]);
    convergenced = stopC <tol;
    
    if DEBUG
        if iter==1 || mod(iter,50)==0 || convergenced
            disp(['iter ' num2str(iter) ',mu=' num2str(mu) ...
                ',rank(Z)=' num2str(rank(Z)) ',stopALM=' num2str(stopC,'%2.3e')]);
        end
    end
    if convergenced
        break;
    else
        YA = YA + mu*dYA;
        YB = YB + mu*dYB; 
        YC = YC + mu*dYC;
        mu = min(max_mu, mu*rho);
    end
end

function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end

function [x] = solve_l2(w,lambda)
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);
end
