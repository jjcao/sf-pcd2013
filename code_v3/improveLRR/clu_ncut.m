function [idx] = clu_ncut(L,K)
% this routine groups the data X into K subspaces by NCut
% inputs:
%       L -- an N*N affinity matrix, N is the number of data points
%       K -- the number of subpaces (i.e., clusters)
L = abs(L)+abs(L');


A=find(sum(L,2)~=0);
B=zeros(size(L));
B(A,A) = eye(length(A)) -L(A,A)./repmat(sum(L(A,:),2),1,length(A));
L=B;
%%%%%%%%
[U,S,V] = svd(L);
V = U(:,end-K+1:end);

idx = kmeans(V,K,'emptyaction','singleton','replicates',10,'display','off');

   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
