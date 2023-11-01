function [NORM2xest,NORM2x] = Performance_Test(K,n,x, xest)

NORM2xest=zeros(n,K);
NORM2x=zeros(n,K);

for j=1:K
   for i=1:n 
    NORM2xest(i,j)=(x(i,j) - xest(i,j))'*(x(i,j) - xest(i,j));
    NORM2x(i,j)=x(i,j)'*x(i,j);
   end
end

    
