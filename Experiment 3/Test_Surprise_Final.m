function [Surprise, Mov_Surprise , AVG_Surprise, AVG_inv_PS, AVG_S,SIGMA]=Test_Surprise_Final(K, m, nu, PS, R)


%% Calculate Bayesian Surprise Bounds (Approximate Surprise, Minimum and Maximum value of Surprise)
Surprise=zeros(m,K); % Bayesian Surprise 
for j=1:K    
    Surprise(:,j)= 0.5 * nu(:,j)'*pinv(PS(:,j))*nu(:,j);
end

%% Calculates Sample-based Bayesian Surprise
AVG_Surprise=zeros(m,K); % Sample-based Bayesian Surprise
sumS=zeros(m,K);
for i=2:K
   sumS(:,i)= Surprise(:,i) + sumS (:,i-1);
   AVG_Surprise(i)= (1/i)* sumS(:,i);  
end


%% Calculates Sample-based Bayesian Surprise Covariance and Innovation Covariance

AVG_S=zeros(m,K); % Sample-based Innovation Covariance
sumI=zeros(m,K);
for i=2:K
   sumI(:,i)= nu(:,i)*nu(:,i)' + sumI (:,i-1);
   AVG_S(i)= (1/i)* sumI(:,i);  
end

AVG_inv_PS=zeros(m,K); % Sample-based Bayesian Surprise Covariance
for i=2:K
   AVG_inv_PS(:,i)=inv(R)-pinv(AVG_S(i));   % R-R*inv((R-AVG_S(:,i)))*R;
end

%% Compute the eigenvector of R,S,PS 

SIGMA=ones(m,K); 
invPS=zeros(m,K);
invR=inv(R);

 for j=2:K
   invPS(:,j)=inv(AVG_inv_PS(:,j));
   [V_PS, D_PS]=eig(invPS(:,j));
    for i=1:m
        SIGMA(i,j)=V_PS(i)'*((invR)*AVG_S(:,j) - eye(m))* V_PS(i) ;
    end
 end

%% Calculating Moving average of Bayesian Surprise

Mp=200; % number of points to average
Mov_Surprise=movmean(Surprise,Mp); % Moving average of Bayesian Surprise







