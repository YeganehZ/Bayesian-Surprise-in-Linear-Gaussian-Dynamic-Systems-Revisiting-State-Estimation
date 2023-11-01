clc
clear all


K=200; % Time index
MC=2000;  % Number of Monte Carlo runs

%% Parameters of State-space Model  

F = [0.82 0; 0  1];  % Transition Matrix
H = [1  1];          % Measurement Matrix
Q = [1  0; 0  1];    % State Noise Covariance
R = 1;               % Measurement Noise Covariance 

m=size(H,1);         % Dimension of the Measurement, innovation 
n=size(H,2);         % Dimension of state vector, state estimate, .... 


%% Defining vectors 

%State Vectors
x=zeros(n,K);    % state vector
xpred=zeros(n,K); % predicted state vector k|k-1
xest=zeros(n,K); % estimate state vector k|k
xpred_S=zeros(n,K); % predicted state vector k|k-1 only surprise
xest_S=zeros(n,K); % estimate state vector k|k  only surprise
xpred_IN=zeros(n,K); % predicted state vector k|k-1 only innovation
xest_IN=zeros(n,K); % estimate state vector k|k  only innovation

% Error Covariance
Ppred=zeros(n,n,K); % predicted error covariance k|k-1
Pest=zeros(n,n,K);  % estimate error covariance at k|k
PSpred=zeros(m,K);  % predicted error surprise-based covariance at k|k
PSest=zeros(m,K);   % estimate error surprise-based covariance k|k-1 
nu=zeros(m,K);     %innovation vector
Spred=zeros(m,K);  % predicted innovation covariance
Sest=zeros(m,K);   % estimated innovation covariance
KG=zeros(n,K);     % Kalman Gain

% Error Covariance for only surprise
Ppred_S=zeros(n,n,K); % predicted error covariance k|k-1
Pest_S=zeros(n,n,K);  % estimate error covariance at k|k
PSpred_S=zeros(m,K);  % predicted error surprise-based covariance at k|k
PSest_S=zeros(m,K);   % estimate error surprise-based covariance k|k-1 
nu_S=zeros(m,K);     %innovation vector
Spred_S=zeros(m,K);  % predicted innovation covariance
Sest_S=zeros(m,K);   % estimated innovation covariance
KG_S=zeros(n,K);     % Kalman Gain


% Error Covariance for only innovation
Ppred_IN=zeros(n,n,K); % predicted error covariance k|k-1
Pest_IN=zeros(n,n,K);  % estimate error covariance at k|k
PSpred_IN=zeros(m,K);  % predicted error surprise-based covariance at k|k
PSest_IN=zeros(m,K);   % estimate error surprise-based covariance k|k-1 
nu_IN=zeros(m,K);     %innovation vector
Spred_IN=zeros(m,K);  % predicted innovation covariance
Sest_IN=zeros(m,K);   % estimated innovation covariance
KG_IN=zeros(n,K);     % Kalman Gain

z=zeros(m,K);       % measurement vector

% Performance test vectors
Cnt=zeros(K,MC);    % C=R*inv(Spred)+ R*inv(PSpred) at each monte-carlo 
Cnt_S=zeros(K,MC);    % C=R*inv(PSpred) at each monte-carlo 
Cnt_IN=zeros(K,MC);    % C=R*inv(Spred) at each monte-carlo 

NORM2est_MC=zeros(n,K,MC);  % defines the norm 2 vector between xest and x for all monte-carlo simulation 
NORM2_MC=zeros(n,K,MC);     % defines the norm 2 vector of x for all monte-carlo simulation 

NORM2est_S_MC=zeros(n,K,MC);  % defines the norm 2 vector between xest and x for all monte-carlo simulation 
NORM2_S_MC=zeros(n,K,MC);     % defines the norm 2 vector of x for all monte-carlo simulation 

NORM2est_IN_MC=zeros(n,K,MC);  % defines the norm 2 vector between xest and x for all monte-carlo simulation 
NORM2_IN_MC=zeros(n,K,MC);     % defines the norm 2 vector of x for all monte-carlo simulation 

RMSRE=zeros(n,K);    % Root Mean Square Relative Error (RMSRE) (without considering consistancy test)
RMSRE_S=zeros(n,K);    % Root Mean Square Relative Error only Surprise (RMSRE) (without considering consistancy test)
RMSRE_IN=zeros(n,K);    % Root Mean Square Relative Error only innovation (RMSRE) (without considering consistancy test)



%% Initial Conditions

x(:,1) =[1; 1];
z(:,1)=H*x(:,1);
xpred(:,1) =[1; 1];
xest(:,1) = [1; 1];
Ppred(:,:,1) =[ 10  0; 0 10 ];
Pest(:,:,1) =[ 10  0; 0 10 ];
Spred(:,1)=R + H*Pest(:,:,1)*H';
Sest(:,1)=R + H*Pest(:,:,1)*H';
PSpred(:,1) = Sest(:,1)*pinv(H*Pest(:,:,1)*H')*R;
PSest(:,1) = Sest(:,1)*pinv(H*Pest(:,:,1)*H')*R;
KG(:,1)=Pest(:,:,1)*H'*inv(R);
nu(:,1)=z(:,1)-H*xpred(:,1);

xpred_S(:,1) =[1; 1];
xest_S(:,1) = [1; 1];
Ppred_S(:,:,1) =[ 10  0; 0 10 ];
Pest_S(:,:,1) =[ 10  0; 0 10 ];
Spred_S(:,1)=R + H*Pest_S(:,:,1)*H';
Sest_S(:,1)=R + H*Pest_S(:,:,1)*H';
PSpred_S(:,1) = Sest_S(:,1)*pinv(H*Pest_S(:,:,1)*H')*R;
PSest_S(:,1) = Sest_S(:,1)*pinv(H*Pest_S(:,:,1)*H')*R;
KG_S(:,1)=Pest_S(:,:,1)*H'*inv(R);
nu_S(:,1)=z(:,1)-H*xpred_S(:,1);

xpred_IN(:,1) =[1; 1];
xest_IN(:,1) = [1; 1];
Ppred_IN(:,:,1) =[ 10  0; 0 10 ];
Pest_IN(:,:,1) =[ 10  0; 0 10 ];
Spred_IN(:,1)=R + H*Pest_IN(:,:,1)*H';
Sest_IN(:,1)=R + H*Pest_IN(:,:,1)*H';
PSpred_IN(:,1) = Sest_IN(:,1)*pinv(H*Pest_IN(:,:,1)*H')*R;
PSest_IN(:,1) = Sest_IN(:,1)*pinv(H*Pest_IN(:,:,1)*H')*R;
KG_IN(:,1)=Pest_IN(:,:,1)*H'*inv(R);
nu_IN(:,1)=z(:,1)-H*xpred_IN(:,1);


for mc=1:MC
%% System Model

Cnt(1,mc)=1;
Cnt_S(1,mc)=1;
Cnt_IN(1,mc)=1;


 for i=2:K
    x(:,i)=F*x(:,i-1)+  mvnrnd([0 0],Q)';
    z(:,i)=H*x(:,i)+ sqrt(R)*randn(1,1); % Measurement noise vector 
 end
 
%% Kalman Filter Estimation and Prediction

for j=2:K   
 [xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j)]=timeupdate_surprise(xest(:,j-1),Pest(:,:,j-1),PSest(:,j-1),Sest(:,j-1),R,H,Q,F); 
 [xest(:,j), Pest(:,:,j), PSest(:,j), Sest(:,j), KG(:,j),nu(:,j),Cnt(j,mc)] = measurementupdate_surpise(xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j),z(:,j),R,H);      
end

for j=2:K   
 [xpred_S(:,j),Ppred_S(:,:,j),PSpred_S(:,j),Spred_S(:,j)]=timeupdate_onlysurprise(xest_S(:,j-1),Pest_S(:,:,j-1),PSest_S(:,j-1),Sest_S(:,j-1),R,H,Q,F); 
 [xest_S(:,j), Pest_S(:,:,j), PSest_S(:,j), Sest_S(:,j), KG_S(:,j),nu_S(:,j),Cnt_S(j,mc)] = measurementupdate_onlysurpise(xpred_S(:,j),Ppred_S(:,:,j),PSpred_S(:,j),Spred_S(:,j),z(:,j),R,H);     
end

for j=2:K   
 [xpred_IN(:,j),Ppred_IN(:,:,j),PSpred_IN(:,j),Spred_IN(:,j)]=timeupdate_onlyinnovation(xest_IN(:,j-1),Pest_IN(:,:,j-1),PSest_IN(:,j-1),Sest_IN(:,j-1),R,H,Q,F); 
 [xest_IN(:,j), Pest_IN(:,:,j), PSest_IN(:,j), Sest_IN(:,j), KG_IN(:,j),nu_IN(:,j),Cnt_IN(j,mc)] = measurementupdate_onlyinnovation(xpred_IN(:,j),Ppred_IN(:,:,j),PSpred_IN(:,j),Spred_IN(:,j),z(:,j),R,H);     
end

%% Calculating the second norm for checking system performance in terms of RMSRE
[NORM2est,NORM2] = Performance_Test(K,n,x, xest);
[NORM2est_S,NORM2_S] = Performance_Test(K,n,x, xest_S);
[NORM2est_IN,NORM2_IN] = Performance_Test(K,n,x, xest_IN);

NORM2est_MC(:,:,mc)=NORM2est;
NORM2_MC(:,:,mc)=NORM2;  

NORM2est_S_MC(:,:,mc)=NORM2est_S;
NORM2_S_MC(:,:,mc)=NORM2_S; 

NORM2est_IN_MC(:,:,mc)=NORM2est_IN;
NORM2_IN_MC(:,:,mc)=NORM2_IN; 

mc
end

%% Calculate RMSRE for both cases (regular Kalman, surprise only Kalman, innovation only Kalman)
NOM1=sum(NORM2est_MC,3); 
DIM1=sum(NORM2_MC,3);
for j=1:K
    for i=1:n
    RMSRE(i,j)=sqrt(NOM1(i,j))/sqrt(DIM1(i,j));
    end
end

NOM2=sum(NORM2est_S_MC,3); 
DIM2=sum(NORM2_S_MC,3);
for j=1:K
    for i=1:n
    RMSRE_S(i,j)=sqrt(NOM2(i,j))/sqrt(DIM2(i,j));
    end
end

NOM3=sum(NORM2est_IN_MC,3); 
DIM3=sum(NORM2_IN_MC,3);
for j=1:K
    for i=1:n
    RMSRE_IN(i,j)=sqrt(NOM3(i,j))/sqrt(DIM3(i,j));
    end
end
% 

figure, 
p= plot(2:K, (RMSRE(2,2:K)), 'b :', 2:K, (RMSRE_S(2,2:K)), 'r --',  2:K, (RMSRE_IN(2,2:K)), 'k -');
p(1).LineWidth = 3;
p(2).LineWidth = 3;
p(3).LineWidth = 3;
leg2 = legend('\textbf{Traditional Kalman Filter}','\textbf{Kalman Filter in terms of Bayesian Surprise}', '\textbf{Kalman Filter in terms of Innovation}');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('\textbf{RMSRE}','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
axis([1 180 1 0]);
hold on





