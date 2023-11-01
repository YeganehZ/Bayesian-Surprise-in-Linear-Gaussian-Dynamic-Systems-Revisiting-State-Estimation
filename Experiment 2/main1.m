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

%% Pre-defining state vectors 

%State and Measurement Vectors
x=zeros(n,K);     % state vector
xpred=zeros(n,K); % predicted state vector x(k|k-1)
xest=zeros(n,K);  % estimate state vector x(k|k)
z=zeros(m,K);     % measurement vector

% Error Covariance
Ppred=zeros(n,n,K); % predicted error covariance P(k|k-1)
Pest=zeros(n,n,K);  % estimate error covariance at P(k|k)
PSpred=zeros(m,K);  % predicted error surprise-based covariance at k|k
PSest=zeros(m,K);   % estimate error surprise-based covariance k|k-1 
Spred=zeros(m,K);   % predicted innovation covariance S(k|k-1)
Sest=zeros(m,K);    % estimated innovation covariance S(k|k)
nu=zeros(m,K);      % innovation vector
KG=zeros(n,K);      % Kalman Gain

% Performance test vectors
check=zeros(2,MC);   % marker to check consistancy at each monte-carlo 
Cnt=zeros(K,MC);    % C=R*inv(Spred)+ R*inv(PSpred) at each monte-carlo 
MahDis_MC=zeros(m,K,MC); % Mahalanobis distance at each monte-carlo 
Surprise_MC=zeros(m,K,MC); % Surprise distance at each monte-carlo 

AVG_inv_PS_MC=zeros(m,K,MC);   % inverse of the Sample-based Innovation Covariance at each monte-carlo
AVG_S_MC=zeros(m,K,MC);    % Sample-based Surprise Covariance at each monte-carlo
AVG_Surprise_MC=zeros(m,K,MC); % Sample-based Bayesian Surprise at each monte-carlo

%% Initial Conditions

x(:,1) =[1; 1];
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
z(:,1)=H*x(:,1);


for mc=1:MC
%% System Model

 Cnt(1,mc)=1;
 for i=2:K
    x(:,i)=F*x(:,i-1)+  mvnrnd([0 0],Q)'; % state equation
    z(:,i)=H*x(:,i)+ sqrt(R)*randn(1,1); % Measurement equation 
 end
 
%% Kalman Filter Estimation and Prediction
for j=2:K    
 [xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j)]=timeupdate_surprise(xest(:,j-1),Pest(:,:,j-1),PSest(:,j-1),Sest(:,j-1),R,H,Q,F); 
 [xest(:,j), Pest(:,:,j), PSest(:,j), Sest(:,j), KG(:,j),nu(:,j),Cnt(j,mc)] = measurementupdate_surpise(xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j),z(:,j),R,H); 
end

%% Tests for checking The consistancy of Kalman Filter

 [Surprise, Mov_Surprise, AVG_Surprise, AVG_inv_PS, AVG_S, SIGMA]=Test_Surprise_Final(K, m, nu, PSpred,R);

  Surprise_MC(:,:,mc)=Surprise;
  AVG_inv_PS_MC(:,:,mc)=AVG_inv_PS;
  AVG_S_MC(:,:,mc)=AVG_S;

  mc 
end


Total_AVG_inv_PS=mean(AVG_inv_PS_MC,3);
Total_AVG_S=mean(AVG_S_MC,3);
Total_Surprise=mean(Surprise_MC,3);

AVG_S_inv(1)=0;
for i=2:K
   AVG_S_inv(i)=pinv(Total_AVG_S(i));
end



f=2:K;
figure, 
p=plot(f, Total_AVG_inv_PS(:,f), 'r :',f, AVG_S_inv(:,f),'b --');
 p(1).LineWidth = 3;
 p(2).LineWidth = 3;
leg2 = legend('$\bar{\mathbf{P}}_{\mathcal{S}^B}(k|k-1)^{-1}$','$\bar{\mathbf{S}}(k|k-1)^{-1}$');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('\textbf{Sample-based Covariance}','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
axis([1 100 0 1])
hold on

figure, 
p=plot(1:K, Total_Surprise, 'r :');
 p(1).LineWidth = 3;
leg2 = legend('$\mathcal{S}^B_k$');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('\textbf{Distances}','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
hold on
