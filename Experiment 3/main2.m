clc
clear all

K=500; % Time index
MC=1;  % Number of Monte Carlo runs

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
check=zeros(2,1);   % marker to check consistancy 
Cnt=ones(K,1);     % C=R*inv(Spred)+ R*inv(PSpred)  

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


%% System Model
 for i=2:K
    x(:,i)=F*x(:,i-1)+  mvnrnd([0 0],Q)'; % state equation
    z(:,i)=H*x(:,i)+ sqrt(R)*randn(1,1); % Measurement equation 
 end
 
%% Kalman Filter Estimation and Prediction

 for j=2:K    
  [xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j)]=timeupdate_surprise(xest(:,j-1),Pest(:,:,j-1),PSest(:,j-1),Sest(:,j-1),R,H,Q,F); 
  [xest(:,j), Pest(:,:,j), PSest(:,j), Sest(:,j), KG(:,j),nu(:,j),Cnt(j,1)] = measurementupdate_surpise(xpred(:,j),Ppred(:,:,j),PSpred(:,j),Spred(:,j),z(:,j),R,H);     
end

%% Tests for checking The consistancy of Kalman Filter

[Surprise, Mov_Surprise, AVG_Surprise, AVG_inv_PS, AVG_S, SIGMA]=Test_Surprise_Final(K, m, nu, PSpred,R);

%% Plot the consistancy test results
figure, 
p=plot(1:K, Surprise, 'k -',1:K, Mov_Surprise,'b --', 1:K, AVG_Surprise, 'r -.', 1:K, 0.5*SIGMA, 'm :' );
 p(1).LineWidth = 2;
 p(2).LineWidth = 3;
 p(3).LineWidth = 3;
 p(4).LineWidth = 3;
leg2 = legend('$\mathcal{S}^{B}_k$','$\bar{\bar{\mathcal{S}}}^{B}_{k}$','$\bar{\mathcal{S}}^{B}_k$','$\textbf{E}[\mathcal{S}^B_{k}]$');
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',22);
xlabel({'\textbf{Time Samples($k$)}'},'Interpreter','latex','FontSize',22,'Color','k'); 
ylabel('\textbf{Bayesian Surprise} ($\mathcal{S}^{B}_k$)','Interpreter','latex','FontSize',22,'Color','k');
set(gca,'FontSize',22)
set(gca,'FontName','Times New Roman')
set(gcf,'color','white')
axis([50 500 0 6])
% grid on;
hold on

