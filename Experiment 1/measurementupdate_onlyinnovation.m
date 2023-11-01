function [Xest, Pest, PSest, Sest, KG,nu,C] = measurementupdate_onlyinnovation(Xpred,Ppred,PSpred,Spred,z,R,H)

KG = Ppred*H'*pinv(Spred); %% Kalman gain 
nu= z - H * Xpred;
C=R*pinv(Spred);
Xest= Xpred + KG * C* nu;
Pest= Ppred - KG *C*Spred* KG';
PSest= PSpred + R ;
Sest= Spred - H*KG*C*Spred*KG'*H';