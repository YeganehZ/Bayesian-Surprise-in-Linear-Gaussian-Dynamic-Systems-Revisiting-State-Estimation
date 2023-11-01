function [Xest, Pest, PSest, Sest, KG,nu,C] = measurementupdate_surpise(Xpred,Ppred,PSpred,Spred,z,R,H)

KG = Ppred*H'*pinv(Spred); %% Kalman gain 
nu= z - H * Xpred;
C=R*pinv(Spred)+ R*pinv(PSpred);
Xest= Xpred + KG * C* nu;
Pest= Ppred - KG *C* Spred * KG';
PSest= PSpred + R ;
Sest= Spred - H*KG*(Spred)*KG'*H';