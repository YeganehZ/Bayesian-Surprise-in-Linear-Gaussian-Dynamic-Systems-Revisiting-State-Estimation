function  [Xpred,Ppred,PSpred,Spred,A,B]=timeupdate_surprise(Xest,Pest,PSest,Sest,R,H,Q,F)

Xpred = F*Xest; 
Ppred = F*Pest*F' + Q;
Spred=  R + H*Ppred*H';

A=H'*pinv(Sest)*PSest*inv(R)*H;
B=H*(Q + F*pinv(A)*F')*H';
PSpred = Spred * pinv(B) * R;


