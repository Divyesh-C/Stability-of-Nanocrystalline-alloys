close all;
%Inputting all the constants involved in calculating the delG value of the mixture.
z=12;
v=0.5;
delHmix=0;
delHseg=0;
Oa=0.0000014*(10^27);
Ga=0.26*(10^-18);
Ob=0.0000481*(10^27);
Gb=0.93*(10^-18);
t=0.5;
k=8.314;
T=300;

%function [delGM, D, X] =CalcFree(z, v, delHmix, delHseg, Oa, Ga, Ob, Gb, t, k, T)
X=zeros(1,1);
D=zeros(1,1);
delGM=zeros(1,1);
    for d=1:1:100
        for Xgb=0.1*(10^-27):0.05*(10^-27):0.9*(10^-27)
                Fgb=1-(((d-t)/d)^3);
                X1=0.1*(10^-27);
                Xc=(X1*(10^-27)-Fgb*Xgb*(10^-27))/(1-Fgb);
                if Xc>0
                        Oc=delHmix/(z*X1*(1-X1));
                        Ogb=2*(Oc-delHseg/z);
                        DelGc=z*Oc*Xc*(1-Xc) + k*T*(Xc*log(Xc)+(1-Xc)*log(1-Xc));
                        DelGgb=z*Ogb*Xgb*(1-Xgb) + Oa*Ga*(1-Xgb)/t + Ob*Gb*Xgb/t + k*T*(Xgb*log(Xgb)+(1-Xgb)*log(1-Xgb));
                        DelGm= (1-Fgb)*DelGc + Fgb*DelGgb + z*v*Fgb*(Xgb-Xc) * ((2*Xgb-1)*Ogb - (Oa*Ga - Ob*Gb)/(z*t));
                        X(end+1,1)=Xgb;
                        D(end+1,1)=d;
                        delGM(end+1,1)=DelGm; 
                        
               end
        end
    end
%L1=length(D);
D=D(2:length(D));
%L2=length(D);
X=X(2:length(X));
%L3=length(X);
delGM=delGM(2:length(delGM));
%L4=length(delGM);
 [delmin,I]=min(real(delGM(:)));
 x=X(I);
 y=D(I);
 disp(delmin);
 disp(x);
 disp(y);
 figure;
 T=table(X,D,real(delGM));
 %disp(T)
 scatter3(X,D,real(delGM));
 title('Free Energy Plot');
 xlabel('Solute Concentration'); 
 ylabel('Grain Size');
 zlabel('Free Energy'); 
%end



