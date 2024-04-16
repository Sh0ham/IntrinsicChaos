function dydt = odeSys(t,y,bG2)
aG1=12; bG1=30; cG1=1; kG1=2;  mG1=2;
aG2=17;         cG2=1; kG2=60; mG2=6;
d=0.3;          z=90;  N_in=20; 

betaG1=(aG1*y(4)^bG1)/(kG1+y(4)^bG1);
betaG2=(aG2*y(4)^bG2)/(kG2+y(4)^bG2);    
muG1=(mG1*y(4))/(cG1+y(4));
muG2=(mG2*y(4))/(cG2+y(4));

dydt=zeros(4,1);
dydt(1)=-betaG1*y(1)+2*z*y(3)-d*y(1);
dydt(2)=betaG1*y(1)-betaG2*y(2)-d*y(2);
dydt(3)=betaG2*y(2)-z*y(3)-d*y(3);
dydt(4)=-muG1*y(1)-muG2*y(2)+d*N_in-d*y(4);
end