function expo = lyapExp(bG2)
aG1=12; bG1=30; cG1=1; kG1=2;  mG1=2;
aG2=17;         cG2=1; kG2=60; mG2=6;
d=0.3;          z=90;  N_in=20; 
expo=zeros(1,5);
h=1e-3; h_norm=10*h; N=h_norm/h; tn=300-h; t=0:h:tn; n=length(t); T=0:h_norm:tn;
Q=1; q1=Q;q2=Q;q3=Q;q4=Q; %Order
cp1=1; cp2=1; cp3=1; cp4=1; %Fractional-order binomial coefficient

for j=1:n
    c1(j)=(1-(1+q1)/j)*cp1;
    c2(j)=(1-(1+q2)/j)*cp2;
    c3(j)=(1-(1+q3)/j)*cp3;
    c4(j)=(1-(1+q4)/j)*cp4;
    cp1=c1(j); cp2=c2(j); cp3=c3(j); cp4=c4(j);
end

x1(1)=1; x2(1)=2; x3(1)=0; x4(1)=0.1; k=1; %Initialization
f11(k)=1; f12(k)=0; f13(k)=0; f14(k)=0;
f21(k)=0; f22(k)=1; f23(k)=0; f24(k)=0;
f31(k)=0; f32(k)=0; f33(k)=1; f34(k)=0;
f41(k)=0; f42(k)=0; f43(k)=0; f44(k)=1;

J=[f11(k),f12(k),f13(k),f14(k);...
   f21(k),f22(k),f23(k),f24(k);...
   f31(k),f32(k),f33(k),f34(k);...
   f41(k),f42(k),f43(k),f44(k)];

SUM=[0;0;0;0;]; LE=[];

for k=2:n
    x1(k)=((-aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1))*x1(k-1) +
          2*z*x3(k-1)-d*x1(k-1))*h^q1-calmem(x1,c1,k);
    x2(k)=((aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1))*x1(k-1) -
           (aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2))*x2(k-1) -
           d*x2(k-1))*h^q2-calmem(x2,c2,k);
    x3(k)=((aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2))*x2(k-1) -
           z*x3(k-1) - d*x3(k-1))*h^q3-calmem(x3,c3,k);
    x4(k)=((-mG1*x4(k-1))/(cG1+x4(k-1))*x1(k-1) -
           (mG2*x4(k-1))/(cG1+x4(k-1))*x2(k-1) -
           d*x4(k-1) + d*N_in)*h^q4-calmem(x4,c4,k);

    f11(k)=((((-aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f11(k-1)) -
           (x1(k-1)*f41(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) + (2*z*f31(k-1)) -
           (d*f11(k-1)))*h^q1-calmem(f11,c1,k);
    f12(k)=((((-aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f12(k-1)) -
           (x1(k-1)*f42(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) + (2*z*f32(k-1)) -
           (d*f12(k-1)))*h^q1-calmem(f12,c1,k);
    f13(k)=((((-aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f13(k-1)) -
           (x1(k-1)*f43(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) + (2*z*f33(k-1)) -
           (d*f13(k-1)))*h^q1-calmem(f13,c1,k);
    f14(k)=((((-aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f14(k-1)) -
           (x1(k-1)*f44(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) +  (2*z*f34(k-1)) -
           (d*f14(k-1)))*h^q1-calmem(f14,c1,k);

    f21(k)=((((aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f11(k-1)) +
           (x1(k-1)*f41(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) - (((aG2*(x4(k-1)^bG2))/
           (kG2+(x4(k-1)^bG2)))*f21(k-1)) - (x2(k-1)*f41(k-1)*
           ((aG2*bG2*kG2*x4(k-1)^(bG2-1))/((kG2+x4(k-1)^bG2)^2))) -
           (d*f21(k-1)))*h^q2-calmem(f21,c2,k);
    f22(k)=((((aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f12(k-1)) +
           (x1(k-1)*f42(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) - (((aG2*(x4(k-1)^bG2))/
           (kG2+(x4(k-1)^bG2)))*f22(k-1)) - (x2(k-1)*f42(k-1)*
           ((aG2*bG2*kG2*x4(k-1)^(bG2-1))/((kG2+x4(k-1)^bG2)^2))) -
           (d*f22(k-1)))*h^q2-calmem(f22,c2,k);
    f23(k)=((((aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f13(k-1)) +
           (x1(k-1)*f43(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) - (((aG2*(x4(k-1)^bG2))/
           (kG2+(x4(k-1)^bG2)))*f23(k-1)) - (x2(k-1)*f43(k-1)*
           ((aG2*bG2*kG2*x4(k-1)^(bG2-1))/((kG2+x4(k-1)^bG2)^2))) -
           (d*f23(k-1)))*h^q2-calmem(f23,c2,k);
    f24(k)=((((aG1*(x4(k-1)^bG1))/(kG1+(x4(k-1)^bG1)))*f14(k-1)) +
           (x1(k-1)*f44(k-1)*((aG1*bG1*kG1*x4(k-1)^(bG1-1))/
           ((kG1+x4(k-1)^bG1)^2))) - (((aG2*(x4(k-1)^bG2))/
           (kG2+(x4(k-1)^bG2)))*f24(k-1)) - (x2(k-1)*f44(k-1)*
           ((aG2*bG2*kG2*x4(k-1)^(bG2-1))/((kG2+x4(k-1)^bG2)^2))) -
           (d*f24(k-1)))*h^q2-calmem(f24,c2,k);

    f31(k)=((((aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2)))*f21(k-1)) +
           (x2(k-1)*f41(k-1)*((aG2*bG2*kG2*x4(k-1)^(bG2-1))/
           ((kG2+x4(k-1)^bG2)^2))) - (z*f31(k-1)) -
           (d*f31(k-1)))*h^q3-calmem(f31,c3,k);
    f32(k)=((((aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2)))*f22(k-1)) +
           (x2(k-1)*f42(k-1)*((aG2*bG2*kG2*x4(k-1)^(bG2-1))/
           ((kG2+x4(k-1)^bG2)^2))) - (z*f32(k-1)) -
           (d*f32(k-1)))*h^q3-calmem(f32,c3,k);
    f33(k)=((((aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2)))*f23(k-1)) +
           (x2(k-1)*f43(k-1)*((aG2*bG2*kG2*x4(k-1)^(bG2-1))/
           ((kG2+x4(k-1)^bG2)^2))) - (z*f33(k-1)) -
           (d*f33(k-1)))*h^q3-calmem(f33,c3,k);
    f34(k)=((((aG2*(x4(k-1)^bG2))/(kG2+(x4(k-1)^bG2)))*f24(k-1)) +
           (x2(k-1)*f44(k-1)*((aG2*bG2*kG2*x4(k-1)^(bG2-1))/
           ((kG2+x4(k-1)^bG2)^2))) - (z*f34(k-1)) -
           (d*f34(k-1)))*h^q3-calmem(f34,c3,k);

    f41(k)=((((-mG1*x4(k-1))/(cG1+x4(k-1)))*f11(k-1)) -
           (x1(k-1)*f41(k-1)*((cG1*mG1)/((cG1+x4(k-1))^2))) -
           (((mG2*x4(k-1))/(cG2+x4(k-1)))*f21(k-1)) -
           (x2(k-1)*f41(k-1)*((cG2*mG2)/((cG2+x4(k-1))^2))) -
           (d*f41(k-1)))*h^q4-calmem(f41,c4,k);
    f42(k)=((((-mG1*x4(k-1))/(cG1+x4(k-1)))*f12(k-1)) -
           (x1(k-1)*f42(k-1)*((cG1*mG1)/((cG1+x4(k-1))^2))) -
           (((mG2*x4(k-1))/(cG2+x4(k-1)))*f22(k-1)) -
           (x2(k-1)*f42(k-1)*((cG2*mG2)/((cG2+x4(k-1))^2))) -
           (d*f42(k-1)))*h^q4-calmem(f42,c4,k);
    f43(k)=((((-mG1*x4(k-1))/(cG1+x4(k-1)))*f13(k-1)) -
           (x1(k-1)*f43(k-1)*((cG1*mG1)/((cG1+x4(k-1))^2))) -
           (((mG2*x4(k-1))/(cG2+x4(k-1)))*f23(k-1)) -
           (x2(k-1)*f43(k-1)*((cG2*mG2)/((cG2+x4(k-1))^2))) -
           (d*f43(k-1)))*h^q4-calmem(f43,c4,k);
    f44(k)=((((-mG1*x4(k-1))/(cG1+x4(k-1)))*f14(k-1)) -
           (x1(k-1)*f44(k-1)*((cG1*mG1)/((cG1+x4(k-1))^2))) -
           (((mG2*x4(k-1))/(cG2+x4(k-1)))*f24(k-1)) -
           (x2(k-1)*f44(k-1)*((cG2*mG2)/((cG2+x4(k-1))^2))) -
           (d*f44(k-1)))*h^q4-calmem(f44,c4,k);

    J=[f11(k),f12(k),f13(k),f14(k);...
       f21(k),f22(k),f23(k),f24(k);...
       f31(k),f32(k),f33(k),f34(k);...
       f41(k),f42(k),f43(k),f44(k)];
    
    if mod(k,N)==0
        [J,E]=GSR(J);
        SUM=SUM+log(E);
        f11(k)=J(1,1); f12(k)=J(1,2); f13(k)=J(1,3); f14(k)=J(1,4);
        f21(k)=J(2,1); f22(k)=J(2,2); f23(k)=J(2,3); f24(k)=J(2,4);
        f31(k)=J(3,1); f32(k)=J(3,2); f33(k)=J(3,3); f34(k)=J(3,4);
        f41(k)=J(4,1); f42(k)=J(4,2); f43(k)=J(4,3); f44(k)=J(4,4);
        LE=[LE,SUM/(k*h)];
        fprintf('bG2 = %f, Completed: %f%%\n',bG2,k/n*100)
    end
end
expo(1:4) = LE(end-3:end,end:end);
expo(5:5) = max(LE(:,end));
end