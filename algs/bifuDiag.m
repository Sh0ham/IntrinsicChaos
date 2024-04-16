function bifuDiag
dots=[]; time=[0 300]; IC=[1;2;0;0.1]; bG2=1;
for i=1:1500
    [t,y]=ode45(@odeSys, time, IC, [], bG2);
    obj=y(:,1)+y(:,2)+y(:,3);
    for j=round(length(obj)/2):length(obj)-1
        if((obj(j)>=obj(j-1))&&(obj(j)>=obj(j+1)))...
        ||((obj(j)<=obj(j-1))&&(obj(j)<=obj(j+1)))
            dots=[dots;i obj(j)];
        end
    end
    bG2=bG2+0.04;
end
plot(dots(:,1)/25,dots(:,2),'k.','MarkerSize',2);
xlabel('$b_{G_2}$','FontSize',48,'Interpreter','latex');
ylabel('$Min/Max\;Abundance\left(ind.\mu^{-1}\right)$',
       'FontSize',48,'Interpreter','latex');
end