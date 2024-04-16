function chemostatModel(bG2,t0,tf)
time=[t0 tf]; IC=[1; 2; 0; 0.1];
[t,y]=ode45(@odeSys, time, IC, [], bG2);
abundance=y(:,1) + y(:,2) + y(:,3);

plot(t,y(:,1),'k-');
hold on;
plot(t,y(:,2),'r-');
plot(t,y(:,3),'-','Color','#7E2F8E');
plot(t,y(:,4),'-','Color','#0072BD')
plot(t,abundance,'-','Color','#D95319');
axis([0 60 0 3]);
title('$Chemostat\;Model$','FontSize',48,'Interpreter','latex');
legend('$G_{1}$','$G_{2}$','$M$','$N$','$Total$',
       'FontSize',32,'Interpreter','latex');
xlabel('$Time$','FontSize',40,'Interpreter','latex');
ylabel('$Abundance\left(ind.\mu^{-1}\right)$',
       'FontSize',40,'Interpreter','latex');

plot3(y(2500:end,1),y(2500:end,2),y(2500:end,3), 'k');
title('$Attractor$','FontSize',48,'Interpreter','latex');
xlabel('$M\;\left(ind.\mu^{-1}\right)$','FontSize',40,'Interpreter','latex');
ylabel('$G_{2}\;\left(ind.\mu^{-1}\right)$','FontSize',40,'Interpreter','latex');
zlabel('$G_{1}\;\left(ind.\mu^{-1}\right)$','FontSize',40,'Interpreter','latex');
end