function LLE
bG2=1; j=1;
for i=1:600
    lambdas(j,1:5) = lyapExp(bG2);
    bG2=bG2+0.1;
    j=j+1;
end

figure(1)
xlim([1 i/100]); hold on; yline(0,'k--')
plot((1:i)/100,lambdas(1:end,1:1),'k.-');
plot((1:i)/100,lambdas(1:end,2:2),'r.-');
plot((1:i)/100,lambdas(1:end,3:3),'.-','Color','#7E2F8E');
legend('','$\lambda_{1}$','$\lambda_{2}$','$\lambda_{3}$',
       'Location','best','FontSize',32,'Interpreter','latex');
xlabel('$b_{G_{2}}$','FontSize',48,'Interpreter','latex');
ylabel('$\lambda_{i}$','FontSize',48,'Interpreter','latex');
title('Dynamics of Lyapunov Exponents Spectrum','FontSize',48,
      'Interpreter','latex');

figure(2)
xlim([1 i/100]); hold on; yline(0,'k--')
set(gca,'FontSize',25);
plot((1:i)/100,lambdas(1:end,1:1),'k.-');
plot((1:i)/100,lambdas(1:end,2:2),'r.-');
plot((1:i)/100,lambdas(1:end,3:3),'.-','Color','#7E2F8E');
plot((1:i)/100,lambdas(1:end,4:4),'.-','Color','#0072BD');
legend('','$\lambda_{1}$','$\lambda_{2}$','$\lambda_{3}$','$\lambda_{4}$',
       'Location','best','FontSize',32,'Interpreter','latex');
xlabel('$b_{G_{2}}$','FontSize',48,'Interpreter','latex');
ylabel('$\lambda_{i}$','FontSize',48,'Interpreter','latex');
title('Dynamics of Lyapunov Exponents Spectrum','FontSize',48,
      'Interpreter','latex');

figure(3)
lmax = lambdas(1:end,5:5);
xlim([1 i/100]); hold on; yline(0,'k--')
set(gca,'FontSize',25);
LMAX_RED = lmax; LMAX_RED(lmax<=0) = NaN;
plot((1:i)/100,lmax,'k.-');
plot((1:i)/100,LMAX_RED,'.-','Color','#EA362A');
xlabel('$b_{G_{2}}$','FontSize',48,'Interpreter','latex');
ylabel('$\lambda$','FontSize',48,'Interpreter','latex');
title('Maximal Lyapunov Exponent','FontSize',48,'Interpreter','latex');
end