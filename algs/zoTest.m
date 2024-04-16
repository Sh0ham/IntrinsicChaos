function zoTest(data)
if size(data,2) == 1; data=data'; end
N=length(data); j=1:N; N0=round(N/10); Xi=1:N0;
D=zeros(1,N0); K=zeros(100,1); c=pi/5+rand(1,100)*3*pi/5;
for i=1:100
    p=cumsum(data.*cos(j*c(i))); q=cumsum(data.*sin(j*c(i)));
    for n=1:N0
        D(n)=mean((p(n+1:N)-p(1:N-n)).^2 + (q(n+1:N)-q(1:N-n)).^2) - ...
        mean(data)^2*((1-cos(n*c(i)))/(1-cos(c(i))));
    end
    K(i)=corr(Xi',D');
end
K = mean(K);
plot(p,q,'k.-'); set(gca,'FontSize',30);
xlabel('$p_c$','FontSize',64,'Interpreter','latex');
ylabel('$q_c$','FontSize',64,'Interpreter','latex');
legend(sprintf('$K_c=$ %.2f',K),'FontSize',48,'Interpreter','latex');
end