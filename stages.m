% stages
alpha = 0.4;

figure
sc = 29;
h1 = plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
hold on 
h2 = plot(data(sc).Pss(:,1,1),data(sc).Pss(:,2,1),'-', 'color',[1,0,0], 'linewidth', 3.5);
h3 = plot(data(sc).Pss(:,1,3),data(sc).Pss(:,2,3),'-', 'color',[0,0,0,alpha], 'linewidth', 3);
%axis([-0.03 0.6 -0.08 0.09])
samples = linspace(1,size(data(sc).Pss,3),4);
for i = 2:4
% plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
% hold on 
% plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'g',  'linewidth', 2.8);
% alpha = alpha+0.15;
plot(data(sc).Pss(:,1,round(samples(i))),data(sc).Pss(:,2,round(samples(i))),'-', 'color',[0,0,0,alpha], 'linewidth', 3);
%axis([-0.03 0.6 -0.08 0.09])
%axis([-0.02 0.6 -0.01 0.21])
%axis([-0.1 0.65 -0.12 0.23])
% axis([-0.02 0.7 -0.15 0.1])
ylabel('Y [m]')
xlabel('X [m]')
% hold off
end
% plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
% hold on 
h4 = plot(data(sc).Pss(:,1,end),data(sc).Pss(:,2,end),'k-',  'linewidth', 3.5);

legend([h1,h2,h3, h4], 'Desired shape', 'Initial shape', 'Intermediary stages', 'Final shape')
legend show
% axis([-0.02 0.7 -0.15 0.1])
%axis([-0.03 0.6 -0.08 0.09])
%axis([-0.02 0.6 -0.01 0.21])
% axis([-0.1 0.65 -0.12 0.23])

ylabel('Y [m]')
xlabel('X [m]')