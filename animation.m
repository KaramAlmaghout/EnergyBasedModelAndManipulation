% animation

figure
% h1 = plot(data(4).desP(:,1),data(4).desP(:,2),'g',  'linewidth', 4);
% hold on 
% h2 = plot(data(4).Pss(:,1,1),data(4).Pss(:,2,1),'-', 'color',[1,0,0], 'linewidth', 2);
% h3 = plot(data(4).Pss(:,1,2),data(4).Pss(:,2,2),'-', 'color',[1,0,0,0.5], 'linewidth', 2);
%axis([-0.03 0.6 -0.08 0.09])
sc = 4
for i = 1:1:9
plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
hold on 
% plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'g',  'linewidth', 2.8);
plot(data(sc).Pss(:,1,i),data(sc).Pss(:,2,i),'-', 'color',[1,0,0,0.5], 'linewidth', 3.5);
%axis([-0.03 0.6 -0.08 0.09])
%axis([-0.02 0.6 -0.01 0.21])
%axis([-0.1 0.65 -0.12 0.23])
axis([-0.02 0.7 -0.15 0.1])
pause(0.4)
ylabel('Y [m]')
xlabel('X [m]')
hold off
end
plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
hold on 
plot(data(sc).Pss(:,1,end),data(sc).Pss(:,2,end),'r-',  'linewidth', 3.5);

% legend([h1,h2,h3, h4], 'Desired shape', 'Initial shape', 'Intermediary stages', 'Final shape')
% legend show
axis([-0.02 0.7 -0.15 0.1])
%axis([-0.03 0.6 -0.08 0.09])
%axis([-0.02 0.6 -0.01 0.21])
axis([-0.1 0.65 -0.12 0.23])

ylabel('Y [m]')
xlabel('X [m]')