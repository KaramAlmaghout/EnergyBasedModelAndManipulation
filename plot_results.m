ln = 8;
N = 6:2:20;
legends = {'e_{mean}', 'e_{std}', 'e_{max}'};
for j = 1:4
    figure
    x = zeros(1,ln);
    y = zeros(ln,2);
    z = zeros(ln,1);
    ynull = zeros(ln,1);
    znull = zeros(ln,2);
    c=j;
    k = 1;
    n=0;
    for i = 1+(c*ln-ln):c*ln
%         i = i+1;
        n = n+1;
        x(k) = N(n);
        y(k,:) = [ data(i).error_mean(end) data(i).error_std(end)]*1000;
        z(k,1) = data(i).error_max(end)*1000;
        k  = k+1;
    end
%     data(i).error_max(end)
    b = bar(x,[y ynull],'grouped');
    axis([7 17 -inf 12])
    ylabel('[mm]');
%     legend(b,'e_{mean}', 'e_{std}')
    yyaxis right
    b1 = bar(x,[znull z],'grouped');
    ylabel('[mm]');
    xlabel('N')
%     legend(b,'e_{max}')
    legend([b(1) b(2) b1],'e_{mean}', 'e_{std}','e_{max}');    
    title(append('case',num2str(j)))
    hold on
    xlim=get(gca,'xlim');
%     ylim([0 0.015]);
   
    
%     legend(b,'e_{max}','e_{mean}', 'e_{std}')
%     axes('Position',[0.3 0.5 0.4 0.4])
%     box on
% % %     axis off
%     plot(data(j*4).xy(:,1),data(j*4).xy(:,2), 'color', [0,1,0,0.4], 'linewidth', 4)
%     set(gca,'visible','off')
grid on
end

% initial and final
for i = 8:8:32
figure
plot(data(i).Pss(:,1,1),data(i).Pss(:,2,1),'r',  'linewidth', 2.8);
hold on 
plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'g',  'linewidth', 2.8);
plot(data(i).desP(:,1),data(i).desP(:,2),'k',  'linewidth', 2.8);
legend('Initial shape', 'Final shape')
hold off
ylabel(['Y [m]'])
xlabel(['X [m]'])
end

% initial and final

for i = 1:8:32
figure
plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'r*',  'linewidth', 2.8);
hold on 
% plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'g',  'linewidth', 2.8);
plot(data(i).desP(:,1),data(i).desP(:,2),'g',  'linewidth', 2.8);
legend('Initial shape', 'Final shape')
hold off
ylabel(['Y [m]'])
xlabel(['X [m]'])
end


for k = 1:8:32
    figure
    plot(data(k+7).desP(:,1),data(k+7).desP(:,2), 'g-', 'linewidth', 4)
    hold on
%     plot(data(k+5).P(:,1),data(k+5).P(:,2), '*-', 'linewidth', 3.6)
%     plot(data(k+4).P(:,1),data(k+4).P(:,2), '*-', 'linewidth', 3.2)
    plot(data(k+5).Pss(:,1,end),data(k+5).Pss(:,2,end), '*-', 'linewidth', 3.2)
    plot(data(k+4).Pss(:,1,end),data(k+4).Pss(:,2,end), '*-', 'linewidth', 2.8)
    plot(data(k+3).Pss(:,1,end),data(k+3).Pss(:,2,end), '*-', 'linewidth', 2.4)
    plot(data(k+2).Pss(:,1,end),data(k+2).Pss(:,2,end), '*-', 'linewidth', 2)
    plot(data(k+1).Pss(:,1,end),data(k+1).Pss(:,2,end), '*-', 'linewidth', 1.6)
    legend('desired shape', 'N = 16','N = 14','N = 12','N = 10','N = 8')
ylabel('Y [m]')
xlabel('X [m]')
end

figure
sc = 4;
h1 = plot(data(sc).desP(:,1),data(sc).desP(:,2),'g',  'linewidth', 4);
hold on 
h2 = plot(data(sc).Pss(:,1,1),data(sc).Pss(:,2,1),'-', 'color',[1,0,0], 'linewidth', 2);
h3 = plot(data(sc).Pss(:,1,2),data(sc).Pss(:,2,2),'-', 'color',[1,0,0,0.5], 'linewidth', 2);

for i = 4:2:size(data(sc).Pss,3)

% plot(data(i).Pss(:,1,end),data(i).Pss(:,2,end),'g',  'linewidth', 2.8);
plot(data(sc).Pss(:,1,i),data(sc).Pss(:,2,i),'-', 'color',[1,0,0,0.5], 'linewidth', 2);

end
h4 = plot(data(sc).Pss(:,1,end),data(sc).Pss(:,2,end),'k-',  'linewidth', 2);

legend([h1,h2,h3, h4], 'Desired shape', 'Initial shape', 'Intermediary stages', 'Final shape')
legend show
hold off
ylabel('Y [m]')
xlabel('X [m]')

figure
for k = 1:8:16
    x = [6, 8, 10, 12, 14, 16, 18, 20];
% %     y = [data(k).iter data(k+1).iter , data(k+2).iter, data(k+3).iter ];
    y = [data(k).time data(k+1).time , data(k+2).time, data(k+3).time , data(k+4).time, data(k+5).time , data(k+6).time, data(k+7).time];
    plot(x,y, '*-', 'linewidth', 3)
    hold on
end
legend('Scenario1','Scenario2')
xlabel('Cable feature points N')
ylabel('Time [sec]')