  close all
    clear 
    clc
    

    
    % Cable length:
    L = 0.7;

    % Cable points:
    N = 6;

    % Segment length
    l_s = L/(N-1);


shape = ["QSW", "I", "U", "QSW"; "HSW1", "S", "QSW", "L"];
for i = 1:size(shape,2)
    [profs,xy] = des_profile(shape(1,i), N, l_s, L);
    P = profs(:,1:2) ;
    P(:,2) = -P(:,2);
    xy(:,2) = -xy(:,2); 
    
    [des_pos, des_xy] = des_profile(shape(2,i), N, l_s, L);
    des_P = des_pos(:,1:2) +[0.100 0.0000]*1;
    des_xy = des_xy+[0.100 0.0000]*1;
%     if i == 2
%         P(:,2) =P(:,2); 
%         P = P +0*[-0.3 0];
%         
%     end
    profiles = inter_proiles(P, des_P, l_s);
    figure
    hold on
    title(append('Case',num2str(i)))
%     h1 =  plot(profiles(:,1,1), profiles(:,2,1), 'Color',[0.5, 0.5, 1, 0.4], 'LineWidth', 2);
%     for j= 2:size(profiles,3)
% 
%         plot(profiles(:,1,j), profiles(:,2,j), 'Color',[0.5, 0.5, 1, 0.4], 'LineWidth', 2)
%         
%     end
    h2 = plot(xy(:,1),xy(:,2),'r',  'LineWidth', 3);
    
    h3 = plot(des_xy(:,1),des_xy(:,2), 'g-', 'LineWidth', 3);
    handlevec = [h2 h3];
    legend(handlevec,'Initial shape','Desired shape');
   
%    legend('Initial shape','Desired shape', 'Intermediary profiles')
legend show
    yt = get(gca, 'YTick');
    xt = get(gca, 'XTick');
    xlabel('X [m]')
    ylabel('Y [m]')
%     axis([min(xt(1),yt(1)) max(xt(end),yt(end)) min(xt(1),yt(1)) max(xt(end),yt(end))])
   

end
