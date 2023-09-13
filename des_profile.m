function [pts, xy] = des_profile(P,N,l_s,L)


ratio = 1;
switch P
    case "I"
         x = linspace(0,1,1000001);
%         y = 0.17*sin(2*x);
         y = zeros(1,length(x));
        
    case "S"
        R = [cosd(60) -sind(60);sind(60) cosd(60)];
% R = [cosd(0) -sind(0);sind(0) cosd(0)];
        L0 =  1.2347;
        ratio = L/L0+0.05;
%         ratio = 1;
        x = linspace(0,pi/(3/ratio),1000001);
%         y = 0.17*sin(2*x);
        y = 0.15*ratio*sin((6/ratio)*x);
        y = -y;
        xy = [x' y']*R';
        x = xy(:,1)';
        y = xy(:,2)';
    case "U"
        L0 = 0.7054;
        ratio = L/L0+0.05;
        x = linspace(0,pi/(7/ratio),1000001);
        y = 0.256*ratio*sin((7/ratio)*x);
        y = -y;
    case "HSW1"
        L0 = 1.0704;
%         ratio = L/L0+0.005;
ratio = L/L0;
        x = linspace(0,pi/(3/ratio),1000001);
        y = 0.1*ratio*sin((3/ratio)*x);
    case "QSW"
        L0 = 1.0704;
%         ratio = L/L0+0.005;
        ratio = L/L0;
        x = linspace(0,pi/(3/ratio),1000001);
        y = 0.2*ratio*sin((1.5/ratio)*x);
    case "TPL"
        
        L0 = 2;
        ratio = L/L0+0.005;
        x = [linspace(0,0,400) linspace(0,0.5,601) linspace(0.5,0.5,400)]*ratio ;
        y = [linspace(0.25,-0.25,400) linspace(0.25,0.25,601) linspace(0.25,-0.5,400)]*ratio;
    case "L1"
        
        L0 = 0.7523;
        ratio = L/L0+0.005;
        x = [linspace(0,0,400) linspace(0,0.7,601)]*ratio ;
        y = [linspace(0.25,-0.2,400) linspace(-0.2,-0.2,601)]*ratio;
    case "L"
        L0 =  9.5708;
        ratio = L/(L0-.5);
%         ratio = 1;
        x1 = -1*ones(1,5000);
        y1 = linspace(4,0,5000);
        r = 1;
        th = pi:pi/5000:3*pi/2;
        x2 = r * cos(th) + 0;
        y2 = r * sin(th) + 0;
        x3 = linspace(x2(end), x2(end)+4, 5000);
        y3 = y2(end)*ones(1,5000);
        x = [x1 x2(2:end) x3(2:end)]*ratio;
        y = [y1 y2(2:end) y3(2:end)]*ratio;
        
        R = [cosd(10) -sind(10);sind(10) cosd(10)];
% R = [cosd(0) -sind(0);sind(0) cosd(0)];
        xy = [x' y']*R';
        x = xy(:,1)';
        y = xy(:,2)';
    case "Sine"
        L0 = 11.3385;
        ratio = 0.05+L/L0;
        x = linspace(0.7,2*pi,1000001)*ratio;
        y = coth(x)*ratio;
%         x = x*ratio;
    case "AFSW"
        L0 = 4.9348;
        ratio = L/L0+0.015;
%         ratio = 1;
        x1 = linspace(0,2*pi/(2/ratio),1000001);
        y1 = 0.2*ratio*sin((3/ratio)*x1);
        x2 = linspace(0,2*pi/(1/ratio),1000001);
        y2 = 0.15*ratio*sin((1.5/ratio)*x2);
        x = [x1(1:500001) x2(500002:end)-x2(500001)+x1(500001)];
        y = [y1(1:500001) y2(500002:end)];
    case "M"
        
        L0 = 16.2052;
        ratio = L/L0;
%         ratio = 1;
        x = linspace(0,(5*pi),1000001);
        y = 0.6*sin(0.6*x);
        x = x*ratio;
        y = y*ratio;
%         x2 = linspace(0,2*pi/(5/ratio),1000001);
%         y2 = 0.15*ratio*sin((7.5/ratio)*x2);
%         x = [x1(1:500001) x2(500002:end)-x2(500001)+x1(500001)];
%         y = [y1(1:500001) y2(500002:end)];
    case "Zshape"
        L0 =  3.7788;
        ratio = L/L0+0.005;
        y1 = linspace(0,0.7*ratio,100000);
        x1 = y1/2;
        y2 = linspace(y1(end),0,100000);
        x2 = x1(end) + 3*(y2(1)-y2);
        y3 = linspace(y2(end),0.7*ratio,100000);
        x3 = x2(end)+y3/2;
        x = [x1(1:end-1) x2(1:end-1) x3];
       y = [y1(1:end-1) y2(1:end-1) y3];
    case "Z"
        L0 =  2.3899;
        ratio = L/L0+0.005;
        y1 = linspace(0.7,0.7,100000)*ratio;
        x1 = linspace(0,0.7,100000)*ratio;
        y2 = linspace(0.7,0,100000)*ratio;
        x2 = linspace(0.7,0,100000)*ratio;
        y3 = linspace(0,0,100000)*ratio;
        x3 = linspace(0,0.7,100000)*ratio;
        x = [x1(1:end-1) x2(1:end-1) x3];
        y = [y1(1:end-1) y2(1:end-1) y3];
end


%%%%%%%
pts = zeros(N,2);
pts(1,:) = [x(1) y(1)];
j = 1;

for i = 2:length(x)
    p1 = pts(j,:);
    p2 = [x(i), y(i)];
    if abs((norm(p2-p1)-l_s)) < 0.0005
        j = j + 1;
        pts(j,:) = p2;
%         p1 = pts(j,:);
    end
    if j == N
        break
    end
end
% if ~all(pts(N,:))
%     pts(N,:) = [x(end), y(end)];
% end
%%%%%%%


% pts = [x(50001),y(50001)];
% p1 = pts(1,:);
% for i = 50002:length(x)
%     p2 = [x(i), y(i)];
%     if abs((norm(p2-p1)-l_s)) < 0.00005
% %         abs((norm(p2-p1)-l_s)) < 0.001
%         pts = [pts;p2];
%         p1 = [x(i), y(i)];
%     end
%     if length(pts) == ((N)/2)
%         break
%     end
% end
% p1 = pts(1,:);
% for i = 50000:-1:1
%     p2 = [x(i), y(i)];
%     if abs((norm(p2-p1)-l_s)) < 0.00005
%         pts = [p2;pts];
%         p1 = p2;
%     end
%     if length(pts) == N
%         break
%     end
% end
% plot(x,y,'g', 'LineWidth', 2)
% grid on
% hold on
% plot(pts(:,1),pts(:,2),'g-*','LineWidth', 5, 'MarkerSize', 12)
for i = 2:size(pts,1)
    dy = pts(i,2)-pts(i-1,2);
    dx = pts(i,1)-pts(i-1,1);
    theta(i-1,1) = atan2(dy,dx);
end
dy = pts(i-1,2)-pts(i,2);
dx = pts(i-1,1)-pts(i,1);
theta(i,1) = atan2(dy,dx);
pts = [pts theta];
xy = [x' y'];
end