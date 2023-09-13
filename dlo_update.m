

function dataN = dlo_update(N, initial, desired)

% Cable length:
L = 0.7;

% Cable points:
% N = 20;
N_cable = 40;

% Segment length
l_s = L/(N-1);

l_s_cable = L/(N_cable-1);
% Cable stiffness parameters

% K_s = 10e6/l_s;
% K_b = 1.625e1/l_s;
% K_sc = 10e6/l_s_cable;
% K_bc = 1.625e2/l_s_cable;
% K_s = 10e4;
% K_b = 2.25e2;
% K_sc = 10e4;
% K_bc = 2.25e2;
%0.7
K_s = 4.8491e3/l_s;
K_b = 0.0447/l_s;
% K_sc = 4.8491e3/l_s;
% K_bc = 0.0447/l_s;
% sochi_data_local3
% 0.5
% K_s = 16e2/l_s;
% K_b = 0.0048/l_s;

% K_sc = 16e2/l_s_cable;
% K_bc = 0.0048/l_s_cable;
% K_d = 0.5/L;
% K_s = 1e4/L;
% K_b = 1.0e3/L;
% K_d = 210/L;
v0 = zeros(N,2);
dP1 = [0.001 0.001]*0;
dP2 = [0.001 0.001]*0;

profs = des_profile(initial, N, l_s, L);
P = profs(:,1:2) ;
P(:,2) = -P(:,2);
% P = P +[0.1 0.150];
P0 = P;

profs = des_profile(initial, N_cable, l_s_cable, L);
P_c = profs(:,1:2) ;
P_c(:,2) = P_c(:,2);
P0_c = P_c;

theta1_0 = atand((P(2,2)-P(1,2))/(P(2,1)-P(1,1)));
thetan_0 = atand((P(end-1,2)-P(end,2))/(P(end-1,1)-P(end,1)));



pp = 1;
P_ = optimvar('P_', [N-4,2]);
cable_prob = optimproblem;
opts=optimoptions(@fmincon,'Display','off');
opts1=optimoptions(@fminunc,'Display','off');
sol0.P_ = P(3:N-2,:);

P_cable = optimvar('P_cable', [N_cable-4,2]);
cable_model_prob = optimproblem;
sol_cable0.P_cable = P_c(3:N_cable-2,:);

            
des_pos = des_profile(desired, N, l_s, L);
des_P = des_pos(:,1:2) +[0.100 0.000]*1;
des_P0 = des_P;

[des_pos_comp, xy] = des_profile(desired, N_cable, l_s_cable, L);
des_pos_comp = des_pos_comp(:,1:2)+[0.100 0.06000]*1;


des_shape = [];

% 
% dis  = des_P0(:,1:2) - P0(:,1:2); % get the difference between current and desired position for each point.
% max_dis = max(max(dis)); % get the maximum displacement
% n0 = ceil(max_dis/0.01); % get number of steps
% step = dis/n0; % get steps for each point
% % generate profiles
% profiles(:,:,1) = P0(:,1:2) + step;
% flag = 0;
% i = 2;
% while max(max(des_P0(:,1:2) - profiles(:,:,i-1)))>0.05
%     profiles(:,:,i) = profiles(:,:,i-1) + step;
%     A = profiles(1:end-1,:,i);
%     B = profiles(2:end,:,i);
%     if flag == 1
%         step = dis/n0;
%         flag = 0;
%     end
%     norm_d = vecnorm((B-A),2,2);
%     if max(norm_d) > 0.05
%         step = step*0.8;
%         flag = 1;
%     else
%         i = i+1;
%     end
% 
% end


profiles = inter_proiles(P, des_P, l_s);
profiles(:,:,end+1) = des_P;
% if size(profiles,3) >= 3
%     des_P = profiles(:,:,3);
%     des_P0 = profiles(:,:,3);
% end
%%%%%%%%%%%%
%%%% Control points
r_control = optimvar('r_control', 8);
prob_control = optimproblem;
sol0_control.r_control(1) = P0(1,1);
sol0_control.r_control(2) = P0(1,2);
sol0_control.r_control(3) = P0(2,1);
sol0_control.r_control(4) = P0(2,2);
sol0_control.r_control(5) = P0(end-1,1);
sol0_control.r_control(6) = P0(end-1,2);
sol0_control.r_control(7) = P0(end,1);
sol0_control.r_control(8) = P0(end,2);
%%%%%%%%%%%%%%%%%%%%
Qa = zeros(N);
lambda = 1500000;
for i = 1:N
    if i <= N/2
        Qa(i,i) = lambda;
        lambda = lambda + 5;
    else
        lambda = lambda - 5;
        Qa(i,i) = lambda;
    end
    
end
failed = false;
error_max = [];
error_avg = [];
error_std = [];
error_max1 = [];
error_avg1 = [];
error_std1 = [];

Pss = P0;
Pcss = P_c;
cnt = 2;
tic
for k =1:size(profiles,3)
    
    s = 1;
    E = 0.01;
    if k == size(profiles,3)
        E = 0.001;
    end
    while mean(vecnorm(P-profiles(:,:,k),2,2)) > 0.005

        des_P = profiles(:,:,k);

        Lp = [sqrt((des_P(1,1)-r_control(1))^2+(des_P(1,2)-r_control(2))^2);sqrt((des_P(2,1)-r_control(3))^2+(des_P(2,2)-r_control(4))^2)];
        
        for i = 3:N-2
            Lp = [Lp; sqrt((des_P(i,1)-P(i,1))^2+(des_P(i,2)-P(i,2))^2)];
        end
    
        Lp = [Lp;sqrt((des_P(end-1,1)-r_control(5))^2+(des_P(end-1,2)-r_control(6))^2);sqrt((des_P(end,1)-r_control(7))^2+(des_P(end,2)-r_control(8))^2)];
        
        %%% length constraint
        Es = (sqrt((r_control(3)-r_control(1))^2 + (r_control(4)-r_control(2))^2) - l_s).^2 ...
            + (sqrt((P(3,1)-r_control(3)).^2 + (P(3,2)-r_control(4))^2) - l_s).^2;
        for n = 4:N-2
            Es = Es + (sqrt((P(n,1)-P(n-1,1)).^2 + (P(n,2)-P(n-1,2)).^2) - l_s).^2;
        end
        Es = Es + (sqrt((r_control(5)-P(end-2,1))^2 + (r_control(6)-P(end-2,2))^2) - l_s)^2 ...
            + (sqrt((r_control(7)-r_control(5))^2 + (r_control(8)-r_control(6))^2) - l_s)^2;
        
       
        %%% Bending Energy
        
        Eb = (atan(((((P(3,2)-r_control(4))*(r_control(3)-r_control(1)))-((P(3,1)-r_control(3))*(r_control(4)-r_control(2)))))/(((P(3,1)-r_control(3))*(r_control(3)-r_control(1)))+((P(3,2)-r_control(4))*(r_control(4)-r_control(2)))))).^2 ...
            + (atan(((((P(4,2)-P(3,2))*(P(3,1)-r_control(3)))-((P(4,1)-P(3,1))*(P(3,2)-r_control(4)))))/(((P(4,1)-P(3,1))*(P(3,1)-r_control(3)))+((P(4,2)-P(3,2))*(P(3,2)-r_control(4)))))).^2;
        for n = 4:N-3
            Eb = Eb + (atan(((((P(n+1,2)-P(n,2))*(P(n,1)-P(n-1,1)))-((P(n+1,1)-P(n,1))*(P(n,2)-P(n-1,2)))))/(((P(n+1,1)-P(n,1))*(P(n,1)-P(n-1,1)))+((P(n+1,2)-P(n,2))*(P(n,2)-P(n-1,2)))))).^2;
        end
        Eb = Eb + (atan(((((r_control(6)-P(end-2,2))*(P(end-2,1)-P(end-3,1)))-((r_control(5)-P(end-2,1))*(P(end-2,2)-P(end-3,2)))))/(((r_control(5)-P(end-2,1))*(P(end-2,1)-P(end-3,1)))+((r_control(6)-P(end-2,2))*(P(end-2,2)-P(end-3,2)))))).^2 ...
            + (atan(((((r_control(8)-r_control(6))*(r_control(5)-P(end-2,1)))-((r_control(7)-r_control(5))*(r_control(6)-P(end-2,2)))))/(((r_control(7)-r_control(5))*(r_control(5)-P(end-2,1)))+((r_control(8)-r_control(6))*(r_control(6)-P(end-2,2)))))).^2;
        
        
        obj_fun = Lp'*Qa*Lp + 0.5*K_s*Es + 0.5*K_b*Eb;
        prob_control.Objective = obj_fun;
        vel_cons1 = (P(1,1)-r_control(1))^2 <= 0.0004;
        vel_cons2 = (P(1,2)-r_control(2))^2 <= 0.0004;
        vel_cons3 = (P(2,1)-r_control(3))^2 <= 0.0004;
        vel_cons4 = (P(2,2)-r_control(4))^2 <= 0.0004;
        vel_cons5 = (P(end-1,1)-r_control(5))^2 <= 0.0004;
        vel_cons6 = (P(end-1,2)-r_control(6))^2 <= 0.0004;
        vel_cons7 = (P(end,1)-r_control(7))^2 <= 0.0004;
        vel_cons8 = (P(end,2)-r_control(8))^2 <= 0.0004;

        dis1 = sqrt((l_s*cos(r_control(3)))^2 + (l_s*sin(r_control(3)))^2);  
        dis2 = sqrt((P(3,1)-r_control(1)+l_s*cos(r_control(3)))^2 + (P(3,2)-r_control(2)+l_s*sin(r_control(3)))^2);

        dis3 = sqrt((r_control(4)-l_s*cos(r_control(6))-P(N-2,1))^2 + ((r_control(5)-l_s*sin(r_control(6)))-P(N-2,2))^2);
        dis4 =  sqrt((-l_s*cos(r_control(6)))^2 + (-l_s*sin(r_control(6)))^2);
        dis_cons1 = dis1  <= l_s + 0.002;
        dis_cons2 = dis2  <= l_s + 0.002;
        dis_cons3 = dis3  <= l_s + 0.002;
        dis_cons4 = dis4  <= l_s + 0.002;

        prob_control.Constraints.cons1 = vel_cons1;
        prob_control.Constraints.cons2 = vel_cons2;
        prob_control.Constraints.cons3 = vel_cons3;
        prob_control.Constraints.cons4 = vel_cons4;
        prob_control.Constraints.cons5 = vel_cons5;
        prob_control.Constraints.cons6 = vel_cons6;
        prob_control.Constraints.cons7 = vel_cons7;
        prob_control.Constraints.cons8 = vel_cons8;
%         
%         
        sol_control = solve(prob_control, sol0_control, 'Options', opts);
        
        
        %%%%%%%%%%%%%%
        
        theta1 = atan2((sol_control.r_control(4)-sol_control.r_control(2)),(sol_control.r_control(3)-sol_control.r_control(1)));
        thetan = atan2((sol_control.r_control(8)-sol_control.r_control(6)),(sol_control.r_control(7)-sol_control.r_control(5)));
        
        
        dP1 = [sol_control.r_control(1) sol_control.r_control(2)] - P(1,:);
        dP2 = [sol_control.r_control(3) sol_control.r_control(4)] - P(2,:);
        dP3 = [sol_control.r_control(5) sol_control.r_control(6)] - P(end-1,:);
        dP4 = [sol_control.r_control(7) sol_control.r_control(8)] - P(end,:);
        
        P(1,:) = P(1,:)+dP1;
      	P(2,:) = P(2,:)+dP2;
        P(end-1,:) = P(end-1,:)+dP3;
        P(end,:) = P(end,:)+dP4;
        
         %%% length constraint
        Es = (sqrt((P(2,1)-P(1,1))^2 + (P(2,2)-P(1,2))^2) - l_s).^2 ...
            + (sqrt((P_(1,1)-P(2,1)).^2 + (P_(1,2)-P(2,2))^2) - l_s).^2;
        for n = 2:N-4
            Es = Es + (sqrt((P_(n,1)-P_(n-1,1)).^2 + (P_(n,2)-P_(n-1,2)).^2) - l_s).^2;
        end
        Es = Es + (sqrt((P(end-1,1)-P_(end,1))^2 + (P(end-1,2)-P_(end,2))^2) - l_s)^2 ...
            + (sqrt((P(end,1)-P(end-1,1))^2 + (P(end,2)-P(end-1,2))^2) - l_s)^2;
        
        %%% Bending Energy
        
        Eb = (atan(((((P_(1,2)-P(2,2))*(P(2,1)-P(1,1)))-((P_(1,1)-P(2,1))*(P(2,2)-P(1,2)))))/(((P_(1,1)-P(2,1))*(P(2,1)-P(1,1)))+((P_(1,2)-P(2,2))*(P(2,2)-P(1,2)))))).^2 ...
            + (atan(((((P_(2,2)-P_(1,2))*(P_(1,1)-P(2,1)))-((P_(2,1)-P_(1,1))*(P_(1,2)-P(2,2)))))/(((P_(2,1)-P_(1,1))*(P_(1,1)-P(2,1)))+((P_(2,2)-P_(1,2))*(P_(1,2)-P(2,2)))))).^2;
        for n = 2:N-5
            Eb = Eb + (atan(((((P_(n+1,2)-P_(n,2))*(P_(n,1)-P_(n-1,1)))-((P_(n+1,1)-P_(n,1))*(P_(n,2)-P_(n-1,2)))))/(((P_(n+1,1)-P_(n,1))*(P_(n,1)-P_(n-1,1)))+((P_(n+1,2)-P_(n,2))*(P_(n,2)-P_(n-1,2)))))).^2;
        end
        Eb = Eb + (atan(((((P(end-1,2)-P_(end,2))*(P_(end,1)-P_(end-1,1)))-((P(end-1,1)-P_(end,1))*(P_(end,2)-P_(end-1,2)))))/(((P(end-1,1)-P_(end,1))*(P_(end,1)-P_(end-1,1)))+((P(end-1,2)-P_(end,2))*(P_(end,2)-P_(end-1,2)))))).^2 ...
            + (atan(((((P(end,2)-P(end-1,2))*(P(end-1,1)-P_(end,1)))-((P(end,1)-P(end-1,1))*(P(end-1,2)-P_(end,2)))))/(((P(end,1)-P(end-1,1))*(P(end-1,1)-P_(end,1)))+((P(end,2)-P(end-1,2))*(P(end-1,2)-P_(end,2)))))).^2;

        cable_prob.Objective = 0.5*K_s*Es + 0.5*K_b*Eb;

        [sol, ee] = solve(cable_prob, sol0,  'Options', opts1);
        
        P(3:N-2,:) = sol.P_;
        
        sol0.P_ = P(3:N-2,:);
        
   
%         dP1 = [sol_control.r_control(1) sol_control.r_control(2)] - P_c(1,:);
%         dP2 = [sol_control.r_control(1)+l_s_cable*cos(theta1) sol_control.r_control(2)+l_s_cable*sin(theta1)] - P_c(2,:);
%         dP3 =  [sol_control.r_control(7)-l_s_cable*cos(thetan) sol_control.r_control(8)-l_s_cable*sin(thetan)] - P_c(end-1,:);
%         dP4 = [sol_control.r_control(7) sol_control.r_control(8)] - P_c(end,:);
%        
%         
%         P_c(1,:) = P(1,:);
%       	P_c(2,:) = P_c(2,:)+dP2;
%         P_c(end-1,:) = P_c(end-1,:)+dP3;
%         P_c(end,:) = P(end,:);
%         
%       %%% length constraint
%         Es = (sqrt((P_c(2,1)-P_c(1,1))^2 + (P_c(2,2)-P_c(1,2))^2) - l_s_cable).^2 ...
%             + (sqrt((P_cable(1,1)-P_c(2,1)).^2 + (P_cable(1,2)-P_c(2,2))^2) - l_s_cable).^2;
%         for n = 2:N_cable-4
%             Es = Es + (sqrt((P_cable(n,1)-P_cable(n-1,1)).^2 + (P_cable(n,2)-P_cable(n-1,2)).^2) - l_s_cable).^2;
%         end
%         Es = Es + (sqrt((P_c(end-1,1)-P_cable(end,1))^2 + (P_c(end-1,2)-P_cable(end,2))^2) - l_s_cable)^2 ...
%             + (sqrt((P_c(end,1)-P_c(end-1,1))^2 + (P_c(end,2)-P_c(end-1,2))^2) - l_s_cable)^2;
%         
%         %%% Bending Energy
%         
%         Eb = (atan(((((P_cable(1,2)-P_c(2,2))*(P_c(2,1)-P_c(1,1)))-((P_cable(1,1)-P_c(2,1))*(P_c(2,2)-P_c(1,2)))))/(((P_cable(1,1)-P_c(2,1))*(P_c(2,1)-P_c(1,1)))+((P_cable(1,2)-P_c(2,2))*(P_c(2,2)-P_c(1,2)))))).^2 ...
%             + (atan(((((P_cable(2,2)-P_cable(1,2))*(P_cable(1,1)-P_c(2,1)))-((P_cable(2,1)-P_cable(1,1))*(P_cable(1,2)-P_c(2,2)))))/(((P_cable(2,1)-P_cable(1,1))*(P_cable(1,1)-P_c(2,1)))+((P_cable(2,2)-P_cable(1,2))*(P_cable(1,2)-P_c(2,2)))))).^2;
%         for n = 2:N_cable-5
%             Eb = Eb + (atan(((((P_cable(n+1,2)-P_cable(n,2))*(P_cable(n,1)-P_cable(n-1,1)))-((P_cable(n+1,1)-P_cable(n,1))*(P_cable(n,2)-P_cable(n-1,2)))))/(((P_cable(n+1,1)-P_cable(n,1))*(P_cable(n,1)-P_cable(n-1,1)))+((P_cable(n+1,2)-P_cable(n,2))*(P_cable(n,2)-P_cable(n-1,2)))))).^2;
%         end
%         Eb = Eb + (atan(((((P_c(end-1,2)-P_cable(end,2))*(P_cable(end,1)-P_cable(end-1,1)))-((P_c(end-1,1)-P_cable(end,1))*(P_cable(end,2)-P_cable(end-1,2)))))/(((P_c(end-1,1)-P_cable(end,1))*(P_cable(end,1)-P_cable(end-1,1)))+((P_c(end-1,2)-P_cable(end,2))*(P_cable(end,2)-P_cable(end-1,2)))))).^2 ...
%             + (atan(((((P_c(end,2)-P_c(end-1,2))*(P_c(end-1,1)-P_cable(end,1)))-((P_c(end,1)-P_c(end-1,1))*(P_c(end-1,2)-P_cable(end,2)))))/(((P_c(end,1)-P_c(end-1,1))*(P_c(end-1,1)-P_cable(end,1)))+((P_c(end,2)-P_c(end-1,2))*(P_c(end-1,2)-P_cable(end,2)))))).^2;
% 
%         
%         cable_model_prob.Objective = 0.5*K_sc*Es + 0.5*K_bc*Eb;
% 
%         [sol_cable, ee] = solve(cable_model_prob, sol_cable0,  'Options', opts1);
%         
%         P_c(3:N_cable-2,:) = sol_cable.P_cable;
%         
%         sol_cable0.P_cable = P_c(3:N_cable-2,:);


        theta1_0 = atan2((P(2,2)-P(1,2)),(P(2,1)-P(1,1)));
        thetan_0 = atan2((P(end-1,2)-P(end,2)),(P(end-1,1)-P(end,1)));
        sol0_control.r_control(1) = sol_control.r_control(1);
        sol0_control.r_control(2) = sol_control.r_control(2);
        sol0_control.r_control(3) = sol_control.r_control(3);
        sol0_control.r_control(4) = sol_control.r_control(4);
        sol0_control.r_control(5) = sol_control.r_control(5);
        sol0_control.r_control(6) = sol_control.r_control(6);
        sol0_control.r_control(7) = sol_control.r_control(7);
        sol0_control.r_control(8) = sol_control.r_control(8);
        

%         for j= 1:size(profiles,3)
% 
%             plot(profiles(:,1,j), profiles(:,2,j), 'b-', 'LineWidth', 1)
%             hold on
%         end
% 
%         plot(des_P0(:,1), des_P0(:,2), 'g*-', 'LineWidth', 3)% hold on
%         hold on
%         plot(P(:,1), P(:,2), 'r*-','linewidth',3)
%         plot(P_c(:,1), P_c(:,2), 'ko-','linewidth',3)
% %         axis([0 1.2 -0.4 1])
%         
%         drawnow
%         hold off
%         error1 = vecnorm(P_c-des_pos_comp,2,2);
        error = vecnorm(P-des_P0,2,2);
        error_max = [error_max max(error)];
        
        error_avg = [error_avg mean(error)];
        error_std = [error_std std(error)];
%         error_max1 = [error_max1 max(error1)];
%         
%         error_avg1 = [error_avg1 mean(error1)];
%         error_std1 = [error_std1 std(error1)];
        
        Pss(:,:,cnt) = P;
        Pcss(:,:,cnt) = P_c;
        cnt = cnt + 1;
%         pause(0.1)
%             failed = true;
%             disp('failed')
            L_f = cumsum(vecnorm(P(2:end,:)-P(1:end-1,:),2,2));
            L_f = L_f(end);
            L_1 = cumsum(vecnorm(des_P0(2:end,:)-des_P0(1:end-1,:),2,2));
            L_1 = L_1(end);
% L_init = L_init(end)

%             L_f = L_f(end)
%             break
%         end
%         k = k+1;
        s = s+1;
        if s == 20
            break;
        end

    end
%     if failed
%         break 
%     end
%     ss = vecnorm(P(1:end-1,:)-P(2:end,:),2,2)
%     aa = cumsum(ss)
%     aa(end)
end
time = toc;

Pss(:,:,cnt) = P;
% Pcss(:,:,cnt) = P_c;
dataN.Pss = Pss;
% dataN.Pcss = Pcss;
dataN.profiles = size(profiles,3);
dataN.error_max = error_max;
dataN.error_mean = error_avg;
% error_std
dataN.error_std = error_std;
% des_pos_comp
% dataN.error_max1 = error_max1;
% dataN.error_mean1 = error_avg1;
% error_std
% dataN.error_std1 = error_std1;
dataN.desP = des_P0;
% dataN.desPxy = xy;
dataN.time = time;
% figure
% plot(des_pos_comp(:,1), des_pos_comp(:,2), 'g-', 'LineWidth', 3)% hold on
% hold on
% plot(P_c(:,1), P_c(:,2), 'r','linewidth',2)
% figure
% plot(error_max, 'k-', 'LineWidth', 3)
% hold on
% plot(error_avg, 'k-', 'LineWidth', 3)
end
    
