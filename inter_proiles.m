function inter_prof = inter_proiles(P, des, l_s)

    N = size(P,1);
%     des(:,2)-P(:,2)
    [max_dis, arg_idx] = max(abs(des(:,2)-P(:,2)));
%     arg_idx
    theta_p = zeros(1,N);
    theta_d = zeros(1,N);
    thetaD = zeros(1,N);
    thetaD_old = zeros(1,N);


    for i =arg_idx+1:size(P,1)
        theta_p(i) = atan2((P(i,2)-P(i-1,2)),(P(i,1)-P(i-1,1)));
        theta_d(i) = atan2((des(i,2)-des(i-1,2)),(des(i,1)-des(i-1,1)));
        thetaD(i) = theta_d(i) - theta_p(i);
        thetaD_old(i) = thetaD(i);
        thetaD(i) = mod(theta_d(i) - theta_p(i),2*pi); 
        if  thetaD(i) > pi
            thetaD(i) = -(2*pi-thetaD(i));
        end
        
%         if thetaD(i) > 2*pi
%             thetaD(i) = thetaD(i) - 2*pi;
%         elseif thetaD(i) < -2*pi
%             thetaD(i) = thetaD(i) + 2*pi;
%          end
    end

    for i = arg_idx-1:-1:1

        theta_p(i) = atan2((P(i,2)-P(i+1,2)),(P(i,1)-P(i+1,1)));
        theta_d(i) = atan2((des(i,2)-des(i+1,2)),(des(i,1)-des(i+1,1)));
        thetaD(i) =  theta_d(i) - theta_p(i);
        thetaD_old(i) = thetaD(i);
        thetaD(i) = mod(theta_d(i) - theta_p(i),2*pi); 
        if  thetaD(i) > pi
            thetaD(i) = -(2*pi-thetaD(i));
        end
%         if thetaD(i) > pi

%             thetaD(i) = thetaD(i) - 2*pi;
%         elseif thetaD(i) < -pi
%             thetaD(i) = thetaD(i) + 2*pi;
%          end
    end
%     rad2deg(theta_p)
    rad2deg(thetaD_old);
    rad2deg(thetaD);
%     thetaD;
    dis = des(arg_idx,2) - P(arg_idx,2);
    rad2deg(theta_p);
    rad2deg(theta_d);
    dis1 = des(arg_idx,:) - P(arg_idx,:);

%     axis([0 0.8 0 0.8])
    

    n0 = abs(floor(dis/(0.2*l_s)));
    step_l = dis1/n0;
    step_a = thetaD/n0;
    rad2deg(step_a);

    inter_prof = zeros(N,2,n0);

    for i = 1:n0
        if i == 1
            inter_prof(arg_idx,:,i) = P(arg_idx,:) + step_l;
        else

            inter_prof(arg_idx,:,i) = inter_prof(arg_idx,:,i-1) + step_l;
        end
        for j = arg_idx+1:N
            inter_prof(j,:,i) = inter_prof(j-1,:,i) +[l_s*cos(theta_p(j)+i*step_a(j)) l_s*sin(theta_p(j)+i*step_a(j))];
        end
        for j = arg_idx-1:-1:1

            inter_prof(j,:,i) = inter_prof(j+1,:,i) +[l_s*cos(theta_p(j)+i*step_a(j)) l_s*sin(theta_p(j)+i*step_a(j))];
        end
    end





end
