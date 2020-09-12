function task3(Number_of_APs,AP,trajectories,Q,inv_R,F,points_x,points_y,AP_IDs)
    ID_trajectory = 5; % choose your trajectory here
    trajectory = cell2mat(trajectories(1,ID_trajectory)); 
    rhoTraining = importdata("GR35/Task3_rhoUEAP_GR35.mat"); % TOA measurements of trajectories
    measurement_of_current_trajectory = cell2mat(rhoTraining(1,ID_trajectory)); %taking TOA of this trajectory
    
    % initialization
    x_initialized = trajectory(1,:); % we read the first line (measurement)

    P_t_minus_t_minus = eye(4); 

    x_hat_t_minus_t_minus = x_initialized'; % transposing x_initilized
    predictions= zeros(4,200); % used for plotting

    % defining C for post accuracy evaluation
    C_stored = cell(1,200);
    sigma_h = zeros(1,200);
    CEP95 = zeros(1,200);

    for i = 1:(size(trajectory,1)) % 200 entries for a single trajectory
        % PREDICTION SECTION OF KALMAN FILTER
        % assuming F, Q constant

        if (i == 1)
            x_hat_t_t_minus = x_hat_t_minus_t_minus;
            P_t_t_minus = P_t_minus_t_minus;    
        else
            x_hat_t_t_minus = F * x_hat_t_minus_t_minus; 
            P_t_t_minus = F * P_t_minus_t_minus * F' + Q; 
        end

        predictions(:,i)=x_hat_t_t_minus(:,1); % for plotting
        
        %UPDATE SECTION OF KALMAN FILTER
        % CALCULATING H
        H = zeros(Number_of_APs,4);
        h = zeros(Number_of_APs,1);
        for j = 1 : Number_of_APs
            h(j,1) = sqrt((AP(j,1) - x_hat_t_t_minus(1,1))^2 + (AP(j,2) - x_hat_t_t_minus(2,1))^2);

            syms ux uy; % computing H through derivatives of h
            f = sqrt((AP(j,1) - ux)^2 + (AP(j,2) - uy)^2);
            a_dx = diff(f, ux);
            a_dy = diff(f, uy);
            a_ux_final = vpa(subs(a_dx, [ux, uy], [trajectory(i,1), trajectory(i,2)]), 5);
            a_uy_final = vpa(subs(a_dy, [ux, uy], [trajectory(i,1), trajectory(i,2)]), 5);

            H(j,1) = a_ux_final; 
            H(j,2) = a_uy_final;
            H(j,3) = 0; % h is constant w.r.t velocity for TOA measurements => derivative w.r.t velocity = 0
            H(j,4) = 0; % same here
        end

        inv_P = inv(P_t_t_minus);
        tmp = (H'* inv_R * H); % optimizing inverse operation. that's why we use \ instead of inv()
        tempA = inv(inv_P + tmp);

        % updating variables
        x_hat_t_t = x_hat_t_t_minus + tempA * H' * inv_R * (measurement_of_current_trajectory(i, :)' - h);
        P_t_t = inv(inv_P + tmp);

        x_hat_t_minus_t_minus = x_hat_t_t;
        P_t_minus_t_minus = P_t_t;

        %computing performance metrics
        %CALCULATING C
        C = inv(H(:,1:2)'*inv_R*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies of the various APs are different
        C_stored(1, i) = mat2cell(C,2); % storing C
        sigma_h(1, i) = sqrt(C(1,1) + C(2,2)); % drms
        CEP95(1, i) = 2 * sigma_h(1, i); % CEP
    end

 %% PLOT Task 3 KALMAN TRACKING
    figure
    plot_dir(predictions(1,:)', predictions(2, :)', i,'g--o','p');
    hold on
    scatter(points_x, points_y,100, 'b','^'); %APs
    text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
    hold on
    plot_dir(trajectory(:,1), trajectory(:,2), i,'r-*','b');
    legend('Predicted UE location','Predicted movement','AP','Real UE location','Real movement','Location','bestoutside')
    hold on
    xlabel('Ux')
    ylabel('Uy')
    title('Task 3:Predicted UE Position Ux,Uy')
    grid on
    grid minor

    %we assume sampling rate=Ts=1s
    timeVector  = linspace(1,size(trajectory,1),size(trajectory,1));

    %plot vx
    figure
    plot(timeVector, predictions(3,:),'c-o');
    hold on
    legend('Predicted Vx')
    hold on
    xlabel('Time')
    ylabel('Vx')
    title('Task 3:Predicted UE Velocity Vx ')
    grid on
    grid minor

    %plot vy
    figure
    plot(timeVector,  predictions(4,:),'g-o');
    hold on
    legend('Predicted Vy')
    hold on
    xlabel('Time')
    ylabel('Vy')
    title('Task 3:Predicted UE Velocity Vy ')
    grid on
    grid minor

    %plot both vx,vy but this time w.r.t. time
    figure
    plot(timeVector,  predictions(3,:),'c-o');
    hold on
    plot(timeVector,  predictions(4,:),'g-o');
    hold on
    legend('Predicted velocity Vx','Predicted velocity Vy')
    hold on
    xlabel('Time')
    ylabel('Velocity')
    title('Task 3:comparison amid Predicted UE Velocity Vx and Vy ')
    grid on
    grid minor
end
