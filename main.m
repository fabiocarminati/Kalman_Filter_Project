clear all, clc, close all
set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

addpath Functions

%% scenario settings (4x4 Km)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;
task1aSwitch=false;
task2Switch=false;
task3Switch=false;
task3TestingSwitch=false;
task4Switch=false;
task5Switch=false;
task6Switch=true;
%% importing data
AP_Data = importdata("GR35/Task1a_rhoUEAP.mat");
Number_of_APs = size(AP_Data, 1); % number of APs is mapped in row length
u = [0, 0]; % reference point

%% position of the Access Points
%parameters.numberOfAP = Number_of_APs;
%[ AP ] = GeneratepositionofAP(parameters);
points_x = [];
points_y = [];
AP = zeros(Number_of_APs, 2);

for i = 1:Number_of_APs
    toa_measurement = AP_Data(i,1);
    aoa_measurement = AP_Data(i,2);
    sign=1;
    if(aoa_measurement>pi/2) %if TOA is not in the range [-pi/2,pi/2] the atan matlab function doesnt work. So we
                               %always want values in that range; if
                               %necessary we chenge aoa_measurement value
                               %and consequantly add a - (minus sign) on
                               %the result of the system
                               
         aoa_measurement=aoa_measurement-pi;
         sign=-1;
     else
         if(aoa_measurement<-pi/2)
             aoa_measurement=aoa_measurement+pi;
             sign=-1;
         end
     end
    
    syms Sx Sy 
    eq_toa = sqrt((Sx -u(1))^2 + (Sy - u(2))^2) == toa_measurement;
    a = (u(2) - Sy) / (u(1) - Sx);
    eq_aoa =atan(a) == aoa_measurement;    
    solution = vpasolve([eq_toa, eq_aoa], [Sx, Sy]);
    %solution = solve([eq_toa, eq_aoa], [Sx, Sy]); %it works but with
    %warning (It wants to use vpasolve)
    
    SxSol = round(solution.Sx);
    points_x = [points_x sign*SxSol];
 
    SySol = round(solution.Sy);
    points_y = [points_y sign*SySol];
    
    AP(i, 1) = sign*SxSol; 
    AP(i, 2) = sign*SySol;
       
end
AP_IDs={'AP1','AP2','AP3','AP4','AP5','AP6','AP7','AP8'}

%% Plot APs
switch task1aSwitch
    case true
        figure
        scatter(points_x, points_y,100, 'b','^'); %APs
        text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
        hold on
        plot(0,0, 'r-x');
        xlabel('APx')
        ylabel('APy')
        title('AP exact positions')
        grid on
        grid minor
    case false
        fprintf('no plot of APs \n');
end
        

%% Task_1b
u = [500 -800];

rho = importdata("GR35/Task1b_rhoUEAP.mat");
sigma_squared = zeros(1,8);

for i = 1:Number_of_APs
    temp = [];
    
    distance_u_ap = sqrt((AP(i,1) - u(1))^2 + (AP(i,2) - u(2))^2);
    
    for j = 1:size(rho, 1) % 5000
        n = rho(j, i) - distance_u_ap;
        temp = [temp n]; % We don't know how to perform the append operation
    end
    
    sigma_squared(i) = cov(temp);

end

R = diag(sigma_squared);
R = round(R); % covariance matrix

%% Task 2 motion models
trajectories = importdata("GR35/Task2_trajectory_GR35.mat");
Ts = 1; %second
TotalSimulationTime = 20; %s
wx = [0];
wy = [0];
for i = 1:100
    measures = cell2mat(trajectories(1,i));
    
    for temp = 2:size(measures,1)
        vx_t_minus_one = measures(temp - 1, 3);
        vx_t = measures(temp, 3);
        
        vy_t_minus_one = measures(temp - 1, 4);
        vy_t = measures(temp, 4);
        %vy_t = measures(1,4);
        
        wx = [wx ((vx_t - vx_t_minus_one) / Ts)]; % 1x19901 because every trajectory has 200 points and we start the inner for from 2 so (100x200(points) - 100)
        wy = [wy ((vy_t - vy_t_minus_one) / Ts)];
    end
end
mean_x = mean(wx);
mean_y = mean(wy);

var_x = var(wx);
var_y = var(wy);

fprintf('mean of all wx: %d\n', mean_x);
fprintf('variance of all wx: %d\n', var_x);
fprintf('mean of all wy: %d\n', mean_y);
fprintf('variance of all wy: %d\n', var_y);

w_a = zeros(2, size(wx,2));
for i = 1:size(wx,2)
    w_a(1,i) = wx(i);
    w_a(2,i) = wy(i);
end

switch task2Switch
    case true 
        F = [eye(2)     , Ts*eye(2);
             zeros(2,2) ,    eye(2)];
        L = [0.5*Ts^2*eye(2);Ts*eye(2)];

        v = zeros(2,TotalSimulationTime);
        v(:,1) = [7;2];
        x = [zeros(2,TotalSimulationTime);v];

        for t = 2:TotalSimulationTime
            %w_a = sigma_a*randn(2,1);
            x(:,t) = F*x(:,t-1) + L*w_a;
        end


        %% plot task 2 trajectories

        figure
        for i = 1:size(trajectories,2)/2
            trajectory = cell2mat(trajectories(i));
            plot(trajectory(:,1), trajectory(:,2));
            hold on;
        end
        xlabel('Valori X')
        ylabel('Valori Y')
        title('User trajectories')
        
        %% plot task 2 W_at

        figure
        delta_v_x = zeros(1, (size(w_a,2)-1)/100);
        delta_v_y = zeros(1, (size(w_a,2)-1)/100);

        for i = 2:size(w_a, 2)/100

        delta_v_x(i-1) = w_a(1,i - 1) - w_a(1,i);
        delta_v_y(i-1) = w_a(2,i - 1) - w_a(2,i);
        end

        plot(1:size(delta_v_x, 2), delta_v_x);

        hold on
        plot(1:size(delta_v_y, 2), delta_v_y);
        hold on

        plot(1:size(delta_v_x, 2), repelem(mean_x, size(delta_v_x, 2)));
        hold on
        % TODO, verify this if possible plot v_x and v_y combined
        plot(1:size(delta_v_y, 2), repelem(mean_y, size(delta_v_y, 2)));
        xlabel('t')
        ylabel('m/(s\^2)')
        title('random acceleration')
case false
        fprintf('no plot task 2 \n');
end

%% TASK 3 (ONLY FIRST TRAJECTORY)
ID_trajectory = 5; % choose your trajectory here
trajectory = cell2mat(trajectories(1,ID_trajectory)); % if you change the trajectory, remember to change also the measurement_of_current_trajectory
rhoTraining = importdata("GR35/Task3_rhoUEAP_GR35.mat"); % TOA measurements of trajectories
measurement_of_current_trajectory = cell2mat(rhoTraining(1,ID_trajectory)); %taking TOA of this trajectory
Ts = 1; % in seconds. time sampling

F = [eye(2)     , Ts*eye(2);
     zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2);Ts*eye(2)];

% initialization
x_initialized = trajectory(1,:); % we read the first line (measurement)

%x_var = cov(trajectory); % finding covariance of x0 (entire trajectory ux,uy,vx,vy)
sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)

P_t_minus_t_minus = eye(4);% cov(trajectory); % initializing P to cov(x0) since we don't know the exact starting point

Q = sigma_squared_a * (L * L'); % finding Q

x_hat_t_minus_t_minus = x_initialized'; % transposing x_initilized
predictions= zeros(4,200); % used for plotting

inv_R = inv(R); % calculating only one time, the inverse of R

% defining C
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
        x_hat_t_t_minus = F * x_hat_t_minus_t_minus; % remember to update at the end x_hat_t_minus_t_minus
        P_t_t_minus = F * P_t_minus_t_minus * F' + Q; % remember to update at the end P_t_minus_t_minus
    end
    
    predictions(:,i)=x_hat_t_t_minus(:,1); % for plotting
    %UPDATE SECTION OF KALMAN FILTER
    
    % CALCULATING H
    H = zeros(Number_of_APs,4);
    h = zeros(Number_of_APs,1);
    for j = 1 : Number_of_APs
        
        h(j,1) = sqrt((AP(j,1) - x_hat_t_t_minus(1,1))^2 + (AP(j,2) - x_hat_t_t_minus(2,1))^2);
        
        syms ux uy; % computing H through derivatives
        f = sqrt((AP(j,1) - ux)^2 + (AP(j,2) - uy)^2);
        a_dx = diff(f, ux);
        a_dy = diff(f, uy);
        a_ux_final = vpa(subs(a_dx, [ux, uy], [trajectory(i,1), trajectory(i,2)]), 5);
        a_uy_final = vpa(subs(a_dy, [ux, uy], [trajectory(i,1), trajectory(i,2)]), 5);
        
        H(j,1) = a_ux_final; % no minus because we use the approach of the derivatives.
        H(j,2) = a_uy_final;
        H(j,3) = 0; % h is constant w.r.t velocity => derivative w.r.t velocity = 0
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
    C = inv(H(:,1:2)'*inv_R*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies are different among themselves;
    C_stored(1, i) = mat2cell(C,2); % storing C
    sigma_h(1, i) = sqrt(C(1,1) + C(2,2)); % drms
    CEP95(1, i) = 2 * sigma_h(1, i); % CEP
    
    %theta = 0.5 * atan((2*C(1,2))/(C(1,1)^2 - C(2,2)^2));
    switch task3Switch
        case true
            fprintf('doing Task 3 \n');
        case false
            fprintf('stopping kalman iteration Task 3 \n');
            break;
    end
     
       
end

%% TEST PER CAPIRE SE NLS funziona->

%CASO A: trajectory(20,:)per stimare u_hat(ux,uy); punto
%esatto=(822.8452  226.7787)-->anche partendo da (0,0) arrivo dopo 10
%iterazioni molto vicino a (819.1280  228.9507)

%CASO b: trajectory(35,:)per stimare u_hat(ux,uy); punto
%esatto=(897.0610  261.7183)-->anche partendo da (0,0) arrivo dopo 10
%iterazioni molto vicino a ( 897.9944 257.2620)

%u_hat(1,1)=trajectory(20,1); %initialization of u_hat coordinates
%u_hat(2,1)=trajectory(20,2);
switch task3TestingSwitch
    case true
        u_hat=zeros(2,1); %
        for w=1:4 %max_iterations =4 for now
            delta_ro_k=zeros(Number_of_APs,1);
            H_k=zeros(8,2);
            for j = 1 : Number_of_APs   
                %COMPUTATION          
                syms ux_hat uy_hat; % computing H through derivatives
                f = sqrt((AP(j,1) - ux_hat)^2 + (AP(j,2) - uy_hat)^2);
                a_dx = diff(f, ux_hat);
                a_dy = diff(f, uy_hat);
                H_k(j,1) = vpa(subs(a_dx, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                H_k(j,2) = vpa(subs(a_dy, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                delta_ro_j_k=measurement_of_current_trajectory(10,j)-sqrt((AP(j,1) - u_hat(1,1))^2 + (AP(j,2) - u_hat(2,1))^2);% (i,j) dovra essere!
                delta_ro_k(j,1)=delta_ro_j_k;
            end
            %INVERSION
            delta_u_k=inv(H_k'*H_k)*H_k'*delta_ro_k;
            %delta_u_k=inv(H_k'*inv_R*H_k)*H_k'*inv_R*delta_ro_k; %WNLS
            %UPDATE SOLUTION
            u_hat=u_hat+delta_u_k;
        end
        %EKF code
    case false
        fprintf('do not test uhat\n');
end
%%

 %% PLOT Task 3 KALMAN TRACKING
switch task3Switch
    case true
        figure
        plot(predictions(1,:), predictions(2,:),'g--o');
        hold on
        scatter(points_x, points_y,100, 'b','^'); %APs
        text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
        hold on
        plot(trajectory(:,1), trajectory(:,2),'r-*');

        legend('Predicted UE location','AP','Real UE location','Location','bestoutside')
        hold on
        xlabel('Ux')
        ylabel('Uy')
        title('Task 3:Predicted UE Position Ux,Uy')
        grid on
        grid minor

        %we assume sampling rate=Ts=1
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
    case false
        fprintf('no plot kalman iteration Task 3 \n');
end   

%% Task 4
switch task4Switch
    case true 
        % we suppose that the motion model is still M3
        rhoTraining4 = importdata("GR35/Task4_rhoUEAP_GR35.mat"); % TOA measurements of trojectories
        % INITIALIZATION 
        Ts = 1; % in seconds. time sampling

        F = [eye(2)     , Ts*eye(2);
             zeros(2,2) ,    eye(2)];
        L = [0.5*Ts^2*eye(2);Ts*eye(2)];

        sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)
        P_t_minus_t_minus = eye(4);% cov(trajectory); % initializing P

        Q = sigma_squared_a * (L * L'); % finding Q

        x_hat_t_minus_t_minus = zeros(4,1); 
        predictions= zeros(4, size(rhoTraining4,1)); % used for plotting

        vx_zero = 50 *1000 / 3600; % in m/s

        u_hat_trajectory = zeros(size(rhoTraining4,1), 2);

        inv_R = inv(R); % calculating only one time, the inverse of R

        % defining C
        C_stored_4 = cell(1,size(rhoTraining4,1));
        sigma_h_4 = zeros(1,size(rhoTraining4,1));
        CEP95_4 = zeros(1,size(rhoTraining4,1));
        % END initialization kalman filter task 4

        u_hat=zeros(2,1); % initialization of u_hat coordinates
        % not the same H of kalman filter(which is a 8*4)
        for i = 1:(size(rhoTraining4,1))
            for w=1:4 % max_iterations
                delta_ro_k=zeros(Number_of_APs,1);
                H_k=zeros(8,2);
                for j = 1 : Number_of_APs   
                    % COMPUTATION          
                    syms ux_hat uy_hat; % computing H through derivatives
                    f = sqrt((AP(j,1) - ux_hat)^2 + (AP(j,2) - uy_hat)^2);
                    a_dx = diff(f, ux_hat);
                    a_dy = diff(f, uy_hat);
                    H_k(j,1) = vpa(subs(a_dx, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                    H_k(j,2) = vpa(subs(a_dy, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                    delta_ro_j_k=rhoTraining4(i,j)-sqrt((AP(j,1) - u_hat(1,1))^2 + (AP(j,2) - u_hat(2,1))^2);% (i,j) dovra essere!
                    delta_ro_k(j,1)=delta_ro_j_k;
                end
                %INVERSION
                delta_u_k=inv(H_k'*H_k)*H_k'*delta_ro_k;
                %UPDATE SOLUTION
                u_hat=u_hat+delta_u_k;

                %saving u_hat for plotting
                u_hat_trajectory(i,:) = u_hat';
            end

            %EKF code using x_hat as current location 

            % PREDICTION SECTION OF KALMAN FILTER
            % assuming F, Q constant

            if (i == 1)
                x_hat_t_minus_t_minus = u_hat; % still missing the velocity part

                x_hat_t_minus_t_minus(3,1) = vx_zero;
                x_hat_t_minus_t_minus(4,1) = 0; % we suppose velocity on y axis = 0

                x_hat_t_t_minus = x_hat_t_minus_t_minus;
                P_t_t_minus = P_t_minus_t_minus;    
            else
                x_hat_t_t_minus = F * x_hat_t_minus_t_minus; % remember to update at the end x_hat_t_minus_t_minus
                P_t_t_minus = F * P_t_minus_t_minus * F' + Q; % remember to update at the end P_t_minus_t_minus
            end

            predictions(:,i)=x_hat_t_t_minus(:,1); % for plotting
            %UPDATE SECTION OF KALMAN FILTER

            % CALCULATING H
            H = zeros(Number_of_APs,4);
            h = zeros(Number_of_APs,1);
            for j = 1 : Number_of_APs

                h(j,1) = sqrt((AP(j,1) - x_hat_t_t_minus(1,1))^2 + (AP(j,2) - x_hat_t_t_minus(2,1))^2);

                syms ux uy; % computing H through derivatives
                f = sqrt((AP(j,1) - ux)^2 + (AP(j,2) - uy)^2);
                a_dx = diff(f, ux);
                a_dy = diff(f, uy);
                a_ux_final = vpa(subs(a_dx, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);
                a_uy_final = vpa(subs(a_dy, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);

                H(j,1) = a_ux_final; % no minus because we use the approach of the derivatives.
                H(j,2) = a_uy_final;
                H(j,3) = 0; % h is constant w.r.t velocity => derivative w.r.t velocity = 0
                H(j,4) = 0; % same here
            end

            inv_P = inv(P_t_t_minus);

            tmp = (H'* inv_R * H); % optimizing inverse operation. that's why we use \ instead of inv()
            tempA = inv(inv_P + tmp);

            % updating variables
            x_hat_t_t = x_hat_t_t_minus + tempA * H' * inv_R * (rhoTraining4(i, :)' - h);
            P_t_t = inv(inv_P + tmp);

            x_hat_t_minus_t_minus = x_hat_t_t;
            P_t_minus_t_minus = P_t_t;

            %computing performance metrics
            %CALCULATING C
            C = inv(H(:,1:2)'*inv_R*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies are different among themselves;
            C_stored_4(1, i) = mat2cell(C,2); % storing C
            sigma_h_4(1, i) = sqrt(C(1,1) + C(2,2)); % drms
            CEP95_4(1, i) = 2 * sigma_h_4(1, i); % CEP

        end
    case false
        fprintf('stopping kalman iteration Task 4 \n');
end

% Task 4 can be improved by using WNLS (delta_u_k=inv(H_k'*inv_R*H_k)*H_k'*inv_R*delta_ro_k) instead of NLS to improve the
% accuracy of u_hat


%% PLOT KALMAN TRACKING Task 4
switch task4Switch
    case true
        figure
        plot(predictions(1,:), predictions(2,:),'g--o');
        hold on
        scatter(points_x, points_y,100, 'b','^'); %APs
        text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
        hold on
        plot(u_hat_trajectory(:,1), u_hat_trajectory(:,2),'r-*');
        legend('Predicted UE location','AP','Real UE location','Location','bestoutside')
        hold on
        xlabel('Ux')
        ylabel('Uy')
        title('Task 4:Predicted UE Position Ux,Uy')
        grid on
        grid minor

        %we assume sampling rate=Ts=1
        timeVector  = linspace(1,size(rhoTraining4,1),size(rhoTraining4,1));

        %plot vx
        figure
        plot(timeVector, predictions(3,:),'c-o');
        hold on
        legend('Predicted Vx')
        hold on
        xlabel('Time')
        ylabel('Vx')
        title('Task 4:Predicted UE Velocity Vx ')
        grid on
        grid minor

        %plot vy
        figure
        plot(timeVector, predictions(4,:),'g-o');
        hold on
        legend('Predicted Vy')
        hold on
        xlabel('Time')
        ylabel('Vy')
        title('Task 4:Predicted UE Velocity Vy ')
        grid on
        grid minor

        %plot both vx,vy but this time w.r.t. time
        figure
        plot(timeVector, predictions(3,:),'c-o');
        hold on
        plot(timeVector, predictions(4,:),'g-o');
        hold on
        legend('Predicted velocity Vx','Predicted velocity Vy')
        hold on
        xlabel('Time')
        ylabel('Velocity')
        title('Task 4:comparison amid Predicted UE Velocity Vx and Vy ')
        grid on
        grid minor
    case false
        fprintf('no plot task 4\n');
end

%% Task 5    
switch task5Switch
    case true
        rhoTraining5 = importdata("GR35/Task5_rhoUEAP_GR35.mat"); % TOA measurements of trojectories
        % INITIALIZATION 
        Ts = 1; % in seconds. time sampling

        F = [eye(2)     , Ts*eye(2);
             zeros(2,2) ,    eye(2)];
        L = [0.5*Ts^2*eye(2);Ts*eye(2)];

        sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)
        P_t_minus_t_minus = eye(4);% cov(trajectory); % initializing P

        Q = sigma_squared_a * (L * L'); % finding Q

        x_hat_t_minus_t_minus = zeros(4,1); 
        predictions= zeros(4, size(rhoTraining5,1)); % used for plotting

        u_hat_trajectory = zeros(size(rhoTraining5,1), 2);

        inv_R = inv(R); % calculating only one time, the inverse of R

        % defining C
        C_stored_5 = cell(1,size(rhoTraining5,1));
        sigma_h_5 = zeros(1,size(rhoTraining5,1));
        CEP95_5 = zeros(1,size(rhoTraining5,1));
        % END initialization kalman filter task 5

        u_hat=zeros(2,1); % initialization of u_hat coordinates
        % not the same H of kalman filter(which is a 8*4)
        for i = 1:(size(rhoTraining5,1))
            for w=1:4 % max_iterations
                delta_ro_k=zeros(Number_of_APs,1);
                H_k=zeros(8,2);
                for j = 1 : Number_of_APs   
                    % COMPUTATION          
                    syms ux_hat uy_hat; % computing H through derivatives
                    f = sqrt((AP(j,1) - ux_hat)^2 + (AP(j,2) - uy_hat)^2);
                    a_dx = diff(f, ux_hat);
                    a_dy = diff(f, uy_hat);
                    H_k(j,1) = vpa(subs(a_dx, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                    H_k(j,2) = vpa(subs(a_dy, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                    delta_ro_j_k=rhoTraining5(i,j)-sqrt((AP(j,1) - u_hat(1,1))^2 + (AP(j,2) - u_hat(2,1))^2);
                    delta_ro_k(j,1)=delta_ro_j_k;
                end
                %INVERSION
                delta_u_k=inv(H_k'*H_k)*H_k'*delta_ro_k;
                %UPDATE SOLUTION
                u_hat=u_hat+delta_u_k;

                %saving u_hat for plotting
                u_hat_trajectory(i,:) = u_hat';
            end

            %EKF code using x_hat as current location 

            % PREDICTION SECTION OF KALMAN FILTER
            % assuming F, Q constant

            if (i == 1)
                x_hat_t_minus_t_minus = u_hat; % still missing the velocity part

                x_hat_t_minus_t_minus(3,1) = 0; % we suppose velocity on y axis and x axis is = 0 since we have no info
                x_hat_t_minus_t_minus(4,1) = 0; 

                x_hat_t_t_minus = x_hat_t_minus_t_minus;
                P_t_t_minus = P_t_minus_t_minus;    
            else
                x_hat_t_t_minus = F * x_hat_t_minus_t_minus; % remember to update at the end x_hat_t_minus_t_minus
                P_t_t_minus = F * P_t_minus_t_minus * F' + Q; % remember to update at the end P_t_minus_t_minus
            end

            predictions(:,i)=x_hat_t_t_minus(:,1); % for plotting
            %UPDATE SECTION OF KALMAN FILTER

            % CALCULATING H
            H = zeros(Number_of_APs,4);
            h = zeros(Number_of_APs,1);
            for j = 1 : Number_of_APs

                h(j,1) = sqrt((AP(j,1) - x_hat_t_t_minus(1,1))^2 + (AP(j,2) - x_hat_t_t_minus(2,1))^2);

                syms ux uy; % computing H through derivatives
                f = sqrt((AP(j,1) - ux)^2 + (AP(j,2) - uy)^2);
                a_dx = diff(f, ux);
                a_dy = diff(f, uy);
                a_ux_final = vpa(subs(a_dx, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);
                a_uy_final = vpa(subs(a_dy, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);

                H(j,1) = a_ux_final; % no minus because we use the approach of the derivatives.
                H(j,2) = a_uy_final;
                H(j,3) = 0; % h is constant w.r.t velocity => derivative w.r.t velocity = 0
                H(j,4) = 0; % same here
            end

            inv_P = inv(P_t_t_minus);

            tmp = (H'* inv_R * H); % optimizing inverse operation. that's why we use \ instead of inv()
            tempA = inv(inv_P + tmp);

            % updating variables
            x_hat_t_t = x_hat_t_t_minus + tempA * H' * inv_R * (rhoTraining5(i, :)' - h);
            P_t_t = inv(inv_P + tmp);

            x_hat_t_minus_t_minus = x_hat_t_t;
            P_t_minus_t_minus = P_t_t;

            %computing performance metrics
            %CALCULATING C
            C = inv(H(:,1:2)'*inv_R*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies are different among themselves;
            C_stored_5(1, i) = mat2cell(C,2); % storing C
            sigma_h_5(1, i) = sqrt(C(1,1) + C(2,2)); % drms
            CEP95_5(1, i) = 2 * sigma_h_5(1, i); % CEP

            % CRITICAL COMPARISON COMMENT:
            % w.r.t Task 4 we noticed by plotting the trajectories that the
            % trajectory of task 5 is more or less the same of the task 4, but with
            % greater variability on the distances between consecutive points (probably due
            % to the missing information about the velocity)

        end  
    case false
        fprintf('stopping kalman iteration Task 5 \n');    
end

%% PLOT KALMAN TRACKING Task 5
switch task5Switch
    case true
        figure
        plot(predictions(1,:), predictions(2,:),'g--o');
        hold on
        scatter(points_x, points_y,100, 'b','^'); %APs
        text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
        hold on
        plot(u_hat_trajectory(:,1), u_hat_trajectory(:,2),'r-*');
        legend('Predicted UE location','AP','Real UE location','Location','bestoutside')
        hold on
        xlabel('Ux')
        ylabel('Uy')
        title('Task 5:Predicted UE Position Ux,Uy')
        grid on
        grid minor

        %we assume sampling rate=Ts=1
        timeVector  = linspace(1,size(rhoTraining5,1),size(rhoTraining5,1));

        %plot vx
        figure
        plot(timeVector, predictions(3,:),'c-o');
        hold on
        legend('Predicted Vx')
        hold on
        xlabel('Time')
        ylabel('Vx')
        title('Task 5:Predicted UE Velocity Vx ')
        grid on
        grid minor

        %plot vy
        figure
        plot(timeVector, predictions(4,:),'g-o');
        hold on
        legend('Predicted Vy')
        hold on
        xlabel('Time')
        ylabel('Vy')
        title('Task 5:Predicted UE Velocity Vy ')
        grid on
        grid minor

        %plot both vx,vy but this time w.r.t. time
        figure
        plot(timeVector,predictions(3,:),'c-o');
        hold on
        plot(timeVector,predictions(4,:),'g-o');
        hold on
        legend('Predicted velocity Vx','Predicted velocity Vy')
        hold on
        xlabel('Time')
        ylabel('Velocity')
        title('Task 5:comparison amid Predicted UE Velocity Vx and Vy ')
        grid on
        grid minor
    case false
        fprintf('no plot of task 5 \n');
end




%% Task 6
rhoTraining6 = importdata("GR35/Task6_rhoUEAP_GR35.mat");
% allora se elemento in rhoTraining(1, :) è NaN, togli il rispettivo AP da
% H e rho e fai la computazione NLS (e.g. AP4 e AP5 hanno TOA NaN => H diventa una 6x4 e rho 6x1)) e anche nel Kalman
%isnan(rhoTraining6(20,2)) % To check if element is NaN

% INITIALIZATION 
Ts = 1; % in seconds. time sampling

F = [eye(2)     , Ts*eye(2);
     zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2);Ts*eye(2)];

sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)
P_t_minus_t_minus = eye(4);% cov(trajectory); % initializing P

Q = sigma_squared_a * (L * L'); % finding Q

x_hat_t_minus_t_minus = zeros(4,1); 
predictions= zeros(4, size(rhoTraining6,1)); % used for plotting

u_hat_trajectory = zeros(size(rhoTraining6,1), 2);


missing_info=[];
% defining C
C_stored_6 = cell(1,size(rhoTraining6,1));
sigma_h_6 = zeros(1,size(rhoTraining6,1));
CEP95_6 = zeros(1,size(rhoTraining6,1));
% END initialization kalman filter task 5
invalidRho=0;
u_hat=zeros(2,1); % initialization of u_hat coordinates
% not the same H of kalman filter(which is a 8*4)
for i = 1:(size(rhoTraining6,1))
    sizeValidMeasurements=0;
    indexValidMeasurements=zeros(1,Number_of_APs);
    for w=1:4 % max_iterations
        invalidRho=0;
        delta_ro_k=zeros(1,1); %this time delta_ro_k is by default of size 1 and could be of size 8 only if all the 8 rhos are not NaN
        H_k=zeros(1,2); %this time H_k is by default of size 1 and could be of size 8 only if all the 8 rhos are not NaN
        firstValidValue=false;
        for j = 1 : Number_of_APs
            if(isnan(rhoTraining6(i,j))==false)
                if(firstValidValue==true)
                    size_H_k=size(H_k,1)+1;
                else
                   size_H_k=1;
                end
                firstValidValue=true;
                indexValidMeasurements(1,j)=j;
                % COMPUTATION          
                syms ux_hat uy_hat; % computing H through derivatives
                f = sqrt((AP(j,1) - ux_hat)^2 + (AP(j,2) - uy_hat)^2);
                a_dx = diff(f, ux_hat);
                a_dy = diff(f, uy_hat);
                H_k(size_H_k,1) = vpa(subs(a_dx, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                H_k(size_H_k,2) = vpa(subs(a_dy, [ux_hat, uy_hat], [u_hat(1,1), u_hat(2,1)]), 5);
                delta_ro_j_k=rhoTraining6(i,j)-sqrt((AP(j,1) - u_hat(1,1))^2 + (AP(j,2) - u_hat(2,1))^2);
                delta_ro_k(size_H_k,1)=delta_ro_j_k;
            else
                invalidRho=invalidRho+1;
            end
        end

        %INVERSION
        if(size_H_k>1 && invalidRho<Number_of_APs)
            delta_u_k=inv(H_k'*H_k)*H_k'*delta_ro_k;
            %UPDATE SOLUTION
            u_hat=u_hat+delta_u_k;
        %if size==1 delta_u_k explodes-->I can say that for that point I do not have enough information->I cannot update x_hat
        end
    end
    %saving u_hat for plotting
    u_hat_trajectory(i,:) = u_hat';
    sizeValidMeasurements=size(H_k,1);

    if(invalidRho==Number_of_APs || invalidRho==Number_of_APs-1)
        tempMatrix(1:2,1:2)= -1;
        C_stored_6(1, i) = mat2cell(tempMatrix);
        sigma_h_6(1, i) = -1;
        CEP95_6(1, i) = -1;
        fprintf('not enough TOA measurements available for this step->no prediction can be done: # valid %d \n',Number_of_APs-invalidRho);
        predictions([1 2],i)=u_hat; %riassegno come prediction ux,uy quella dello step precedente e con vx,vy=0
    else    
        %EKF code using x_hat as current location 

        % PREDICTION SECTION OF KALMAN FILTER
        % assuming F, Q constant

        if (i == 1)
            x_hat_t_minus_t_minus = u_hat; % still missing the velocity part
            x_hat_t_minus_t_minus(3,1) = 0; % we suppose velocity on y axis and x axis is = 0 since we have no info
            x_hat_t_minus_t_minus(4,1) = 0;   
            x_hat_t_t_minus = x_hat_t_minus_t_minus;
            P_t_t_minus = P_t_minus_t_minus;    
        else
            x_hat_t_t_minus = F * x_hat_t_minus_t_minus; % remember to update at the end x_hat_t_minus_t_minus
            P_t_t_minus = F * P_t_minus_t_minus * F' + Q; % remember to update at the end P_t_minus_t_minus
        end

        predictions(:,i)=x_hat_t_t_minus(:,1); % for plotting
        %UPDATE SECTION OF KALMAN FILTER

        % CALCULATING H
        H = zeros(sizeValidMeasurements,4);
        h = zeros(sizeValidMeasurements,1);
        rhoValid=zeros(1,sizeValidMeasurements);
        sizeRhoValid=1;
        R_Partial=zeros(sizeValidMeasurements,sizeValidMeasurements);
        for j = 1 : Number_of_APs
            if(indexValidMeasurements(1,j)>0)
                h(sizeRhoValid,1) = sqrt((AP(j,1) - x_hat_t_t_minus(1,1))^2 + (AP(j,2) - x_hat_t_t_minus(2,1))^2);
                rhoValid(1,sizeRhoValid)=rhoTraining6(i, j);
                syms ux uy; % computing H through derivatives
                f = sqrt((AP(j,1) - ux)^2 + (AP(j,2) - uy)^2);
                a_dx = diff(f, ux);
                a_dy = diff(f, uy);
                a_ux_final = vpa(subs(a_dx, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);
                a_uy_final = vpa(subs(a_dy, [ux, uy], [u_hat(1,1), u_hat(2,1)]), 5);

                H(sizeRhoValid,1) = a_ux_final; % no minus because we use the approach of the derivatives.
                H(sizeRhoValid,2) = a_uy_final;
                H(sizeRhoValid,3) = 0; % h is constant w.r.t velocity => derivative w.r.t velocity = 0
                H(sizeRhoValid,4) = 0; % same here

                R_Partial(sizeRhoValid,sizeRhoValid)=R(j,j);
                sizeRhoValid=sizeRhoValid+1;
            end

        end
        
        inv_R_Partial=inv(R_Partial);

        inv_P = inv(P_t_t_minus);
        tmp = (H'* inv_R_Partial * H); % optimizing inverse operation. that's why we use \ instead of inv()
        tempA = inv(inv_P + tmp);

        % updating variables
        x_hat_t_t = x_hat_t_t_minus + tempA * H' * inv_R_Partial * (rhoValid' - h);
        P_t_t = inv(inv_P + tmp);
        x_hat_t_minus_t_minus = x_hat_t_t;
        P_t_minus_t_minus = P_t_t;

         
        % H,R attenzione
        % computing performance metrics
        % CALCULATING C
        C = inv(H(:,1:2)'*inv_R_Partial*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies are different among themselves;
        C_stored_6(1, i) = mat2cell(C,2); % storing C
        sigma_h_6(1, i) = sqrt(C(1,1) + C(2,2)); % drms
        CEP95_6(1, i) = 2 * sigma_h_6(1, i); % CEP
        

        % CRITICAL COMPARISON COMMENT:
        % w.r.t Task 4 we noticed by plotting the trajectories that the
        % trajectory of task 5 is more or less the same of the task 4, but with
        % greater variability on the distances between consecutive points (probably due
        % to the missing information about the velocity)
    end
    
    %TODO: Add the switch 
    if(mod(i, 3)== 0) % stopping before completing
        fprintf('stopping kalman iteration Task 6\n');
        break;
    end
       
end

%% PLOT KALMAN TRACKING Task 6
switch task6Switch    
    case true
        plot_dir(predictions(1,:)', predictions(2, :)', i);
        
        figure
        plot(predictions(1,:), predictions(2,:),'g--o');
        hold on
        scatter(points_x, points_y,100, 'b','^'); %APs
        text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
        hold on
        plot(u_hat_trajectory(:,1), u_hat_trajectory(:,2),'r-*');
        legend('Predicted UE location','AP','Real UE location','Location','bestoutside')
        hold on
        xlabel('Ux')
        ylabel('Uy')
        title('Task 6:Predicted UE Position Ux,Uy')
        grid on
        grid minor

        %we assume sampling rate=Ts=1
        timeVector  = linspace(1,size(rhoTraining6,1),size(rhoTraining6,1));

        %plot vx
        figure
        plot(timeVector, predictions(3,:),'c-o');
        hold on
        legend('Predicted Vx')
        hold on
        xlabel('Time')
        ylabel('Vx')
        title('Task 6:Predicted UE Velocity Vx ')
        grid on
        grid minor

        %plot vy
        figure
        plot(timeVector, predictions(4,:),'g-o');
        hold on
        legend('Predicted Vy')
        hold on
        xlabel('Time')
        ylabel('Vy')
        title('Task 6:Predicted UE Velocity Vy ')
        grid on
        grid minor

        %plot both vx,vy but this time w.r.t. time
        figure
        plot(timeVector, predictions(3,:),'c-o');
        hold on
        plot(timeVector, predictions(4,:),'g-o');
        hold on
        legend('Predicted velocity Vx','Predicted velocity Vy')
        hold on
        xlabel('Time')
        ylabel('Velocity')
        title('Task 6:comparison amid Predicted UE Velocity Vx and Vy ')
        grid on
        grid minor
    case false
        fprintf('no plot of task 6 \n');
end

