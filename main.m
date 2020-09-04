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


%PARAMETERS INITIALIZATION for task 3,4,5,6
sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)
inv_R = inv(R); % calculating only one time, the inverse of R
Ts = 1; % in seconds. time sampling

F = [eye(2)     , Ts*eye(2);
         zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2);Ts*eye(2)];
Q = sigma_squared_a * (L * L'); % finding Q

%% TASK 3 (ONLY FIRST TRAJECTORY)
switch task3Switch
    case true 
        task3(Number_of_APs,AP,trajectories,Q,inv_R,F,points_x,points_y,AP_IDs);
    case false
        fprintf('no kalman filter Task 3\n');
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
        test3NLS(Number_of_APs,AP,trajectories);
    case false
        fprintf('Do not test the trajectory prediction \n');
end

%% Task 4
switch task4Switch
    case true 
        task4(Number_of_APs,AP,Q,inv_R,F,points_x,points_y,AP_IDs);
    case false
        fprintf('no kalman filter Task 4\n');
end

% Task 4 can be improved by using WNLS (delta_u_k=inv(H_k'*inv_R*H_k)*H_k'*inv_R*delta_ro_k) instead of NLS to improve the
% accuracy of u_hat

%% Task 5    
switch task5Switch
    case true
        task5(Number_of_APs,AP,Q,inv_R,F,points_x,points_y,AP_IDs);
    case false
         fprintf('no kalman filter Task 5\n');   
end

%% Task 6
switch task6Switch    
    case true
        task6(Number_of_APs,AP,Q,inv_R,F,points_x,points_y,AP_IDs,R);
    case false
        fprintf('no kalman filter Task 6\n');
end