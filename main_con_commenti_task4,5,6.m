%% Fabio Carminati, Emanuele Gallone LNSM Project tasks 1-6
clear all, clc, close all
set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)

addpath Functions

%% scenario settings (4x4 Km)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;
% With these switches enable a custom run of the code(i.e. decide each time which tasks to run)
% Since tasks 1,2 are mandatory for other tasks they are always computed and their switches refers just to their plot 

task1aSwitch=false;
task2Switch=false;
task3Switch=false;
task3TestingSwitch=true; 
task4Switch=false;
task5Switch=false;
task6Switch=true;
%% importing data
AP_Data = importdata("GR35/Task1a_rhoUEAP.mat");
Number_of_APs = size(AP_Data, 1); % number of APs is mapped in row length
u = [0, 0]; % reference point

%% Task_1a: find position of the Access Points
%The AOA and TOA measurements for each AP are placed in a system in order to find the
%its coordinates
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
    
    %define the system
    syms Sx Sy 
    eq_toa = sqrt((Sx -u(1))^2 + (Sy - u(2))^2) == toa_measurement;
    a = (u(2) - Sy) / (u(1) - Sx);
    eq_aoa =atan(a) == aoa_measurement;    
    solution = vpasolve([eq_toa, eq_aoa], [Sx, Sy]);
    
    SxSol = round(solution.Sx);
    points_x = [points_x sign*SxSol];
 
    SySol = round(solution.Sy);
    points_y = [points_y sign*SySol];
    
    %assign coordinates
    AP(i, 1) = sign*SxSol; 
    AP(i, 2) = sign*SySol;
       
end
AP_IDs={'AP1','AP2','AP3','AP4','AP5','AP6','AP7','AP8'};

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
        

%% Task_1b:find covariance matrix of TOA measurements (R)
u = [500 -800];
rho = importdata("GR35/Task1b_rhoUEAP.mat");
sigma_squared = zeros(1,8);

for i = 1:Number_of_APs
    temp = [];
    distance_u_ap = sqrt((AP(i,1) - u(1))^2 + (AP(i,2) - u(2))^2);
    for j = 1:size(rho, 1) 
        n = rho(j, i) - distance_u_ap;
        temp = [temp n]; 
    end
    sigma_squared(i) = cov(temp);

end

R = diag(sigma_squared);
R = round(R); % covariance matrix

%% Task 2 motion models statistics and parameters
% We observe the trajectory dataset and since
%1)Ux,Ux,Vx,Vy are given
%2)There is a random acceleration process whose mean  is 0 while the variance is not 0  
% we came to the conclusion that this is a 
% M3 motion model: Random Force model
%{
    TotalSimulationTime = 20; %s
    F = [eye(2)     , Ts*eye(2);
         zeros(2,2) ,    eye(2)];
    L = [0.5*Ts^2*eye(2);Ts*eye(2)];

    v = zeros(2,TotalSimulationTime);
    v(:,1) = [7;2];
    x = [zeros(2,TotalSimulationTime);v];
    x(:,t) = F*x(:,t-1) + L*w;
%}

trajectories = importdata("GR35/Task2_trajectory_GR35.mat");
Ts = 1; %samplign interval in seconds
wx = [0]; %Vector where we placed the random acceleration measurements over x 
wy = [0]; %Vector where we placed the random acceleration measurements over y
for i = 1:100 %100 is the total number of trajectories
    measures = cell2mat(trajectories(1,i));
    
    for temp = 2:size(measures,1) %compute wx,t and vy,t vectors for each measure of that trajectory
        vx_t_minus_one = measures(temp - 1, 3);
        vx_t = measures(temp, 3);
        
        vy_t_minus_one = measures(temp - 1, 4);
        vy_t = measures(temp, 4);
          
        wx = [wx ((vx_t - vx_t_minus_one) / Ts)];
        wy = [wy ((vy_t - vy_t_minus_one) / Ts)];
    end
end
%Compute mean and variance to check whether they are compliant with M3 
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
        vector_legend=[];
        %% plot task 2 trajectories (in this case 100/10=10 in order to be able to clearly dinstinguish each trajectory in the plot)
        figure
        for i = 1:size(trajectories,2)/10
            trajectory = cell2mat(trajectories(i));
            plot(trajectory(:,1), trajectory(:,2));
            hold on;
            text(trajectory(1,1)+10,trajectory(1,2)+10,string(i),'FontSize',11)
        end
        xlabel('Ux[m]')
        ylabel('Uy[m]')
        title('M3 Trajectories')
        grid on
        
        %% plot task 2 velocity variation (w_at)
        figure
        delta_v_x = zeros(1, (size(w_a,2)-1)/100);
        delta_v_y = zeros(1, (size(w_a,2)-1)/100);

        for i = 1:size(w_a, 2)/100
            delta_v_x(i) = w_a(1,i);
            delta_v_y(i) = w_a(2,i);
        end

        plot(1:size(delta_v_x, 2), delta_v_x);
        hold on
        plot(1:size(delta_v_y, 2), delta_v_y);
        hold on
        xlabel('Time')
        ylabel('Variation')
        title('Random acceleration process for the first trajectory')
        legend('Wx,t','Wy,t')
        grid on
        grid minor
case false
        fprintf('No plot task 2 \n');
end


%PARAMETERS INITIALIZATION for task 3,4,5,6
sigma_squared_a = mean([var_x var_y]); % finding scalar sigma (x_var(3,3) is variance of vx while x_var(4,4) is variance of vy)
inv_R = inv(R); % calculating only one time, the inverse of R
Ts = 1; % in seconds. time sampling

F = [eye(2)     , Ts*eye(2);
         zeros(2,2) ,    eye(2)];
L = [0.5*Ts^2*eye(2);Ts*eye(2)];
Q = sigma_squared_a * (L * L'); % finding Q

%% TASK 3:Kalman Filter on known trajectory 
%{ 
Some assumptions:
1)X is initialized with the first line of the chosen trajectory (measurement)
(1 trajectory at a time)
2)P is initialized as an eye matrix 4*4
3)assuming F, Q constant over time
%}
switch task3Switch
    case true 
        task3(Number_of_APs,AP,trajectories,Q,inv_R,F,points_x,points_y,AP_IDs);
    case false
        fprintf('no kalman filter Task 3\n');
 end

%% TEST to check whether our NLS implementation works:

switch task3TestingSwitch
    case true
        test3NLS(Number_of_APs,AP,trajectories);
    case false
        fprintf('Do not test the trajectory prediction \n');
end

%% Task 4:Kalman Filter over an unknown trajectory
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