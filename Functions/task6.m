function task6(Number_of_APs,AP,Q,inv_R,F,points_x,points_y,AP_IDs,R) 
    rhoTraining6 = importdata("GR35/Task6_rhoUEAP_GR35.mat");
    % allora se elemento in rhoTraining(1, :) è NaN, togli il rispettivo AP da
    % H e rho e fai la computazione NLS (e.g. AP4 e AP5 hanno TOA NaN => H diventa una 6x4 e rho 6x1)) e anche nel Kalman
    %isnan(rhoTraining6(20,2)) % To check if element is NaN

    % INITIALIZATION 
    P_t_minus_t_minus = eye(4);% cov(trajectory); % initializing P

    x_hat_t_minus_t_minus = zeros(4,1); 
    predictions= zeros(4, size(rhoTraining6,1)); % used for plotting

    u_hat_trajectory = zeros(size(rhoTraining6,1), 2);

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
%             tempMatrix(1:2,1:2)= -1;
%             C_stored_6(1, i) = mat2cell(tempMatrix); non gira FIXME
%             sigma_h_6(1, i) = -1;
%             CEP95_6(1, i) = -1;
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
%             C = inv(H(:,1:2)'*inv_R_Partial*H(:,1:2)); % calculating lower bound since: R not equal to sigma * I => accuracies are different among themselves;
%             C_stored_6(1, i) = mat2cell(C,2); % storing C
%             sigma_h_6(1, i) = sqrt(C(1,1) + C(2,2)); % drms
%             CEP95_6(1, i) = 2 * sigma_h_6(1, i); % CEP


            % CRITICAL COMPARISON COMMENT:
            % w.r.t Task 4 we noticed by plotting the trajectories that the
            % trajectory of task 5 is more or less the same of the task 4, but with
            % greater variability on the distances between consecutive points (probably due
            % to the missing information about the velocity)
        end

        %TODO:remove this break 
%         if(mod(i, 3)== 0) % stopping before completing
%             fprintf('stopping in the middle kalman iteration Task 6\n');
%             break;
%         end

    end
    
    %% PLOT KALMAN TRACKING Task 6
    figure
    plot_dir(predictions(1,:)', predictions(2, :)', i,'g--o','p');
    hold on
    scatter(points_x, points_y,100, 'b','^'); %APs
    text(points_x, points_y,AP_IDs,'Color','blue','FontSize',12)
    hold on
    plot_dir(u_hat_trajectory(:,1), u_hat_trajectory(:,2), i,'r-*','b');
    legend('Predicted UE location','Predicted movement','AP','Real UE location','Real movement','Location','bestoutside')
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
end