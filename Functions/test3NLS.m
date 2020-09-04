function test3NLS(Number_of_APs,AP,trajectories)
    ID_trajectory = 5; % choose your trajectory here
    trajectory = cell2mat(trajectories(1,ID_trajectory)); % if you change the trajectory, remember to change also the measurement_of_current_trajectory
    rhoTraining = importdata("GR35/Task3_rhoUEAP_GR35.mat"); % TOA measurements of trajectories
    measurement_of_current_trajectory = cell2mat(rhoTraining(1,ID_trajectory)); %taking TOA of this trajectory
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
end