function u_hat = test3NLS(Number_of_APs,AP,measurement_of_current_trajectory,row_)
   
    u_hat=zeros(2,1); 
    for w=1:4 
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
            delta_ro_j_k=measurement_of_current_trajectory(row_,j)-sqrt((AP(j,1) - u_hat(1,1))^2 + (AP(j,2) - u_hat(2,1))^2);
            delta_ro_k(j,1)=delta_ro_j_k;
        end
        %INVERSION
        delta_u_k=inv(H_k'*H_k)*H_k'*delta_ro_k;
        %UPDATE SOLUTION
        u_hat=u_hat+delta_u_k;
    end
end