function [u_new, v_new, a_new] = newmark_step(M, K, F, u, v, a, dt, beta, gamma)
    % Newmark-β单步计算
    u_pred = u + dt*v + dt^2/2*((1-2*beta)*a);
    v_pred = v + dt*((1-gamma)*a);
    
    K_eff = K + 1/(beta*dt^2)*M;
    F_eff = F - K*(u_pred - u) - M*v_pred;
    
    a_new = K_eff\F_eff;
    u_new = u_pred + beta*dt^2*a_new;
    v_new = v_pred + gamma*dt*a_new;
end 