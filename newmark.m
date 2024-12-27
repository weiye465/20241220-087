% Newmark-β参数
beta = 1/4;  % 常用值
gamma = 1/2;  % 常用值
dt = T/20;    % 时间步长
t_end = 20*T; % 总计算时间

% 初始条件
u = zeros(2*nodes, 1);
v = zeros(2*nodes, 1);
a = zeros(2*nodes, 1);

% 时间步推进
for t = 0:dt:t_end
    % 1. 计算波浪力
    F = calculate_wave_forces(t);
    
    % 2. 预测位移和速度
    u_pred = u + dt*v + dt^2/2*((1-2*beta)*a);
    v_pred = v + dt*((1-gamma)*a);
    
    % 3. 求解加速度
    K_eff = K + 1/(beta*dt^2)*M;
    F_eff = F - K*(u_pred - u) - M*v_pred;
    a_new = K_eff\F_eff;
    
    % 4. 修正位移和速度
    u = u_pred + beta*dt^2*a_new;
    v = v_pred + gamma*dt*a_new;
    a = a_new;
    
    % 5. 记录指定位置的响应
    record_response(t, u, v, a);
end 