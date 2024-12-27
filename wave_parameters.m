% 波浪参数计算
g = 9.81;          % 重力加速度
omega = 2*pi/T;    % 角频率
k = omega^2/g;     % 初始波数估计

% 波数迭代求解(色散关系)
for i = 1:10
    k_new = omega^2/(g*tanh(k*h));
    if abs(k_new - k) < 1e-6
        break;
    end
    k = k_new;
end

% 波面升高函数
function eta = get_wave_elevation(x, t)
    eta = H/2 * cos(k*x - omega*t);
end 