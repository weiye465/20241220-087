function F = calculate_wave_forces(t, z_coords, free_dofs)
    % 声明全局变量
    global H omega k D CM CD rho h T
    
    persistent force_log_fid
    if isempty(force_log_fid)
        force_log_fid = fopen('wave_forces_log.txt', 'w');
        fprintf(force_log_fid, '=== 波浪力计算日志 ===\n\n');
    end
    
    F = zeros(length(free_dofs), 1);
    
    % 记录当前时刻
    if mod(t, T) == 0
        fprintf(force_log_fid, '\n时刻 t = %.2f s:\n', t);
    end
    
    for i = 1:length(z_coords)
        if z_coords(i) <= 0 && z_coords(i) >= -h
            z = z_coords(i);
            % 水质点运动
            u = H/2 * omega * cosh(k*(z+h))/sinh(k*h) * cos(omega*t);
            du_dt = -H/2 * omega^2 * cosh(k*(z+h))/sinh(k*h) * sin(omega*t);
            
            force = calculate_morison_force(u, du_dt, D);
            
            % 每个周期开始时记录典型位置的波浪力
            if mod(t, T) == 0 && (abs(z) == 0 || abs(z+h/4) < 0.1 || abs(z+h/2) < 0.1 || abs(z+3*h/4) < 0.1)
                fprintf(force_log_fid, 'z = %.2f m: 速度 = %.4f m/s, 加速度 = %.4f m/s², 波浪力 = %.4f N/m\n', ...
                    z, u, du_dt, force);
            end
            
            dof = 2*i-1;
            if dof > 2
                F(dof-2) = F(dof-2) + force;
            end
        end
    end
end

% 在同一个文件中定义 Morison 方程计算函数
function f = calculate_morison_force(u, du_dt, D)
    % 声明全局变量
    global CM CD rho
    
    % u: 水质点速度
    % du_dt: 水质点加速度
    % 计算单位长度上的波浪荷载
    f_inertia = CM * rho * pi * D^2/4 * du_dt;
    f_drag = CD * rho * D/2 * abs(u) * u;
    f = f_inertia + f_drag;
end 