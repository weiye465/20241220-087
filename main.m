% 主程序：海洋单桩结构动力响应分析
clear all; close all; clc;

% 创建日志文件
fid = fopen('calculation_log.txt', 'w');
fprintf(fid, '=== 海洋单桩结构动力响应分析日志 ===\n\n');

% 声明全局变量
global H omega k D CM CD rho h T

% 1. 加载基本参数
parameters;
fprintf(fid, '基本参数：\n');
fprintf(fid, '水深 h = %.2f m\n', h);
fprintf(fid, '质量 m = %.2f kg\n', m);
fprintf(fid, '直径 D = %.2f m\n', D);
fprintf(fid, '抗弯刚度 EI = %.2e N·m²\n', EI);

fprintf(fid, '基本参数：\n');
fprintf(fid, '水深 h = %.2f m\n', h);
fprintf(fid, '质量 m = %.2f kg\n', m);
fprintf(fid, '直径 D = %.2f m\n', D);
fprintf(fid, '抗弯刚度 EI = %.2e N·m²\n', EI);

fprintf(fid, '基本参数：\n');
fprintf(fid, '水深 h = %.2f m\n', h);
fprintf(fid, '质量 m = %.2f kg\n', m);
fprintf(fid, '直径 D = %.2f m\n', D);
fprintf(fid, '抗弯刚度 EI = %.2e N·m²\n', EI);

% 计算波浪参数
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
fprintf(fid, '\n波浪参数：\n');
fprintf(fid, '角频率 omega = %.4f rad/s\n', omega);
fprintf(fid, '波数 k = %.4f m^-1\n', k);
fprintf(fid, '迭代次数 = %d\n', i);

% 2. 结构离散化
discretization;
fprintf(fid, '\n结构离散化：\n');
fprintf(fid, '单元数量 = %d\n', n_elements);
fprintf(fid, '单元长度 = %.3f m\n', element_length);
fprintf(fid, '节点数量 = %d\n', nodes);

% 3. 组装总体矩阵
% 初始化总体矩阵
M = zeros(2*nodes, 2*nodes);
K = zeros(2*nodes, 2*nodes);

% 计算单元矩阵
me = element_length/420 * [... 
    156, 22*element_length, 54, -13*element_length;
    22*element_length, 4*element_length^2, 13*element_length, -3*element_length^2;
    54, 13*element_length, 156, -22*element_length;
    -13*element_length, -3*element_length^2, -22*element_length, 4*element_length^2];

ke = EI/(element_length^3) * [...
    12, 6*element_length, -12, 6*element_length;
    6*element_length, 4*element_length^2, -6*element_length, 2*element_length^2;
    -12, -6*element_length, 12, -6*element_length;
    6*element_length, 2*element_length^2, -6*element_length, 4*element_length^2];

% 组装总体矩阵
for e = 1:n_elements
    % 获取单元节点编号
    node1 = e;
    node2 = e + 1;
    
    % 计算自由度编号
    dofs = [2*node1-1, 2*node1, 2*node2-1, 2*node2];
    
    % 组装质量矩和刚度矩阵
    M(dofs, dofs) = M(dofs, dofs) + me;
    K(dofs, dofs) = K(dofs, dofs) + ke;
end

% 4. 施加边界条件（底部固定）
fixed_dofs = [1, 2];  % 底部节点的自由度
free_dofs = 3:2*nodes;  % 其余自由度

% 修改矩阵以应用边界条件
M = M(free_dofs, free_dofs);
K = K(free_dofs, free_dofs);

% 5. 时间步长设置
dt = T/20;
t_end = 20*T;
time = 0:dt:t_end;
nt = length(time);

% 6. 记录特定位置的响应
% 找到最接近指定高度的节点
z_targets = [0, -(1/4)*h, -(1/2)*h, -(3/4)*h];  % 从上到下排列
target_nodes = zeros(1, length(z_targets));
for i = 1:length(z_targets)
    [~, target_nodes(i)] = min(abs(z_coords - z_targets(i)));
end

% 初始化响应记录数组
displacement = zeros(length(z_targets), nt);
velocity = zeros(length(z_targets), nt);
acceleration = zeros(length(z_targets), nt);

% 7. Newmark时间积分
beta = 1/4;
gamma = 1/2;

% 初始条件
u = zeros(length(free_dofs), 1);
v = zeros(length(free_dofs), 1);
a = zeros(length(free_dofs), 1);

% 时间步推进
max_disp = zeros(length(z_targets), 1);
max_vel = zeros(length(z_targets), 1);
max_acc = zeros(length(z_targets), 1);

% 记录波浪力
wave_forces = zeros(length(z_targets), nt);

for it = 1:nt
    t = time(it);
    
    % 计算波浪力
    F = calculate_wave_forces(t, z_coords, free_dofs);
    
    % 记录特定位置的波浪力
    for i = 1:length(target_nodes)
        if z_coords(target_nodes(i)) <= 0 && z_coords(target_nodes(i)) >= -h
            z = z_coords(target_nodes(i));
            u = H/2 * omega * cosh(k*(z+h))/sinh(k*h) * cos(omega*t);
            du_dt = H/2 * omega^2 * cosh(k*(z+h))/sinh(k*h) * sin(omega*t);
            % 直接计算 Morison 力
            f_inertia = CM * rho * pi * D^2/4 * du_dt;
            f_drag = CD * rho * D/2 * abs(u) * u;
            wave_forces(i,it) = f_inertia + f_drag;
        end
    end
    
    % Newmark积分
    [u, v, a] = newmark_step(M, K, F, u, v, a, dt, beta, gamma);
    
    % 记录特定位置的响应
    for i = 1:length(target_nodes)
        node_dof = 2*(target_nodes(i)-1) - 1;  % 节点的水平位移自由度
        if node_dof > 0 && node_dof <= length(u)
            displacement(i,it) = u(node_dof);
            velocity(i,it) = v(node_dof);
            acceleration(i,it) = a(node_dof);
            
            % 更新最大响应值
            max_disp(i) = max(max_disp(i), abs(u(node_dof)));
            max_vel(i) = max(max_vel(i), abs(v(node_dof)));
            max_acc(i) = max(max_acc(i), abs(a(node_dof)));
        end
    end
end

% 输出最大响应值
fprintf(fid, '\n最大响应值：\n');
for i = 1:length(z_targets)
    fprintf(fid, '\nz = %.2f m 处：\n', z_targets(i));
    fprintf(fid, '最大位移 = %.6f m\n', max_disp(i));
    fprintf(fid, '最大速度 = %.6f m/s\n', max_vel(i));
    fprintf(fid, '最大加速度 = %.6f m/s²\n', max_acc(i));
end

% 8. 绘制结果
plot_results(time, displacement, velocity, acceleration, z_targets);
% 绘制波浪力
plot_wave_forces(time, wave_forces, z_targets, h);

% 关闭日志文件
fclose(fid); 

% 添加 Morison 方程计算函数
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