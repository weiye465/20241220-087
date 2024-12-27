function plot_wave_forces(time, wave_forces, z_targets, h)
    figure('Position', [150 150 800 400]);
    
    % 绘制波浪力时程
    plot(time, wave_forces, 'LineWidth', 1.5);
    title('波浪力时程');
    xlabel('时间 (s)');
    ylabel('波浪力 (N/m)');
    legend(sprintf('z = %.2f m', z_targets(1)), ...
           sprintf('z = %.2f m', z_targets(2)), ...
           sprintf('z = %.2f m', z_targets(3)), ...
           sprintf('z = %.2f m', z_targets(4)));
    grid on;
    
    % 绘制某一时刻的波浪力分布
    figure('Position', [200 200 400 600]);
    t_idx = floor(length(time)/4);  % 选择1/4周期时刻
    plot(wave_forces(:,t_idx), z_targets, 'ro-', 'LineWidth', 1.5);
    set(gca, 'YDir', 'reverse');  % 让y轴向下为正
    title(sprintf('t = %.2f s时刻的波浪力分布', time(t_idx)));
    xlabel('波浪力 (N/m)');
    ylabel('水深 (m)');
    grid on;
    ylim([-h 0]);  % 设置y轴范围从水面(0)到海底(-h)
end 