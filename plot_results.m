function plot_results(time, displacement, velocity, acceleration, z_targets)
    % 绘制结果
    figure('Position', [100 100 800 600]);
    
    % 位移响应
    subplot(3,1,1);
    plot(time, displacement);
    title('位移响应');
    xlabel('时间 (s)');
    ylabel('位移 (m)');
    legend(sprintf('z = %.2f m', z_targets(1)), ...
           sprintf('z = %.2f m', z_targets(2)), ...
           sprintf('z = %.2f m', z_targets(3)), ...
           sprintf('z = %.2f m', z_targets(4)));
    grid on;
    
    % 速度响应
    subplot(3,1,2);
    plot(time, velocity);
    title('速度响应');
    xlabel('时间 (s)');
    ylabel('速度 (m/s)');
    grid on;
    
    % 加速度响应
    subplot(3,1,3);
    plot(time, acceleration);
    title('加速度响应');
    xlabel('时间 (s)');
    ylabel('加速度 (m/s^2)');
    grid on;
end 