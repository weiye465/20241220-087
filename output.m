% 绘制响应时程
figure;
subplot(3,1,1);
plot(time, displacement);
title('位移响应');

subplot(3,1,2);
plot(time, velocity);
title('速度响应');

subplot(3,1,3);
plot(time, acceleration);
title('加速度响应'); 