ode = @(t,y) y;
y0 = 1;
T = 10;
dt = 1;

figure();
hold on;

[t_course, y_course] = fweuler(ode, y0, [0, T], dt);
plot(t_course, y_course, 'DisplayName', 'Course');

for k=1:numel(t_course)-1
  [~, y_fine] = fweuler(ode, y_course(k), t_course(k:k+1)', 1/1000);
  y_course(k+1) = y_fine(end);
end

plot(t_course, y_course, '--', 'DisplayName', '1-Parareal');
plot(t_course, exp(t_course), 'DisplayName', 'True');

hold off;
title('One Parareal comparison'); ylabel('y(t)'); xlabel('t in [0,1]')
legend();
