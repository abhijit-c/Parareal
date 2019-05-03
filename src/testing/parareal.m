method = @(ode, y0, t_range, dt) fweuler(ode, y0, t_range, dt);
y0 = 0.1;
ode = @(t,y) y;
T = 5;
dt = .1;

num_para_its = 5;
num_temporal_its = 50;
U = zeros(num_temporal_its+1, num_para_its+1);
U(1,:) = y0;
[t, U(:,1)] = method(ode, y0, [0, T], dt);

for k=1:num_para_its
  for n=1:num_temporal_its % Correction iteration (Parallel capable)
    [~, G] = method(ode, U(n,k+1), t(n:n+1)', dt);
    [~, F] = method(ode, U(n,k), t(n:n+1)', dt^2);
    U(n+1,k+1) = G(end) + F(end) - U(n,k);
  end
end

[tf, Fine] = method(ode, y0, [0, T], dt^2);
figure();

hold on;
for k=0:num_para_its
  plot(t, U(:,k+1), 'DisplayName', sprintf('%d-parareal',k));
end
plot(tf, Fine, '--', 'DisplayName', 'Full Fine Solution');
plot(t, y0*exp(t), 'k-s', 'DisplayName', 'True Solution');
hold off;

title('Parareal on y'' = y, y(0) = 1'); ylabel('y(t)'); xlabel('t in [0,1]');
legend('Location', 'northeast');
