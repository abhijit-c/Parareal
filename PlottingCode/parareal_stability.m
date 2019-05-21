% k = 0
N = 100;
[M, T] = meshgrid(linspace(-5, 5, N), linspace(0, 5, N));

roots = abs( (M.*T + 1) );
root_cond = zeros(N, N);

for j=1:N
  for k=1:N
    root_cond(j,k) = 0;
    if roots(j,k) <= 1
      root_cond(j,k) = 1;
    end
  end
end

figure();
surf(M, T, root_cond); colormap('gray');
view(0, 90);
title('Forward Euler Parareal, $k = 0$, $\delta t = \Delta t^2$', 'Interpreter', 'latex');
xlabel('\mu'); ylabel('\Delta t');
%set(gca, 'fontsize', 16);

% k = 1
N = 100;
[M, T] = meshgrid(linspace(-5, 5, N), linspace(0, 5, N));

roots = (1+M.*T).^2 + 2.*M.*T.*(M.*(T.^3)+2.*T-1).*(1+M.*T);
root_cond = zeros(N, N);

for j=1:N
  for k=1:N
    root_cond(j,k) = 0;
    if roots(j,k) <= 1
      root_cond(j,k) = 1;
    end
  end
end

figure();
surf(M, T, root_cond); colormap('gray');
view(0, 90);
title('Forward Euler Parareal, $k = 1$, $\delta t = \Delta t^2$', 'Interpreter', 'latex');
xlabel('\mu'); ylabel('\Delta t');
%set(gca, 'fontsize', 16);

% k = 2
N = 100;
[M, T] = meshgrid(linspace(-5, 5, N), linspace(0, 5, N));

roots = abs( (M.*(T.^2) + 1).^4 );
root_cond = zeros(N, N);

for j=1:N
  for k=1:N
    root_cond(j,k) = 0;
    if roots(j,k) <= 1
      root_cond(j,k) = 1;
    end
  end
end

figure();
surf(M, T, root_cond); colormap('gray');
view(0, 90);
title('Forward Euler Parareal, $k = 2$, $\delta t = \Delta t^2$', 'Interpreter', 'latex');
xlabel('\mu'); ylabel('\Delta t');
%set(gca, 'fontsize', 16);
