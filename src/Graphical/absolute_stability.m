N = 100;
[A, B] = meshgrid(linspace(0, 1, N), linspace(0, 1, N));

fws = @(z) (1 + z); bws = @(z) (1-z); trap = @(z) (1+0.5*z)./(1-0.5*z);

roots = abs(H(fws, fws, 10, 1, -5, A, B));

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
surf(A, B, root_cond);
view(0, 90);
title('Absolute Stability region of Parareal.');
xlabel('\delta t'); ylabel('\Delta t');
set(gca, 'fontsize', 16);

function v = H(G, F, n, k, mu, dt, Dt)
  s = Dt/dt;
  v = 0.0;
  for j=0:k
    v = v + nchoosek(n,j) .* ...
            (F(mu*dt).^s - G(mu*Dt)).^j .* ...
            G(mu*Dt).^(n-j);
  end
end

