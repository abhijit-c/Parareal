function parareal_stability(n, k)
  N = 100;
  [MDT, mdt] = meshgrid(linspace(-4, 0, N), linspace(-2, 2, N));

  fe = @(z) (1 + z);
  be = @(z) (1 + z)^(-1);

  s = MDT./mdt;
  roots = zeros(N,N);
  for j = 0:k
    roots = roots + nchoosek(n,j) .* ...
            (fe(mdt).^s-be(MDT)).^j .* ...
            be(MDT).^(n-j) ;
  end

  roots = abs(roots);
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
  surf(MDT, mdt, root_cond); colormap('gray');
  view(0, 90);
  title('Absolute Stability region of Parareal.');
  xlabel('\mu \Delta t'); ylabel('\mu \delta t');
  %set(gca, 'fontsize', 16);
end
