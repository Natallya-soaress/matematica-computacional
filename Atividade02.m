function Atividade02()

  clc
  xr1old = [1; 2; 3];
  [iter, xr, XX, FXX] = newtonRaphson(xr1old);

  fprintf('X1= %.5f \n', xr(1));
  fprintf('X2= %.5f \n', xr(2));
  fprintf('X3= %.5f \n', xr(3));
  fprintf('Número de iterações: %i \n', iter);

  graficoConvergencia(iter, XX, FXX);
endfunction

function y = f(x)
  y = [(( 2 * x(1)) - x(2) - cos(x(1))); (-x(1)) + (2 * x(2)) - (x(3)) - (cos(x(2))); (-x(2)) + x(3) - (cos(x(3)))];
endfunction

function y = df(x)
  y = [2 + sin(x(1)), -1, 0; -1, 2+sin(x(2)), -1; 0, -1, 1 + sin(x(3))];
endfunction

function [iter, xr, XX, FXX] = newtonRaphson(xr1old)

  tol = 1e-5;
  imax = 1000;
  XX = zeros(imax, 3);
  FXX = zeros(imax, 3);

  for iter = 1:imax
    xr = xr1old - inv(df(xr1old)) * f(xr1old);
    XX(iter, :) = xr;
    FXX(iter, :) = f(xr);

   if max(abs(xr - xr1old)) <= tol
      break
   endif
   xr1old = xr;
  endfor

  XX(iter+1:end, :) = [];
  FXX(iter+1:end, :) = [];
endfunction

function graficoConvergencia(iter, XX, FXX)
  figure(1);
  subplot(2, 1, 1)
  plot(1:iter, XX, 'linewidth', 2);
  set(gca, 'fontsize', 20);
  grid on;
  xlabel('Iteração');
  ylabel('X');
  title(sprintf('Convergência de X'));
  legend('X1', 'X2', 'X3');

  subplot(2, 1, 2)
  plot(1:iter, FXX, 'linewidth', 2);
  set(gca, 'fontsize', 20);
  grid on;
  xlabel('Iteração');
  ylabel('F(x)');
  title(sprintf('Convergência de F(X)'));
  legend('F(X1)', 'F(X2)', 'F(X3)');
endfunction
