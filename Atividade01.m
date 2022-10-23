function Atividade01()
  clc
  xl = 1;
  xu = 2;
  [iter, xr, XR, FXR] = bissecao(xl, xu);

  fprintf('Raiz: %4.6f \n', xr)
  fprintf('Número de iterações: %i \n', iter)

  graficoAnimado(xl, xu, XR, FXR, iter);
  graficoConvergencia(iter, XR, FXR);
endfunction

function y = f(x)
  y = ((x^4) - (8 * (x^3)) + (23 * (x^2)) + (16 * x) - 50);
endfunction

function [iter, xr, XR, FXR] = bissecao(xl, xu)

  tol = 1e-5;
  imax = 1000;
  XR = zeros(imax, 1);
  FXR = zeros(imax, 1);
  xrold = inf;

  for iter = 1:imax
    xr = (xl+xu)/2;
    XR(iter) = xr;
    FXR(iter) = f(xr);

    if abs(xr - xrold) < tol
      break
    endif

    xrold = xr;

    if f(xr) * f(xl) < 0
      xu = xr;
    elseif f(xr) * f(xl) > 0
      xl = xr;
    endif
  endfor

  XR(iter+1:end) = [];
  FXR(iter+1:end) = [];

endfunction

function graficoAnimado(xl, xu, XR, FXR, imax)

  vetorX = xl:0.1:xu;
  vetorY = zeros(size(vetorX));

  for cont = 1:numel(vetorX)
    vetorY(cont) = f(vetorX(cont));
  endfor

  for iter = 1:imax
    figure(1);
    plot(vetorX, vetorY, 'linewidth', 2, 'color', 'black');
    hold on;
    plot(XR(iter), FXR(iter), 'markersize', 50, 'color', 'magenta');
    hold off;
    set(gca, 'fontsize', 20);
    grid on;
    title(sprintf('Iteração: %i, f(%.6f)= %.6f', iter, XR(iter), FXR(iter)));
    pause(0.1);
  endfor
endfunction

function graficoConvergencia(iter, XR, FXR)
  figure(2);
  subplot(2, 1, 1);
  plot(1:iter, XR, 'linewidth', 2, 'color', 'magenta');
  xlim([1 iter]);
  ylim([min(XR) max(XR)])
  set(gca, 'fontsize', 20);
  grid on;
  xlabel('Iteração');
  ylabel('x');
  title(sprintf('Convergência de x'));

  subplot(2, 1, 2);
  plot(1:iter, FXR, 'linewidth', 2, 'color', 'magenta');
  xlim([1 iter]);
  ylim([min(FXR) max(FXR)])
  set(gca, 'fontsize', 20);
  grid on;
  xlabel('Iteração');
  ylabel('f(x)');
  title(sprintf('Convergência de f(x)'));
endfunction
