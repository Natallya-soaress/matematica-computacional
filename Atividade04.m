function Atividade04()

  clc;
  xrold = [2.5; 2.5];
  [iter, xr, XX, FXX] = metodoBFGS(xrold);
  
  flag = 1; % Alterar a flag para mudar o método chamado (1 = gradiente, 2 = Newton, 3 = BFGS)
  
  if(flag == 1)
    [iter, xr, XX, FXX] = metodoGradiente(xrold);
    salto = 0.1;
    fprintf('Metodo do gradiente: \n');
  elseif (flag == 2)
    [iter, xr, XX, FXX] = metodoNewton(xrold);
    salto = 0.1;
    fprintf('Metodo de Newton: \n');
  elseif(flag == 3)
    [iter, xr, XX, FXX] = metodoBFGS(xrold);
    salto = 1000;
    fprintf('Metodo de BFGS: \n');
  endif

  fprintf('X1= %.5f \n', xr(1));
  fprintf('X2= %.5f \n', xr(2));
  fprintf('Número de iterações: %i \n', iter);
  graficoConvergencia(iter, XX, FXX);
  graficoAnimado(iter, XX, FXX, salto);
  
endfunction

function y = f(x)
  y = (((x(1)^2) + x(2) - 11)^2) + ((x(1) + (x(2)^2) - 7)^2);
endfunction

function y = gradiente(x)
  y = [((4 * (x(1)^3)) - (42 * x(1)) + (4 * x(1) * x(2)) + (2 * (x(2)^2)) - 14); ((2 * x(1)^2) - 22 + (4 * x(1) * x(2)) + (4 * (x(2)^3)) - (26 * x(2)))];
endfunction

function y = hessiana(x)
  y = [((12 * (x(1)^2)) - 42  + (4 * x(2))), ((4 * x(1)) + (4 * x(2))); ((4 * x(1)) + (4 * x(2))), ((4 * x(1)) + (12 * (x(2)^2)) - 26)];
endfunction

function [iter, xr, XX, FXX] = metodoGradiente(xrold)

  tol = 1e-5;
  imax = 100;
  alfa = 0.01;
  XX = zeros(imax, 2);
  FXX = zeros(imax, 2);

  for iter = 1:imax
    xr = xrold - (alfa * gradiente(xrold));
    XX(iter, :) = xr;
    FXX(iter, :) = f(xr);

   if max(abs(xr - xrold)) <= tol
      break
   endif
   xrold = xr;
  endfor

  XX(iter+1:end, :) = [];
  FXX(iter+1:end, :) = [];
endfunction

function [iter, xr, XX, FXX] = metodoNewton(xrold)

  tol = 1e-5;
  imax = 100;
  alfa = 1;
  XX = zeros(imax, 2);
  FXX = zeros(imax, 2);

  for iter = 1:imax
    xr = xrold - (alfa * inv(hessiana(xrold)) * gradiente(xrold));
    XX(iter, :) = xr;
    FXX(iter, :) = f(xr);

   if max(abs(xr - xrold)) <= tol
      break
   endif
   xrold = xr;
  endfor

  XX(iter+1:end, :) = [];
  FXX(iter+1:end, :) = [];
endfunction

function [iter, xr, XX, FXX] = metodoBFGS(xrold)

  tol = 1e-5;
  imax = 100;
  alfa = 1.5;
  D = eye(2);

  XX = zeros(imax, 2);
  FXX = zeros(imax, 2);

  for iter = 1:imax
    
    S = - alfa * D * gradiente(xrold);
    Y = gradiente(xrold + S) - gradiente(xrold);
    
    xr = xrold - (D * gradiente(xrold));
    
    D = D + (((((transpose(S) * Y)  + (transpose(Y) * D * Y)) * (S * transpose(S)))/((transpose(S) * Y)^2)) - (((D * Y * transpose(S)) + (S * transpose(Y) * D))/(transpose(S) * Y)));
    
    XX(iter, :) = xr;
    FXX(iter, :) = f(xr);

   if max(abs(xr - xrold)) <= tol
      break
   endif
   xrold = xr;
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
  legend('X1', 'X2');

  subplot(2, 1, 2)
  plot(1:iter, FXX, 'linewidth', 2);
  set(gca, 'fontsize', 20);
  grid on;
  xlabel('Iteração');
  ylabel('F(x)');
  title(sprintf('Convergência de F(X)'));
  legend('F(X1)', 'F(X2)');
endfunction

function graficoAnimado(iter, todosX, todosY, salto)
  
  xl = min(min(todosX));
  xu = max(max(todosX));
  
  x = xl:salto:xu;
  y = x;
  F = zeros(length(x), length(y), 1);
  
  for xx = 1:numel(x)
    for yy = 1:numel(y)
      F(xx, yy, 1) = f([x(xx); y(yy)]);
    end
  end  
  
  for cont = 1:iter
    figure(2);
    surf(x, y, F(:, :, 1));
    hold on;
    plot3(todosX(1:cont, 2), todosX(1:cont, 1), todosY(1:cont, 1), 'color', 'r', 'linewidth', 2)
    plot3(todosX(cont, 2), todosX(cont, 1), todosY(cont, 1), '*', 'markersize', 20, 'color', 'r', 'markerfacecolor', [1 1 1], 'linewidth', 2)
    hold off;
    set(gca, 'fontsize', 10);
    grid on;
    title(sprintf('Iteração: %i, f(%.6f, %.6f)= %.6f', cont, todosX(cont, 1), todosX(cont, 2), todosY(cont)));
    pause(0.1);
  endfor
endfunction
