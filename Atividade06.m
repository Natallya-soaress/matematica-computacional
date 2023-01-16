function Atividade06()
  
  clc;
  limiteInferior = 0;
  limiteSuperior = (pi/3);
  n = 4;
 
  dx = simpson(n, limiteInferior, limiteSuperior);
  dx2 = quad("f", limiteInferior, limiteSuperior);
  fprintf("A derivada exata é: %.6f \n", dx2);
  fprintf("A derivada calculada usando a regra 1/3 de Simpson é: %.6f \n", dx);
  
  tol = 1e-6;
  erro = quad("f", limiteInferior, limiteSuperior);
  d = 0;
  i = 1;
  
  while(abs(erro - d) > tol)
 
   d = simpson(i, limiteInferior, limiteSuperior);
   segmentos = i;
   i += 1;
   
endwhile
    fprintf("A quantidade de segmentos necessários para atingir a tolerância de 6 casas decimais é: %.0f \n", segmentos);

endfunction

function y = f(x)
  y = sin(x);
endfunction

function [dx] = simpson(n, limiteInferior, limiteSuperior)
  
  S1 = 0;
  S2 = 0;

  x = linspace(limiteInferior, limiteSuperior, n+1);
  x = x(2:end);
  
  for i=1:2:n-1
    S1 += f(x(i));
  endfor
  
  for i=2:2:n-2
    S2 += f(x(i));
  endfor
  
  dx = (limiteSuperior - limiteInferior) * (((f(limiteInferior)) + (4 * S1) + (2 * S2) + f(limiteSuperior))/(3 * n));
  
endfunction
