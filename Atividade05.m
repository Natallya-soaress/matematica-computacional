function Atividade05()
  clc;
  X = [0.1; 0.2; 0.4; 0.6; 0.9; 1.3; 1.5; 1.7; 1.8];
  Y = [0.75; 1.25; 1.45; 1.25; 0.85; 0.55; 0.35; 0.28; 0.18];
  
  XX = 1.4;
  ordem = 8;
  if(ordem >= numel(X))
    fprintf("Não é possível calcular o polinômio!!\n");
  else
    a = regressaoPolinomial(X, Y, ordem);
    printaPolinomio(a);
    grafico(X, Y, XX, a);
    [r, r2] = calculaR(X, Y, a);
    fprintf("O coeficiente de determinação (r2) é: %f\nO coeficiente de correlação (r) é: %f\n", r2, r);
  endif
  
endfunction

function [a] = regressaoPolinomial(X, Y, ordem)
  
  Aa = zeros(ordem+1, ordem+2);
  
  for i=1:ordem+1
    for j=1:i
      k = i + j - 2;
      soma = 0;
      for l=1:numel(X)
        soma = soma + X(l)^k;
      endfor
        Aa(i, j) = soma;
        Aa(j, i) = soma;
      endfor
      soma = 0;
      
      for l=1:numel(X)
        soma = soma + Y(l) * X(l)^(i-1);
      endfor
      Aa(i, ordem+2) = soma;
  endfor
  
  A = Aa(:, 1:end-1);
  B = Aa(:, end);
  a = A\B
endfunction

function [pol] = calculaFuncao(a, x)
  pol = 0;
  for i=numel(a):-1:1
    pol = pol + (a(i) * x^(i-1));
  endfor
endfunction

function printaPolinomio(a)
  pol = "";
  i = numel(a);
  while(i > 1)
    pol = [pol, num2str(a(i)), " * x^", num2str(i-1), " + "]; 
    i--;
  endwhile
  pol = [pol, num2str(a(i)), " * x"]; 
  fprintf("O polinômio fitado é: %s \n", pol);
endfunction

function grafico(X, Y, XX, a)
  ponto = [XX, calculaFuncao(a, XX)];
  i = 1;
  
  for j = min(X):0.1:max(X)
    allX(i) = j;
    allFX(i) = calculaFuncao(a, j);
    i = i + 1;
  endfor
  
  clf;
  figure(1);
  hold on;
  plot(X, Y, 'o', 'markersize', 10, 'color', 'magenta')
  plot(allX, allFX, 'linewidth', 2, 'color', 'black')
  plot(ponto(1), ponto(2), '.', 'color', 'red', 'markersize', 30)
  hold off;
  set(gca,'fontsize', 20);
  grid on;
  xlabel('X');
  ylabel('F(x)');
  title(sprintf('Convergencia de F(X)'));
endfunction

function [r, r2] = calculaR(X, Y, a)
  St = 0;
  Sr = 0;
  for i = 1:numel(X)
    St = St + (Y(i) - (sum(Y)/numel(Y)))^2;
    Sr = Sr + (Y(i) - calculaFuncao(a, X(i)))^2;
  endfor
  r2 = (St-Sr)/(St);
  r = sqrt(r2);
endfunction
  







