function Atividade03()
  
  clc;
  Aa = [3, 5, -5, 21;
        -4, 8, -5, 1;
        2, -5, 6, -16];
        
  Aa = pivotamentoParcial(Aa);
  [L, U, X] = fatLU(Aa);
  
  fprintf('Solução do sistema: \n');
  fprintf('X= %.1f \n', X(1));
  fprintf('Y= %.1f \n', X(2));
  fprintf('Z= %.1f \n', X(3));
  fprintf('\n');
  fprintf('Matriz L: \n\n');
  disp(L);
  fprintf('\n');
  fprintf('Matriz U: \n\n');
  disp(U);
  
endfunction

function A = pivotamentoParcial(Aa);
  
  dimensao = size(Aa);
  
  for i = 1:dimensao(1)-1
    [~, pos] = max(abs(Aa(:, i)));
    if pos + i - 1 != i
      aux = Aa(i, :);
      Aa(i, :) = Aa(pos + i - 1, :);
      Aa(pos + i - 1, :) = aux;
    endif
  endfor
  A = Aa;
endfunction

function [L, U, X] = fatLU(Aa)

  dimensao = size(Aa);
  A = Aa(1:3, 1:3); 
  L = eye(size(A));
  U = A;
  b = Aa(1:3, 4);
  
  for i =1:dimensao(1)
    for j=i+1:dimensao(1)
      fator = U(j, i)/U(i, i);
      L(j, i) = fator;
      U(j, :) = U(j, :) - fator * U(i, :);
      endfor
  endfor
  
  Lb = L;
  Lb(:, 4) = b;
  dimensaoLb = size(Lb);
  aux = 0;

  for i= 1: dimensaoLb(1)
    d(i) = Lb(i, dimensaoLb(2)) - sum(Lb(i, 1:i-1))/ Lb(i, i); 
    Lb(:, i) = Lb(:, i) * d(i);
  endfor
  
  Ud = U;
  Ud(:, 4) = d;
  dimensaoUd = size(Ud);
  i = 3;
  aux = 0;
  
  while i > 0
    X(i) = (Ud(i,dimensaoUd(2)) - sum(Ud(i, i+1:dimensaoUd(2)-1)))/Ud(i, i);
    Ud(:, i) = Ud(:, i) *  X(i);
    i -= 1;
  endwhile
endfunction


























