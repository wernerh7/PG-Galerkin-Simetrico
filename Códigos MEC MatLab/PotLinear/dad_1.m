% Entrada de dados  

% Matriz para defini��o de pontos 

PONTOS  = [1   0  0 
    2 .5 -.5
          3   1  0
          4 .5 .5];
          
% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio] 
% Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
%                          inicial para o ponto final) 
%                   < 0 -> O centro � � direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento � uma linha reta

SEGMENTOS = [1 1 2  .5
     2 2 3  .5
      3 3 4  .5
       4 4 1  .5];

% Matriz para defini��o da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

MALHA = [1  6
          2 6
          3 6
          4 6];
      
 k=1; % Condutividade t�rmica do material
 
   % Condi��es de contorno nos segmentos
  % CCSeg=[no do segmento, tipo da CDC, valor da CDC]
  % tipo da CDC = 0 => a temperatura � conhecida
  % tipo da CDC = 1 => O fluxo � conhecido
 CCSeg=[1 0 0
     2 1 0
     3 1 0
     4 1 0];
 
% Defini��o dos pontos internos 
NPX=9; % N�mero de pontos internos na dire��o X
NPY=9; % N�mero de pontos internos na dire��o Y

% fc = fonte de calor concentrada
% fc = [valor da fonte, coordenada x da fonte, coordenada y da fonte];
fc=[0 .0 .0]; 
