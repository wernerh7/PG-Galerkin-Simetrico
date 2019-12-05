% Entrada de dados para an�lise de temperatura pelo
% m�todo dos elementos de contorno

% Matriz para defini��o de pontos que definem a geometria
% PONTO = [n�mero do ponto, coord. x do ponto, coord. y do ponto]

PONTOS  = [1   0   0 ;
    2   1   0 ;
    3   2   1 ;
    4   2   2 ;
    5   1   3 ;
    6   0   3];

% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio, tipo do elemento]
% Raio do segmento: > 0 -> O centro � � esquerda do segmento (do ponto
%                          inicial para o ponto final)
%                   < 0 -> O centro � � direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento � uma linha reta

SEGMENTOS = [1 1 2 0;
    2 2 3 0;
    3 3 4 0;
    4 4 5 0;
    5 5 6 0;
    6 6 1 0];

% Matriz para defini��o da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

MALHA = [1  12 ;
    2  12 ;
    3  12;
    4  12;
    5 12 ;
    6  24];

% Condi��es de contorno nos segmentos
% CCSeg=[no do segmento, tipo da CDC, valor da CDC]
% tipo da CDC = 0 => a temperatura � conhecida
% tipo da CDC = 1 => O fluxo � conhecido
CCSeg=[1 0 7
    2 1 0
    3 0 5
    4 1 3
    5 0 6
    6 1 -2];

% Condutividade T�rmica do material
k=1;

NPX=21;
NPY=21;