% Entrada de dados para análise de temperatura pelo
% método dos elementos de contorno


% Matriz para definição de pontos que definem a geometria
% PONTOS = [número do ponto, coord. x do ponto, coord. y do ponto]

PONTOS  = [1 0 0 ;
    2 1 0 ;
    3 1 1 ;
    4 0 1 ];

% Segmentos que definem a geometria
%  SEGMENTOS=[No do segmento, No do ponto inicial, No do ponto final,
%                                                  Raio, tipo do elemento]
% Raio do segmento: > 0 -> O centro à esquerda do segmento (do ponto
%                          inicial para o ponto final)
%                   < 0 -> O centro à direita do segmento (do ponto
%                          inicial para o ponto final)
%                   = 0 -> O segmento é uma linha reta


SEGMENTOS = [1 1 2 0;
    2 2 3 0;
    3 3 4 0;
    4 4 1 0];

% Matriz para definição da malha

% MALHA =[numero do segmento, numero de elementos no segmento]

ne=11;
MALHA = [1 ne;
    2 ne;
    3 ne;
    4 ne];

% Condições de contorno nos segmentos
% CCSeg=[no do segmento, tipo da CDC, valor da CDC]
% tipo da CDC = 0 => a temperatura é conhecida
% tipo da CDC = 1 => O fluxo é conhecido
CCSeg=[1 1 1
    2 0 1
    3 1 0
    4 0 0];

% Condutividade Térmica do material
k=1;

npi=21;
NPX=npi;
NPY=npi;