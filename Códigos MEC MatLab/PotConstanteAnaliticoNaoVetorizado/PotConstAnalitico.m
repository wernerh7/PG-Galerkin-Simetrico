clear; % Limpa todas as variáveis
clc; % Apaga todo o texto das linhas de comando
close all % Fecha todas as figuras abertas

% Arquivos de entrada de dados
dad_2a; % Arquivos de dados de problemas gerais (dad_1, dad_2, ..., dad_6)

%%
disp('1. Formatando os dados');
% formata os dados
[XY,NOS,ELEM,CDC,normal]=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg);

% mostra os elementos
mostra_elementos(SEGMENTOS,PONTOS,ELEM,XY)

% mostra as condições de contorno
mostra_cdc(SEGMENTOS,PONTOS,XY,NOS,CDC,ELEM,CCSeg) 

disp('2. Gerando pontos internos');
PONTOS_INT=gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); % gera os pontos internos 

%% Elementos de contorno padrão
disp('Resolvendo o elementos de contorno')
% Calcula a matriz A e o vetor b
%[A,b]=mecmatrizvetor(NOS(:,2:3)',XY(:,2:3)',CDC(:,2:3)',normal, ...
 %   ELEM(:,2:3)',k); % usa a equação tradicional
[A,b]=mecmatrizvetor_traction(NOS(:,2:3)',XY(:,2:3)',CDC(:,2:3)', ...
     normal,ELEM(:,2:3)',k); % usa a equação hipersingular
x =A\b';
% x2 =A2\b2';

% Separa temperatura de fluxo
[T,q]= monta_Teq(CDC,x);

disp('4. Calculando temperatura e fluxo nos pontos internos');

% Calcula a temperatura e derivada de temperatura nos pontos internos
[Ti,dTdx,dTdy]=domain_field(PONTOS_INT(:,2:3)',XY(:,2:3)',CDC(:,2:3)', ...
    ELEM(:,2:3)',normal,x,k);

disp('5. Gerando o mapa de cor');

mapa_de_cor(T,dTdx,dTdy,NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,Ti) 