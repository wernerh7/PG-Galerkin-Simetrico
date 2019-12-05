clear; % Limpa todas as vari�veis
clc; % Apaga todo o texto das linhas de comando
close all % Fecha todas as figuras abertas

% Arquivos de entrada de dados
dad_2a; % Arquivos de dados de problemas gerais (dad_1, dad_2, ..., dad_6)

%%
disp('1. Formatando os dados');
% formata os dados
[XY,NOS,ELEM,CDC,normal,conect]=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg);

% mostra os elementos
mostra_elementos(SEGMENTOS,PONTOS,ELEM,XY)

% mostra as condi��es de contorno
mostra_cdc(SEGMENTOS,PONTOS,XY,NOS,CDC,ELEM,CCSeg) 

disp('2. Gerando pontos internos');
PONTOS_INT=gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); % gera os pontos internos 

disp('3. Construindo as matrizes.');

% mecmatrizvetor: baseada na function de mesmo nome do programa
% PotConstAnalitico (elementos de contorno de colocação com elementos
% constantes e integração analítica)

% calc_GeH: baseada na function de mesmo nome do programa
% PotLinearAnalitico (elementos de contorno de colocação com elementos
% constantes e integração analítica)

% calc_GeH: baseada na function de mesmo nome do programa
% PotLinearAnalitico (elementos de contorno de colocação com elementos
% constantes e integração analítica)

% hyp, hyp_c, hyp_a, hyp_ar: functions do livro do Glaucio Paulino
% (Symmetric Galerkin Boundary Element Method) para o cálculo das matrizes
% da equação hipersingular considerando elementos lineares contínuos. Veja
% o link: https://paulino.ce.gatech.edu/SGBEM_book/index.html.

[G1,~]=mecmatrizvetor(XY(:,2:3)',normal,ELEM(:,2:3)',1);
[H1,G2,~] = calc_GeH(XY(:,2:3)',ELEM(:,2:3)',1);
[H2,g] = hyp(XY(:,2),XY(:,3),ELEM(:,2:3)',conect);
[H2,g] = hyp_c(H2,g,XY(:,2),XY(:,3),ELEM(:,2:3)');
[H2,g] = hyp_a(H2,g,XY(:,2),XY(:,3),ELEM(:,2:3)',conect);
[H2,g] = hyp_ar(H2,g,XY(:,2),XY(:,3),ELEM(:,2:3)',conect);

disp('4. Aplicando as condições de contorno');

T_PR=calc_T_PR(CDC(:,2:3),ELEM(:,2:3)');

% Aplica as condicoes de contorno
 [A11,A12,A21,A22,b1,b2]= aplica_CDC(G1,H1,G2,H2,CDC(:,2:3),T_PR); 
 A=[A11, A12;-A21 ,-A22];
 b=[b1' ;-b2'];
 
 disp('5. Resolvendo o sistema linear');

 x=A\b;
%  Separa temperatura de fluxo
[T,q]= monta_Teq(CDC,T_PR,x);


disp('6. Calculando a temperatura nos pontos internos');
T_pint = calc_T_pint(PONTOS_INT,XY,ELEM,k,T,q);

disp('7. Gerando o mapa de cor');

mapa_de_cor(T,XY,PONTOS_INT,PONTOS,SEGMENTOS,T_pint);

disp('8. Verificando a simetria da matriz A: ')
disp(' ')
disp(['   max(max(abs(A-A^T))): ',num2str(max(max(abs(A-A'))))])
                            
