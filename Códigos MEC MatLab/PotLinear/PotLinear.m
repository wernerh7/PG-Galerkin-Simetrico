% Introdução ao Método dos Elementos de Contorno
% Programa de elementos de contorno aplicado a problemas de condução de
% calor bidimensional. Pode haver fontes de calor concentradas
% Tipo de elementos: lineares contínuos
% Por: Éder Lima de Albuquerque
% Brasília, janeiro de 2013

clc; % Limpa as linhas de comando
clear; % Apaga as variáveis da memória
close all % Fecha todas as janelas gráficas que estão abertas
dad_2a; % Arquivo de entrada de dados

[NOS,ELEM,CDC]=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg);% formata os dados

mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS) % mostra a geometria do problema
mostra_cdc(SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg) % mostra as condições
                                          % de contorno
disp('1. Gerando pontos internos');
PONTOS_INT=gera_p_in(NPX,NPY,PONTOS,SEGMENTOS); % gera os pontos internos 

disp('2. Construindo as matrizes H e G');
[G,H,g]=monta_GeH(ELEM,NOS,k,fc); % Monta as matrizes H e G

disp('3. Aplicando as condições de contorno');
[A,b,T_PR]=aplica_CDC(G,H,CDC,ELEM); % Aplica as condições de contorno

disp('4. Resolvendo o sistema linear');
% x=A\(b-g); % Calcula as variáveis desconhecidas

x=A\(b); % Calcula as variáveis desconhecidas

disp('5. Ordenando fluxo e temperatura');
[T,q,q_vet] = Monta_Teq(x,CDC,T_PR); % Separa temperatura e fluxo

disp('6. Calculando a temperatura nos pontos internos');
T_pint=calc_T_pint(ELEM,NOS,PONTOS_INT,T,q_vet,k,fc); % Calcula 
            % a temperatura nos pontos internos
[dTdx,dTdy]=calc_dT_pint(ELEM,NOS,PONTOS_INT,T,q_vet,k,fc); % Calcula as  
                                    %  derivadas da temperatura nos 
                                    %  pontos internos

disp('7. Criando o mapa de cor');
mapa_de_cor(T,dTdx,dTdy,NOS,ELEM,PONTOS_INT,PONTOS,SEGMENTOS,T_pint) % Faz 
                             % um mapa de cor
                            