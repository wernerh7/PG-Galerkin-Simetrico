function mapa_de_cor(temp,dTdx,dTdy,NOS,ELEM,PONTOS_INT,PONTOS, ...
    SEGMENTOS,temp_int)
% Rotina para criação do mapa de cores para temperatura
%
%   Autor: Frederico Lourenço
%   Data de criação julho de 1999
%   Revisão 0.0

% Atribuição das coordenadas e temperatura nos nós e pontos
% internos à matriz Tf[x,y,T]
n_nos=length(ELEM(:,1));
npi=length(PONTOS_INT(:,1));
n_linhas = length(SEGMENTOS(:,1));
xmin = min(PONTOS(:,2));
xmax = max(PONTOS(:,2));
ymin = min(PONTOS(:,3));
ymax = max(PONTOS(:,3));
lx = xmax - xmin;		% Largura do retângulo que contém a geometria
ly = ymax - ymin;		% Altura do retângulo que contém a geometria
lmax=sqrt(lx^2+ly^2);

Tf = zeros(n_nos+npi,3);
Tf(1:n_nos,:) = [NOS(:,2),NOS(:,3), temp'];
Tf(n_nos+1:n_nos+npi,:) = [PONTOS_INT(:,2),PONTOS_INT(:,3), temp_int'];


% Criação dos vetores dos valores de x, de y e de T
xi = Tf(:,1);
yi = Tf(:,2);
ci = Tf(:,3);

% Chamada da função de triangularização (geração de malha
% na domínio para mostrar o mapa de cores)
tri = delaunay(xi,yi);

% Geração de cada triângulo da malha com interpolação bilinear
% dos resultados obtidos
for t = 1 : n_linhas		% Percorre todas as linhas
  xl1(t) = PONTOS(SEGMENTOS(t,2),2);
  xl2(t) = PONTOS(SEGMENTOS(t,3),2);
  yl1(t) = PONTOS(SEGMENTOS(t,2),3);
  yl2(t) = PONTOS(SEGMENTOS(t,3),3);
  raio(t)= SEGMENTOS(t,4);  
end;


for t = 1 : length(tri)	% Percorre todos os triângulos da malha
    
    x = [xi(tri(t,1));xi(tri(t,2));xi(tri(t,3))];
    y = [yi(tri(t,1));yi(tri(t,2));yi(tri(t,3))];
    
    % Verificando se o triângulo pertence ao domínio
    xm = (x(1)+x(2)+x(3))/3;
    ym = (y(1)+y(2)+y(3))/3;
    [xw,ponto] = testa_ponto(xm,ym,xl1,yl1,xl2,yl2,lx,raio);
    
    if strcmp(ponto,'interno')
        
        x = [xi(tri(t,1));xi(tri(t,2));xi(tri(t,3))];
        y = [yi(tri(t,1));yi(tri(t,2));yi(tri(t,3))];
        c = [ci(tri(t,1));ci(tri(t,2));ci(tri(t,3))];
        patch(x,y,c)
    end;
    
end;
hold on
view(0,90)
colorbar 			% Cria a barra de cores
hold on
view(0,90)
colorbar 			% Cria a barra de cores
title('Distribuição de temperatura');
shading interp
dTdxmax=max(dTdx);
dTdymax=max(dTdy);
dTmax=sqrt(dTdxmax^2+dTdymax^2);
escala=0.1*lmax/dTmax;
quiver(PONTOS_INT(:,2),PONTOS_INT(:,3),escala*dTdx',escala*dTdy')




