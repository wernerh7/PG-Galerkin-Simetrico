function mapa_de_cor(temp,NOS,PONTOS_INT,PONTOS,SEGMENTOS,T_in)
% Rotina para cria��o do mapa de cores para temperatura
%
%   Autor: Frederico Louren�o
%   Data de cria��o julho de 1999
%   Revis�o 0.0

% Atribui��o das coordenadas e temperatura nos n�s e pontos
% internos � matriz Tf[x,y,T]
n_linhas = length(SEGMENTOS(:,1));
xmin = min(PONTOS(:,2));
xmax = max(PONTOS(:,2));
ymax = max(PONTOS(:,3));
ymin = min(PONTOS(:,3));

lx = xmax - xmin;		% Largura do ret�ngulo que cont�m a geometria
ly = ymax - ymin;		% Altura do retï¿½ngulo que contï¿½m a geometria
lmax=sqrt(lx^2+ly^2);

X=[NOS(:,2:3) ; PONTOS_INT(:,2:3)];

% Chamada da funï¿½ï¿½o de triangularizaï¿½ï¿½o (geraï¿½ï¿½o de malha
% na domï¿½nio para mostrar o mapa de cores)
tri = delaunay(X);

% Gera��o de cada tri�ngulo da malha com interpola��o bilinear
% dos resultados obtidos
for t = 1 : n_linhas		% Percorre todas as linhas
  xl1(t) = PONTOS(SEGMENTOS(t,2),2);
  xl2(t) = PONTOS(SEGMENTOS(t,3),2);
  yl1(t) = PONTOS(SEGMENTOS(t,2),3);
  yl2(t) = PONTOS(SEGMENTOS(t,3),3);
  raio(t)= SEGMENTOS(t,4);  
end


% Gera��o de cada tri�ngulo da malha com interpola��o bilinear
% das temperaturas
icount=0;
for t = 1 : length(tri)	% Percorre todos os tri�ngulos da malha
    x = [X(tri(t,1),1);X(tri(t,2),1);X(tri(t,3),1)];
    y = [X(tri(t,1),2);X(tri(t,2),2);X(tri(t,3),2)];
    
  
  % Verificando se o tri�ngulo pertence ao dom�nio
  xm = (x(1)+x(2)+x(3))/3;
  ym = (y(1)+y(2)+y(3))/3;
  [xm,ponto] = testa_ponto(xm,ym,xl1,yl1,xl2,yl2,lx,raio);

    if strcmp(ponto,'interno')
        icount=icount+1;
        MALHA(icount,:)=tri(t,:);
    end

end
T=[temp';T_in];
patch('faces',MALHA,'vertices',X, ...
    'FaceVertexCData',T,'CDataMapping','scaled','FaceColor','interp');
hold on
view(0,90)
colorbar 			% Cria a barra de cores
title('Distribuicao de temperatura');
shading interp