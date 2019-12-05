function g=calc_g_sing(x1,y1,x2,y2,k)

% Cálculo do comprimento do elemento
L = sqrt((x1-x2)^2+(y1-y2)^2);
g = L/(4*pi*k)*(3/2-log(L));
return

% Nesse ponto estão calculados os componentes singulares de [G]
% para os dados de entrada