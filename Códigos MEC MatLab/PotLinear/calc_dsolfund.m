function  [dTdx,dTdy,dqdx,dqdy]=calc_dsolfund(x,y,xd,yd,nx,ny,k) 
%Calcula as solu��es fundamentais

r=sqrt((x-xd)^2+(y-yd)^2); % Raio (dist�ncia entre ponto fonte e 
                           % ponto campo)
rx=(x-xd); % Componente x do raio
ry=(y-yd); % Componente y do raio

dTdx=1/(2*pi*k*r^2)*rx; % Solu��o fundamental da temperatura
dTdy=1/(2*pi*k*r^2)*ry; % Solu��o fundamental da temperatura

dqdx=(nx*(rx^2 - ry^2) + 2*ny*rx*ry)/(2*pi*r^4);
dqdy=(ny*(-rx^2 +ry^2) + 2*nx*rx*ry)/(2*pi*r^4);