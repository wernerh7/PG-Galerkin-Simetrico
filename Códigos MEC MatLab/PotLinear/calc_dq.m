function [dqdx,dqdy]=calc_dq(xd,yd,fc,k)

valor_fonte=fc(1);
x_fonte=fc(2);
y_fonte=fc(3);

%Calcula as solu��es fundamentais


rx=(x_fonte-xd); % Componente x do raio
ry=(y_fonte-yd); % Componente y do raio
r=sqrt(rx^2+ry^2);
dTdx=1/(2*pi*k*r^2)*rx; % Solu��o fundamental da temperatura
dTdy=1/(2*pi*k*r^2)*ry; % Solu��o fundamental da temperatura

dqdx=valor_fonte*dTdx;
dqdy=valor_fonte*dTdy;