function q=calc_q(xd,yd,fc,k)

valor_fonte=fc(1);
x_fonte=fc(2);
y_fonte=fc(3);

% Calcula as solu��es fundamentais

R=sqrt((x_fonte-xd)^2+(y_fonte-yd)^2); % Raio (dist�ncia entre ponto fonte e 
                           % ponto campo)
Tast=-1/(2*pi*k)*log(R); % Solu��o fundamental da temperatura

q=valor_fonte*Tast;
