function q=calc_q(xd,yd,fc,k)

valor_fonte=fc(1);
x_fonte=fc(2);
y_fonte=fc(3);

% Calcula as soluções fundamentais

R=sqrt((x_fonte-xd)^2+(y_fonte-yd)^2); % Raio (distância entre ponto fonte e 
                           % ponto campo)
Tast=-1/(2*pi*k)*log(R); % Solução fundamental da temperatura

q=valor_fonte*Tast;
