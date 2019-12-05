function [G,H]=mecmatrizvetor(y,dnorm,node,kmat)

pi2 = 2*pi;
n=size(y,2); % N�mero de n�s = N�mero de elementos

G=zeros(n,n); % Matriz A
H=zeros(n,n); % Matriz A
npg=30;
[xi,w]=Gauss_Legendre(-1,1,npg);

% x(2 x n)  = coordenadas dos extremos dos elementos
% (in�cio e final dos elementos). x=[x1 do n� 1, x1 do n� 2, ...
%       x2 do n� 1, x2 do n� 2, ...
% y(2 x n) = coordenadas dos n�s (pontos fontes)
% bc(2,n) = matriz de condi��es de contorno
% bc(1,i) = 0 => o deslocamento na dire��o x do n� i � conhecido
% bc(2,i) = 0 => o deslocamento na dire��o y do n� i � conhecido
% bc(1,i) = 1 => a for�a de superf�cie na dire��o x do n� i � conhecida
% bc(2,i) = 1 => a for�a de superf�cie na dire��o y do n� i � conhecida
for i=1:n % Loop on source points (Row)
    no1=node(1,i);
    no2=node(2,i);
    xx1=y(1,no1);
    yy1=y(2,no1);
    xx2=y(1,no2);
    yy2=y(2,no2);
    dgamadqsi=sqrt((xx2-xx1)^2+(yy2-yy1)^2)/2;
    for j=1:n % Loop on elements (Column)
        al = sqrt((y(1,node(2,j))-y(1,node(1,j)))^2 ...
            +(y(2,node(2,j))-y(2,node(1,j)))^2); % Comprimento do elemento
        for kk=1:npg % Laço sobre os pontos de integração
            N1=1/2*(1-xi(kk));% Função de forma N1 => linear contínua
            N2=1/2*(1+xi(kk));% Função de forma N2 => linear contínua
            x1=N1*xx1+N2*xx2; % Calcula a coordenada x do ponto de integração
            x2=N1*yy1+N2*yy2; % Calcula a coordenada y do ponto de integração
            % Compute parameters used in the formulas for the two intergals
            x11 = y(1,node(1,j))-x1; % diferen�a de coordenada x do
            % in�cio do elemento at� o ponto fonte
            x21 = y(2,node(1,j))-x2;% diferen�a de coordenada y do
            % in�cio do elemento at� o ponto fonte
            x12 = y(1,node(2,j))-x1;% diferen�a de coordenada x do
            % final do elemento at� o ponto fonte
            x22 = y(2,node(2,j))-x2;% diferen�a de coordenada y do
            % final do elemento at� o ponto fonte
            r1 = sqrt(x11^2 + x21^2); % Dist�ncia do ponto fonte ao in�cio
            % do elemento
            r2 = sqrt(x12^2 + x22^2); % Dist�ncia do ponto fonte ao final
            % do elemento
            % dnorm = vetor normal n1 e n2
            d = x11*dnorm(1,j) + x21*dnorm(2,j); % = delta x1 n1+delta y1 n2
            % = dr1/dn = h
            t1 = -x11*dnorm(2,j) + x21*dnorm(1,j); %=-delta x1 n2+delta y1 n1
            % = dr1/dt
            t2 = -x12*dnorm(2,j) + x22*dnorm(1,j); %=-delta x2 n2+delta y2 n1
            % = dr1/dt
            ds = abs(d);
            dtheta = atan2(ds*al,ds^2+t1*t2);
            % Elementos da matriz G
            g = (-dtheta*ds + al + t1*log(r1)-t2*log(r2))/(pi2*kmat);
            if(d<0 )
                dtheta = -dtheta;
            end
            % Elementos da matriz H
            h = dtheta/pi2;
            if(i==j)
                h = -0.5; % Diagonal da matriz H
            end
            H(i,j)=H(i,j)+h*dgamadqsi*w(kk); % Integral da matriz H
            G(i,j)=G(i,j)+g*dgamadqsi*w(kk); % Integral da matriz G
        end
    end
end
return