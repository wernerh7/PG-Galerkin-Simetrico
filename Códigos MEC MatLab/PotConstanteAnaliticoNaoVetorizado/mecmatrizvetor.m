function [A,b]=mecmatrizvetor(x,y,bc,dnorm,node,kmat)

pi2 = 2*pi; 
n=size(x,2); % Número de nós = Número de elementos

A=zeros(n,n); % Matriz A
b=zeros(1,n); % Vetor b
% x(2 x n)  = coordenadas dos extremos dos elementos 
% (início e final dos elementos). x=[x1 do nó 1, x1 do nï¿½ 2, ...
                             %       x2 do nó 1, x2 do nï¿½ 2, ...
% y(2 x n) = coordenadas dos nós (pontos fontes)
% bc(2,n) = matriz de condições de contorno
% bc(1,i) = 0 => o deslocamento na direção x do nó i é conhecido
% bc(2,i) = 0 => o deslocamento na direção y do nó i é conhecido
% bc(1,i) = 1 => a força de superfície na direção x do nó i é conhecida
% bc(2,i) = 1 => a força de superfície na direção y do nó i é conhecida

for j=1:n % Loop on elements (Column)
    al = sqrt((y(1,node(2,j))-y(1,node(1,j)))^2 ...
        +(y(2,node(2,j))-y(2,node(1,j)))^2); % Comprimento do elemento
    for i=1:n % Loop on source points (Row)
        % Compute parameters used in the formulas for the two intergals
        x11 = y(1,node(1,j))-x(1,i); % diferença de coordenada x do 
        % início do elemento até o ponto fonte
        x21 = y(2,node(1,j))-x(2,i);% diferença de coordenada y do 
        % início do elemento até o ponto fonte
        x12 = y(1,node(2,j))-x(1,i);% diferença de coordenada x do 
        % final do elemento até o ponto fonte
        x22 = y(2,node(2,j))-x(2,i);% diferença de coordenada y do 
        % final do elemento até o ponto fonte
        r1 = sqrt(x11^2 + x21^2); % Distância do ponto fonte ao início 
                       % do elemento
        r2 = sqrt(x12^2 + x22^2); % Distância do ponto fonte ao final 
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
        if(bc(1,j)==0) % A temperatura é conhecida
            A(i,j) = A(i,j)-g; % Os elementos de G vão para a matriz A
            b(i) = b(i) - h*bc(2,j);% Os elementos de H vão para o vetor b
        else % O fluxo é conhecido
            A(i,j) = A(i,j) + h; % Os elementos de H vão para a matriz A
            b(i) = b(i) + g*bc(2,j);% Os elementos de G vão para o vetor b
        end
    end
end
return