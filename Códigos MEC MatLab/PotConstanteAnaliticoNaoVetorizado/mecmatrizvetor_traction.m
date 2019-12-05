function [A,b]=mecmatrizvetor_traction(x,y,bc,dnorm,node,kmat)

% -------------------------------------------------------------------------------
%  Fucntion used to evaluate the temperature/heat flux using the direct
%     evaluation of the integral representation (dfdx and dfdy not implemented)
%
% -------------------------------------------------------------------------------

pi2 = pi*2;
n=size(x,2); % Número de nós = Número de elementos

A=zeros(n,n); % Matriz A
b=zeros(1,n); % Vetor b

n=size(y,2);
for i=1:n
    mx=dnorm(1,i); % normal x no ponto fonte
    my=dnorm(2,i); % normal y no ponto fonte
    
    for j=1:n
        
        % Comprimento do elemento
        al= sqrt((y(1,node(2,j))-y(1,node(1,j)))^2+ ...
            (y(2,node(2,j))-y(2,node(1,j)))^2);
        nx=dnorm(1,j);
        ny=dnorm(2,j);
        tx=-ny;
        ty=nx;
        
        x11 = y(1,node(1,j)) - x(1,i);
        x21 = y(2,node(1,j)) - x(2,i);
        x12 = y(1,node(2,j)) - x(1,i);
        x22 = y(2,node(2,j)) - x(2,i);
        r1 =  sqrt(x11^2 + x21^2); % Distância para o início do elemento
        r2 =  sqrt(x12^2 + x22^2); % Distância para o final do elemento
        
        % Projeção do vetor distância no vetor normal ao elemento
        d  =  x11*dnorm(1,j) + x21*dnorm(2,j); % Figura A.1, página 178
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j); % Distância T1 da figura A.1
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j); % Distância T2 da figura A.1
        
        ds = abs(d);
        dtheta = atan2(ds*al,ds^2+t1*t2);
        %         if(d<0.)
        %             dtheta = -dtheta;
        %         end
        
        kkx = (dtheta*nx + log(r2/r1)*tx)/(pi2*kmat); % Equação
        % (A.7) com nx=1, ny=0, tx=0, ty=1.
        kky = (dtheta*ny + log(r2/r1)*ty)/(pi2*kmat); % Equação
        % (A.7) com nx=0, ny=1, tx=-1, ty=0.
        kk=kkx*mx+kky*my;
        
        hhx = (-(t2/r2^2-t1/r1^2)*nx ...
            + d*(1/r2^2-1/r1^2)*tx)/pi2; % Equação (A.8) com nx=1,
        % ny=0, tx=0, ty=1.
        hhy = (-(t2/r2^2-t1/r1^2)*ny ...
            + d*(1/r2^2-1/r1^2)*ty)/pi2; % Equação (A.8) com nx=0,
        % ny=1, tx=-1, ty=0.
        
        hh=hhx*mx+hhy*my;
        if(i==j)
            kk = -0.5; % Diagonal da matriz H
        end
        kk=-kk;
        if(bc(1,j)==0) % A temperatura é conhecida
            A(i,j) = -kk; % Os elementos de G vão para a matriz A
            b(i) = b(i) - hh*bc(2,j);% Os elementos de H vão para o vetor b
        else % O fluxo é conhecido
            A(i,j) = hh; % Os elementos de H vão para a matriz A
            b(i) = b(i) + kk*bc(2,j);% Os elementos de G vão para o vetor b
        end
        
    end
end
return