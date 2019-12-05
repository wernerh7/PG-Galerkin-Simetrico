function [f,fx,fy]=domain_field(xfield,y,bc,node,dnorm,u,kmat)

% -------------------------------------------------------------------------------
%  Fucntion used to evaluate the temperature/heat flux using the direct
%     evaluation of the integral representation (dfdx and dfdy not implemented)
%
% -------------------------------------------------------------------------------

pi2 = pi*2;
nfield=size(xfield,2);
f=zeros(1,nfield);
fx=zeros(1,nfield);
fy=zeros(1,nfield);
n=size(y,2);

for j=1:n
    if(bc(1,j)==0) % A temperatura � conhecida
        f0  = bc(2,j);
        df0 = u(j);
    else % O fluxo � conhecido
        f0  = u(j);
        df0 = bc(2,j);
    end

    % Comprimento do elemento
    al= sqrt((y(1,node(2,j))-y(1,node(1,j)))^2+ ...
        (y(2,node(2,j))-y(2,node(1,j)))^2); 

    for i=1:nfield
        % Dist�ncias entre o ponto interno e os extremos dos elementos
        x11 = y(1,node(1,j)) - xfield(1,i);
        x21 = y(2,node(1,j)) - xfield(2,i);
        x12 = y(1,node(2,j)) - xfield(1,i);
        x22 = y(2,node(2,j)) - xfield(2,i);
        r1 =  sqrt(x11^2 + x21^2); % Dist�ncia para o in�cio do elemento
        r2 =  sqrt(x12^2 + x22^2); % Dist�ncia para o final do elemento
        
        % Proje��o do vetor dist�ncia no vetor normal ao elemento
        d  =  x11*dnorm(1,j) + x21*dnorm(2,j); % Figura A.1, p�gina 178
        t1 = -x11*dnorm(2,j) + x21*dnorm(1,j); % Dist�ncia T1 da figura A.1
        t2 = -x12*dnorm(2,j) + x22*dnorm(1,j); % Dist�ncia T2 da figura A.1

        ds = abs(d);
        dtheta = atan2(ds*al,ds^2+t1*t2);
 
        g = -(-dtheta*ds + al + t1*log(r1)-t2*log(r2))/(pi2*kmat); % Equa��o (A.5)

        kkx = (dtheta*dnorm(1,j) - log(r2/r1)*dnorm(2,j))/(pi2*kmat); % Equa��o
          % (A.7) com nx=1, ny=0, tx=0, ty=1.
        kky = (dtheta*dnorm(2,j) + log(r2/r1)*dnorm(1,j))/(pi2*kmat); % Equa��o
          % (A.7) com nx=0, ny=1, tx=-1, ty=0.
          
        hhx = -(-(t2/r2^2-t1/r1^2)*dnorm(1,j) ...
            - d*(1/r2^2-1/r1^2)*dnorm(2,j))/pi2; % Equa��o (A.8) com nx=1,
                 % ny=0, tx=0, ty=1. 
        hhy = -(-(t2/r2^2-t1/r1^2)*dnorm(2,j) ...
            + d*(1/r2^2-1/r1^2)*dnorm(1,j))/pi2; % Equa��o (A.8) com nx=0, 
                 % ny=1, tx=-1, ty=0.
             
             
        if(d<0.)
            dtheta = -dtheta;
        end
             
        h = -dtheta/pi2; % Equa��o (A.6)

        f(i) = f(i) + g*df0 - h*f0; % Integral (2.12) com o termo de 
           % dom�nio igual a zero.
        fx(i) = fx(i) + kkx*df0 - hhx*f0; % Integral (2.12) com o termo de 
           % dom�nio igual a zero.
        fy(i) = fy(i) + kky*df0 - hhy*f0; % Integral (2.12) com o termo de 
           % dom�nio igual a zero.

    end
end
return