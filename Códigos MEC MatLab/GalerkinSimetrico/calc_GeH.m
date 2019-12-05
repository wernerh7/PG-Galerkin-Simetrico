function [H1,G2,G] = calc_GeH(NOS,ELEM,k)

n=size(NOS,2); % N�mero de n�s = N�mero de elementos

% Inicializa��o das matrizes H e G
H = zeros(n,n);
G = zeros(n,n);
M1 = zeros(n,n);
M2 = zeros(n,n);

npg=10;

[xi,w]=Gauss_Legendre(0,1,npg);

for i=1:n % Loop on source points (Row)
    % Coordenadas dos dois n�s (x1,y1,x2,y2)
    % N�mero global dos dois n�s do elemento j
    no1d = ELEM(1,i);
    no2d = ELEM(2,i);
    
    % Coordenadas dos dois n�s (x1,y1,x2,y2)
    xx1 = NOS(1,no1d);	yy1 = NOS(2,no1d);
    xx2 = NOS(1,no2d);	yy2 = NOS(2,no2d);
    dgamadxi=sqrt((xx2-xx1)^2+(yy2-yy1)^2);
    

    cij1 = -[1./2. 1./2.]*dgamadxi;
%     cij2=  -[1./2. 1./2.]*dgamadxi;
    cij2=  -0*[1. 0.]*dgamadxi;
    
%     H(i,no1d)=H(i,no1d)+cij(1);
%     H(i,no2d)=H(i,no2d)+cij(2);
    M1(i,no1d)=M1(i,no1d)+cij1(1);
    M1(i,no2d)=M1(i,no2d)+cij1(2);
    
    M2(i,no1d)=M2(i,no1d)+cij2(1);
    M2(i,no2d)=M2(i,no2d)+cij2(2);
    
    
    for j=1:n % Loop on elements (Column)
        
        % N�mero global dos dois n�s do elemento j
        no1 = ELEM(1,j);
        no2 = ELEM(2,j);
        
        % Coordenadas dos dois n�s (x1,y1,x2,y2)
        x1 = NOS(1,no1);	y1 = NOS(2,no1);
        x2 = NOS(1,no2);	y2 = NOS(2,no2);
        
        % C�lculo das submatrizes h_el e g_el
        
        al = sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
        sx = (x2-x1)/al;
        sy = (y2-y1)/al;
        nx = sy;  % Componente x do vetor normal ao elemento
        ny = -sx; % Componente y do vetor normal ao elemento
        noglobal1 = no1; %�ndice da matriz global H
        noglobal2 = no2; %�ndice da matriz global H
        for kk=1:npg % Laço sobre os pontos de integração
            N1=1-xi(kk);% Função de forma N1 => linear contínua
            N2=xi(kk);% Função de forma N2 => linear contínua
            xd=N1*xx1+N2*xx2; % Calcula a coordenada x do ponto de integração
            yd=N1*yy1+N2*yy2; % Calcula a coordenada y do ponto de integração
            % Compute parameters used in the formulas for the two intergals
            x11 = x1-xd; % diferen�a de coordenada x do
            % in�cio do elemento at� o ponto fonte
            x21 = y1-yd;% diferen�a de coordenada y do
            % in�cio do elemento at� o ponto fonte
            x12 = x2-xd;% diferen�a de coordenada x do
            % final do elemento at� o ponto fonte
            x22 = y2-yd;% diferen�a de coordenada y do
            % final do elemento at� o ponto fonte
            r1 = sqrt(x11^2 + x21^2); % Dist�ncia do ponto fonte ao in�cio
            % do elemento
            r2 = sqrt(x12^2 + x22^2); % Dist�ncia do ponto fonte ao final
            % do elemento
            % dnorm = vetor normal n1 e n2
            d = x11*nx + x21*ny; % = delta x1 n1+delta y1 n2
            % = dr1/dn = h
            %         t1 = real(-x11*ny + x21*nx); %=-delta x1 n2+delta y1 n1
            t1 = -x11*ny + x21*nx; %=-delta x1 n2+delta y1 n1
            % = dr1/dt
            %         t2 = real(-x12*ny + x22*nx);%=-delta x2 n2+delta y2 n1
            t2 = -x12*ny + x22*nx;%=-delta x2 n2+delta y2 n1
            ds = abs(d);
%             dtheta = atan2(ds*al,ds^2+t1*t2);
            a=[x11,x21,0];
            b=[x12,x22,0];
            dtheta = atan2(norm(cross(a,b)),dot(a,b));
            
            if(t1~=0)
                t1log=t1*log(r1);
            else
                t1log=0;
            end
            if(t2~=0)
                t2log=t2*log(r2);
            else
                t2log=0;
            end
            dslogr2r1=d*log(r2/r1);
            dslogr1r2=d*log(r1/r2);
            if(r2==0)
                r22log=0;
                dslogr2r1=0;
                dslogr1r2=0;
            else
                r22log=(r2)^2*log(r2);
            end
            if(r1==0)
                r12log=0;
                dslogr2r1=0;
                dslogr1r2=0;
            else
                r12log=(r1)^2*log(r1);
            end
            % Elementos da matriz G
            g1=(4*t2*(-ds*dtheta+t1log-t2log+al) ...
                -2*r12log+r1^2+2*r22log-r2^2)/(8*pi*k*al);
            
            g2=(4*t1*(ds*dtheta-t1log+t2log-al) ...
                +2*r12log-2*r22log+r2^2-r1^2)/(8*pi*k*al);
%             if(d<0)
            if(d<0&&abs(d)>1e-8)                
                dtheta = -dtheta;
            end
            % Elementos da matriz H
            h1=(t2*dtheta+dslogr1r2)/(2*pi*al);
            h2=(dslogr2r1-t1*dtheta)/(2*pi*al);
            H(i,noglobal1) = H(no1d,noglobal1) + dgamadxi*w(kk)*h1;
            H(i,noglobal2) = H(no1d,noglobal2) + dgamadxi*w(kk)*h2;
            G(i,noglobal1) = G(no1d,noglobal1) + dgamadxi*w(kk)*g1;
            G(i,noglobal2) = G(no1d,noglobal2) + dgamadxi*w(kk)*g2;
        end
    end
end
H1=H+M1; % H1 multiplica T (soma dos elementos de uma linha igual a zero)
G2=H'+M2; % H2 multiplica q (soma dos elementos de uma linha NÂO é zero)
return