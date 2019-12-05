function T_pint = calc_T_pint(PONTOS_int,NOS,ELEM,k,T,q)
n=size(NOS,1); % N�mero de n�s = N�mero de elementos
n_pint=size(PONTOS_int,1);
% Inicializa��o das matrizes H e G
T_pint=zeros(n_pint,1);
sum=0;
for i=1:n_pint % Loop on source points (Row)
    % Coordenadas dos dois n�s (x1,y1,x2,y2)
    % N�mero global dos dois n�s do elemento j
    xd = PONTOS_int(i,2);
    yd = PONTOS_int(i,3);
    for j=1:n % Loop on elements (Column)
        % N�mero global dos dois n�s do elemento j
        no1 = ELEM(j,2);
        no2 = ELEM(j,3);
        % Coordenadas dos dois n�s (x1,y1,x2,y2)
        x1 = NOS(no1,2);	y1 = NOS(no1,3);
        x2 = NOS(no2,2);	y2 = NOS(no2,3);
        % C�lculo das submatrizes h_el e g_el
        al = sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
        sx = (x2-x1)/al;
        sy = (y2-y1)/al;
        nx = sy;  % Componente x do vetor normal ao elemento
        ny = -sx; % Componente y do vetor normal ao elemento
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
        dtheta = atan2(ds*al,ds^2+t1*t2);
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
        if(d<0&&abs(d)>1e-8)
            dtheta = -dtheta;
        end
        % Elementos da matriz H
        h1=(t2*dtheta+dslogr1r2)/(2*pi*al);
        h2=(dslogr2r1-t1*dtheta)/(2*pi*al);
        T_pint(i) =T_pint(i) + h1*T(no1);
        T_pint(i) =T_pint(i) + h2*T(no2);
        sum=sum+h1+h2;
        T_pint(i) =T_pint(i) - g1*q(j);
        T_pint(i) =T_pint(i) - g2*q(j);
    end
end
return