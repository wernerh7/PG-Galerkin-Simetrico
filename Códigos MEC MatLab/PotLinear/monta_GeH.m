function [G,H,q]=monta_GeH(ELEM,NOS,k,fc)

n_el = length(ELEM(:,1));	% Número total de elementos
n_nos = length(NOS(:,1));	% Número total de nós

% Inicialização das matrizes H e G
H = zeros(n_nos,n_nos);
G = zeros(n_nos,2*n_el);
q=zeros(n_nos,1);

npg=16;
[qsi,w]=Gauss_Legendre(-1,1,npg);


for i = 1 : n_nos	% Percorre os pontos fontes

    % Coordenadas (xf,yf) dos pontos fontes
    xd = NOS(i,2);
    yd = NOS(i,3);

    for j = 1 : n_el	% Percorre os elementos

        % Numeração dos dois nós do elemento i
        no1 = ELEM(j,2);
        no2 = ELEM(j,3);

        % Coordenadas dos três nós (x1,y1,x2,y2)
        x1 = NOS(no1,2);	y1 = NOS(no1,3);
        x2 = NOS(no2,2);	y2 = NOS(no2,3);

        % Cálculo das submatrizes h_el e g_el

        % Integração não singular
        [g,h] = calc_gh_nsing(x1,y1,x2,y2,xd,yd,k,qsi,w);

        if ((i == no1) || (i == no2))  % O nó j pertence ao elemento i
            gii=calc_g_sing(x1,y1,x2,y2,k); % Integração singular com ponto fonte na extremidade do elemento qsi0=-1 ou qsi0=1)
            if(i==no1)
                g(1)=gii;
            else
                g(2)=gii;
            end
        end
        for nolocal = 1 : 2
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
            H(i,noglobal) = H(i,noglobal) + h(nolocal);
        end;
        G(i,2*j-1:2*j) = g;
    end
    if(fc(1,1)~=0)
        q(i)=calc_q(xd,yd,fc,k);
    else
        q(i)=0;
    end
end


% Calculo dos termos da diagonal da matriz H (consideraçao de corpo a
% temperatura constante)
for m = 1 : n_nos
    H(m,m) = 0; %zera a diagonal principal
    for n = 1 : n_nos
        if n ~= m
            H(m,m) = H(m,m) - H(m,n);
        end;
    end;
end;
return
