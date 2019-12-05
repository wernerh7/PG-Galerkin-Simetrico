function PONTOS_INT=gera_p_in(NPX,NPY,PONTO,SEGMENTOS)
% Programa para cria��o de pontos internos a uma  geometria
% gen�rica formada por retas
%
%   Autor: Frederico Louren�o
%   Data de cria��o setembro de 1999
%   Revis�o 0.0

% Defini��o da �rea m�xima para cria��o de pontos internos
xmin = min(PONTO(:,2));
xmax = max(PONTO(:,2));
ymin = min(PONTO(:,3));
ymax = max(PONTO(:,3));
lx = xmax - xmin;		% Largura do ret�ngulo que cont�m a geometria
ly = ymax - ymin;		% Altura do ret�ngulo que cont�m a geometria
n_SEGMENTOSs = length(SEGMENTOS(:,1));

% Defini��o da maior SEGMENTOS do problema
max_dl = 0;
for lin = 1 : length(SEGMENTOS(:,1))
    p1 = SEGMENTOS(lin,2);
    p2 = SEGMENTOS(lin,3);
    xp1 = PONTO(p1,2);	yp1 = PONTO(p1,3);
    xp2 = PONTO(p2,2);	yp2 = PONTO(p2,3);
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
    if dl > max_dl
        max_dl = dl;
    end;
end;


d_min = 0.03*max_dl;	% Dist�ncia m�nima dos pontos internos ao contorno
npx = NPX+1;				% N. de pontos na horizontal
npy = NPY+1;				% N. de pontos na vertical

PONTOS_INT =[];
% Atribui��o dos pontos finais e iniciais das SEGMENTOSs aos
% vetores xl1, xl2, yl1 e yl2
for t = 1 : n_SEGMENTOSs		% Percorre todas as SEGMENTOSs
    xl1(t) = PONTO(SEGMENTOS(t,2),2);
    xl2(t) = PONTO(SEGMENTOS(t,3),2);
    yl1(t) = PONTO(SEGMENTOS(t,2),3);
    yl2(t) = PONTO(SEGMENTOS(t,3),3);
    raio(t)= SEGMENTOS(t,4);    
end;

npi = 0;	% Inicializa��o no n. de pontos internos

for i = 1 : npy
    % Cria��o do candidato a ponto interno (xpi,ypi)
    ypi = ymin + (ly/npy)*i;	% y dentro do ret�ngulo
    for j = 1 : npx
        xpi = xmin + (lx/npx)*j;	% x dentro do ret�ngulo
        
        % In�cio dos testes para valida��o do ponto interno
        
        % 1. Verificando se o ponto est� dentro da geometria
        [xpi,ponto] = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio);
        % 2. Verificando se o ponto est� muito pr�ximo do contorno
        if strcmp(ponto,'interno')
            aceita = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio);
        else
            aceita = 'n�o';
        end;
        
        % Armazenando os dados do ponto interno
        if strcmp(aceita,'sim')	% O ponto est� dentro da geometria e
            npi = npi + 1;		% e respeita a dist�ncia ao contorno
            PONTOS_INT(npi,:) = [npi,xpi,ypi];
        end;
    end;
end;

% Nesse ponto est�o calculados os pontos internos




