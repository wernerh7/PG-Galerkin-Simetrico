function [A,b,T_PR]=aplica_CDC(G,H,CDC,ELEM)
% Aplica as condi��es de contorno (troca as colunas das matrizes G e H)

n_el=size(CDC,1);
todos_valores=zeros(1,2*n_el);
cont=0;

% Cria a vari�vel T_PR que cont�m os n�s onde a temperatura � prescrita
% (conhecida).
% T_PR tem 5 colunas e o n�mero de linhas � igual ao n�mero de n�s onde a
% temperatura � conhecida.
% T_PR=[a1,a2,a3,a4,a5]
% a1 = n�mero do n� com temperatura prescrita.
% a2 = n�mero do elemento com temperatura prescrita ao qual este n�
% pertence.
% a3 = n�mero local do n� neste elemento.
% a4: caso a temperatura tamb�m seja prescrita no segundo elemento a que
% este n� pertence, ent�o, a4 conter� o n�mero deste elemento, caso
% contr�rio, conter� zero.
% a5: caso a4 seja diferente de zero, a5 conter� o n�mero local do n� no
% segundo elemento, caso contr�rio, conter� zero.

for el=1:n_el % for sobre os elementos para criar T_PR
    for no=1:2 % for sobre os n�s do elemento el
        no_global=ELEM(el,no+1); % n�mero global do n�
        tipoCDC=CDC(el,2*no);  % tipo da condi��o de contorno
        valorCDC=CDC(el,2*no+1); % valor da condi��o de contorno
        todos_valores(2*el-2+no)=valorCDC; % armazerna o valor da condi��o
                   % de contorno no vetor todos_valores
        if(tipoCDC==0) % se a temperatura � conhecida
            compartilha=0; % por enquanto n�o se sabe se a tempertura � 
            % tamb�m conhecida no segundo elemento a que este n� pertence
            i=1;
            while (~compartilha&&i<=cont&&cont>0) % verifica se o n� global
                % j� est� presente em T_PR. Quando ele encontra,
                % compartilha se torna igual a 1 e o while p�ra.
                if(no_global==T_PR(i,1)) % Se sim, o n� global j� est� 
                    % presente em T_PR. As colunas 4 e 5 s�o preenchidas e
                    % compartilha se torna igual a 1.
                    compartilha=1;
                    T_PR(i,4)=el;
                    T_PR(i,5)=no;
                end
                i=i+1;
            end
            if(~compartilha) % Se compartilha continua zero, ent�o o n� 
                % ainda n�o foi inserido em T_PR. Neste caso, as tr�s
                % primeiras colunas de T_PR s�o preenchidas.
                cont=cont+1;
                T_PR(cont,1)=no_global;
                T_PR(cont,2)=el;
                T_PR(cont,3)=no;
            end
        end
    end
end

n_temp_pr = length(T_PR(:,1)); % N�mero de n�s com temperatura conhecidas
for i=1 : n_temp_pr % for sobre os n�s com temperatura conhecida
    i_no=T_PR(i,1); % n�mero do n� com temperatura conhecida
    i_el=T_PR(i,2); % primeiro elemento com temperatura prescrita 
                        % que cont�m o n�
    i_no_loc=T_PR(i,3); % n�mero local do n� neste elemento
    ind_H=i_no; % �ndice da coluna da matriz H que ser� trocada
    ind_G=2*i_el+i_no_loc-2; % �ndice da coluna da matriz G que ser� 
                 %  trocada
    troca = G(:,ind_G); % armazena a coluna de G que ser� trocada
    G(:,ind_G) = -H(:,ind_H); % substitui a coluna de H na de G
    H(:,ind_H) = -troca; % Substitui a de G na de H
    if(T_PR(i,4)~=0) % Se a temperatura tamb�m � conhecida no segundo 
         %   elemento
        i_el=T_PR(i,4); % N�mero do segundo elemento 
        i_no_loc=T_PR(i,5); % n�mero local deste n� no segundo elemento
        ind_G=2*i_el+i_no_loc-2; % �ndice da coluna G que ser� atribu�do 
                  % zero
        H(:,ind_H)=H(:,ind_H)-G(:,ind_G); % soma o valor da coluna G na 
                  % coluna H
        G(:,ind_G)=0; % atribui zero na coluna G
    end;
end;
b=G*todos_valores'; % C�lculo do vetor b
A=H;