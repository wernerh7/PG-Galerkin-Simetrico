function T_PR=calc_T_PR(CDC,ELEM)
% Aplica as condi��es de contorno (troca as colunas das matrizes G e H)

n_el=size(CDC,1);
cont=0;

% Cria a vari�vel T_PR que cont�m os n�s onde a temperatura � prescrita
% (conhecida).
% T_PR tem 5 colunas e o n�mero de linhas � igual ao n�mero de n�s onde a
% temperatura � conhecida.
% T_PR=[a1,a2,a3,a4,a5,a6]
% a1 = n�mero do n� com temperatura prescrita.
% a2 = n�mero do elemento com temperatura prescrita ao qual este n�
% pertence.
% a3 = n�mero local do n� neste elemento.
% a4: caso a temperatura tamb�m seja prescrita no segundo elemento a que
% este n� pertence, ent�o, a4 conter� o n�mero deste elemento, caso
% contr�rio, conter� zero.
% a5: caso a4 seja diferente de zero, a5 conter� o n�mero local do n� no
% segundo elemento, caso contr�rio, conter� zero.
% a6: valor da temperatura prescrita
for el=1:n_el % for sobre os elementos para criar T_PR
    tipoCDC=CDC(el,1);  % tipo da condi��o de contorno
    valorCDC=CDC(el,2);  % valor da condi��o de contorno
    for no=1:2 % for sobre os n�s do elemento el
        no_global=ELEM(no,el); % n�mero global do n�
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
                T_PR(i,8)=tipoCDC;
                T_PR(i,9)=valorCDC;
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
            T_PR(i,6)=tipoCDC;
            T_PR(i,7)=valorCDC;
        end
    end
end
return