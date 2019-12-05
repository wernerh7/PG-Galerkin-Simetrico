function [T,q]= monta_Teq(CDC,T_PR,x)
% Separa fluxo e temperatura

% ncdc = n�mero de linhas da matriz CDC
% T = vetor que cont�m as temperaturas nos n�s
% q = vetor que cont�m o fluxo nos n�s

ncdc = length(CDC(:,1));
ieq=0;
for i=1:ncdc % La�o sobre as condi��es de contorno
    tipoCDC=CDC(i,2); % Tipo da condi��o de contorno
    valorCDC=CDC(i,3); % Valor da condi��o de contorno
    if tipoCDC == 1 % Fluxo � conhecido
        q(i) = valorCDC; % O fluxo � a condi�ao de contorno
    else % A temperatura � conhecida
        ieq=ieq+1;
        valorcalculado=x(ieq); % Valor que antes era desconhecido
        q(i) = valorcalculado; % O fluxo � o valor calculado
    end
end
for i=1:ncdc % Laço sobre as condições de contorno
    if(T_PR(i,6)==1&&T_PR(i,8)==1)% O fluxo é conhecido (a temperatura foi calculada)
        ieq=ieq+1;
        valorcalculado=x(ieq); % Valor que antes era desconhecido
        T(i) = valorcalculado; % O fluxo � o valor calculado
    else
        if(T_PR(i,6)==0)
            valorCDC=T_PR(i,7);
        else
            valorCDC=T_PR(i,9);
        end
        T(i) = valorCDC; % O fluxo � o valor calculado        
    end
end
return


