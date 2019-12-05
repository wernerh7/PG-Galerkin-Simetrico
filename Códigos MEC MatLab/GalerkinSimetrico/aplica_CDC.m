function [A11,A12,A21,A22,b1,b2]= aplica_CDC(G1,H1,G2,H2,CDC,T_PR)
%Aplica as condições de contorno trocando as colunas das matrizes H e G
ncdc = length(CDC(:,1)); % número de linhas da matriz CDC
ieq=0; % número da linha da matriz A
A11=[];
A12=[];
A21=[];
A22=[];
B11=[];
B12=[];
B21=[];
B22=[];
b1=[];
b2=[];
% gdl onde as temperaturas são desconhecidas (equação hipersingular)
for i=1:ncdc % Laço sobre as condições de contorno
    if(T_PR(i,6)==1&&T_PR(i,8)==1)% O fluxo é conhecido (calcular temperatura)
        ieq=ieq+1;
        b2(ieq)=0;
        % verifica temperatura
        icol1=0;
        icol2=0;
        for j=1:ncdc
            if(T_PR(j,6)==1&&T_PR(j,8)==1)% O fluxo é conhecido
                icol1=icol1+1;
                A22(ieq,icol1)=H2(i,j); % A matriz A recebe a coluna da matriz H
            else
                if(T_PR(j,6)==0)
                    valorCDC=T_PR(j,7);
                else
                    valorCDC=T_PR(j,9);
                end
                icol2=icol2+1;
                B22(ieq,icol2)=-H2(i,j); % A matriz B recebe a coluna da matriz H
                b2(ieq)=b2(ieq)-H2(i,j)*valorCDC;
            end
        end
        % verifica fluxo
        icol1=0;
        icol2=0;
        for j=1:ncdc
            if(CDC(j,1)==1)% O fluxo é conhecido
                valorCDC=CDC(j,2);
                icol1=icol1+1;
                B21(ieq,icol1)=G2(i,j); % A matriz B recebe a coluna da matriz H
                b2(ieq)=b2(ieq)+G2(i,j)*valorCDC;
            else
                icol2=icol2+1;
                A21(ieq,icol2)=-G2(i,j); % A matriz B recebe a coluna da matriz H
            end
        end
        
    end
end

% gdl onde os fluxos são desconhecidos (equação singular)
ieq=0;
for i=1:ncdc % Laço sobre as condições de contorno
    if(CDC(i,1)==0)% O fluxo é desconhecido (calcular temperatura)
        ieq=ieq+1;
        b1(ieq)=0;
        % verifica temperatura
        icol1=0;
        icol2=0;
        for j=1:ncdc
            if(T_PR(j,6)==1&&T_PR(j,8)==1)% O fluxo é conhecido
                icol1=icol1+1;
                A12(ieq,icol1)=H1(i,j); % A matriz B recebe a coluna da matriz H
            else
                if(T_PR(j,6)==0)
                    valorCDC=T_PR(j,7);
                else
                    valorCDC=T_PR(j,9);
                end
                icol2=icol2+1;
                B12(ieq,icol2)=-H1(i,j); % A matriz B recebe a coluna da matriz H
                b1(ieq)=b1(ieq)-H1(i,j)*valorCDC;
            end
        end
        % verifica fluxo
        icol1=0;
        icol2=0;
        for j=1:ncdc
            if(CDC(j,1)==1)% O fluxo é conhecido
                valorCDC=CDC(j,2);
                icol1=icol1+1;
                B11(ieq,icol1)=G1(i,j); % A matriz B recebe a coluna da matriz H
                b1(ieq)=b1(ieq)+G1(i,j)*valorCDC;
            else
                icol2=icol2+1;
                A11(ieq,icol2)=-G1(i,j); % A matriz B recebe a coluna da matriz H
            end
        end
        
    end
end
return