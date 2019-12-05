function [NOS,ELEM,CDC]=format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)

% Programa para formatação dos dados de entrada
%
%   Autor: Éder Lima de Albuquerque


cont_nos = 0;  % Counter to the physical nodes
cont_el = 0;	% Counter to the elements (the number of physical and geometric elements is the same).
num_lin = length(SEGMENTOS(:,1));	% Número de linhas no contorno
p_ini = SEGMENTOS(1,2);

%______________________________________________________________________
% Definition of the biggest dimension of the problem
max_dl = 0;
for lin = 1 : num_lin
    p1 = SEGMENTOS(lin,2);
    p2 = SEGMENTOS(lin,3);
    xp1 = PONTOS(p1,2);	yp1 = PONTOS(p1,3);
    xp2 = PONTOS(p2,2);	yp2 = PONTOS(p2,3);
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
    if dl > max_dl
        max_dl = dl;
    end;
end;
%_____________________________________________________________________



 no_ini=1;
 t=1;
p2=0;
no1_prox=0;
while(t<num_lin)  	% While over all lines
    while(p2~=p_ini)
        num_el_lin = MALHA(t,2);	% Number of the elements in the line t
        
        
        % Coordinates of the initial and final PONTOS of each line
        % (x1l,y1l,x2l,y2l)
        p1  = SEGMENTOS(t,2);
        p2  = SEGMENTOS(t,3);
        x1l = PONTOS(p1,2);
        y1l = PONTOS(p1,3);
        x2l = PONTOS(p2,2);
        y2l = PONTOS(p2,3);
        
        
        % 1. Generation of the matrices NOS, NOS_GEO, NOS_DRM, ELEM e ELEM_GEO
        if(SEGMENTOS(t,4)==0) % The segment is a straight line
            % Increment in x and y direction
            delta_x = x2l - x1l;
            delta_y = y2l - y1l;
        else %The segment is an arc
            % Compute the center of the arc and its coordinates
            r = SEGMENTOS(t,4);
            [xc,yc]=calcula_centro(x1l,y1l,x2l,y2l,r);
            % Distance between p1 and c (r1) and between p2 and c (r2)
            r1 = sqrt((x1l-xc)^2+(y1l-yc)^2);
            r2 = sqrt((x2l-xc)^2+(y2l-yc)^2);
            if abs(r1-r2)<.00001*max_dl
                % Compute the angle between the lines from point c to p1 (tet1) and c to p2 (tet2)
                [tet1,tet2] = calcula_arco(x1l,y1l,x2l,y2l,xc,yc);
                if tet2 < tet1
                    tet2 = tet2 + 2*pi;
                end;
                
                % Angle of the sector defined by the arc
                if SEGMENTOS(t,4) > 0
                    tet = abs(tet2-tet1);
                    sig = 1;
                else
                    tet = 2*pi-abs(tet2-tet1);
                    sig = -1;
                end;
                
                % Angle between two nodes of the line
                divtet = tet/(2*num_el_lin);
            else
                error('Error in the data input file: Wrong central point');
            end;
        end
        
        
        % Generation of elements and nodes
        for i = 1 : num_el_lin
            
            if(SEGMENTOS(t,4)==0) % The segment is a straight line
                x_i = x1l + delta_x/num_el_lin*(i-1);			% initial x coordinate of the element
                y_i = y1l + delta_y/num_el_lin*(i-1);			% initial y coordinate of the element
                x_m = x1l + delta_x/num_el_lin*(i-.5);	% midpoint x coordinate of the element
                y_m = y1l + delta_y/num_el_lin*(i-.5);	% midpoint y coordinate of the element
                lx=x_m-x_i;                              % distance in x direction between geometric nodes 1 and 2
                ly=y_m-y_i;                              % distance in y direction between geometirc nodes 1 and 2
                
            else  % The segment is an arc
                % Compute the node coordinates
                x_i = xc+r1*cos(tet1+2*(i-1)*sig*divtet);
                y_i = yc+r1*sin(tet1+2*(i-1)*sig*divtet);
                x_m = xc+r1*cos(tet1+(2*i-1)*sig*divtet);
                y_m = yc+r1*sin(tet1+(2*i-1)*sig*divtet);
            end;
            cont_el = cont_el + 1;
            % Set the coordinate of the physical and DRM nodes
            if(no1_prox==0) % the first node needs to be created
                cont_nos = cont_nos + 1;
                NOS(cont_nos,:)=[cont_nos,x_i,y_i];
                no1=cont_nos;
            else
                no1=no1_prox;
            end;
            if(p2~=p_ini || i<num_el_lin)
                cont_nos = cont_nos + 1;
                if(SEGMENTOS(t,4)==0) % Straight line
                    x_f=x_m+lx;
                    y_f=y_m+ly;
                else                 % arc
                    x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
                    y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
                end;
                NOS(cont_nos,:)=[cont_nos,x_f,y_f];
                no3=cont_nos;
                no1_prox=no3;
            else
                if(no_ini==0)
                    cont_nos = cont_nos + 1;
                    if(SEGMENTOS(t,4)==0) % Straight line
                        x_f=x_m+lx;
                        y_f=y_m+ly;
                    else                 % Arc
                        x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
                        y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
                    end;
                    NOS(cont_nos,:)=[cont_nos,x_f,y_f];
                    no3=cont_nos;
                    no1_prox=0;
                else
                    no3=no_ini;
                    no1_prox=0;
                end;
            end;
            ELEM(cont_el,:)=[cont_el,no1,no3];
        end;% end of for i = 1 : num_el_lin
        if p2 == p_ini
            if t < num_lin
                p_ini = SEGMENTOS(t+1,2);
                if(SEGMENTOS(t+1,3)==2)
                    no_ini = 0;
                else
                    no_ini = cont_nos+1;
                end
            end;
        end;
        t=t+1;
    end;                                  %end of while p2
end;                                   % end of while(t<num_lin)
%

% Geração da matriz CDC (Condições de Contorno)
% CDC = [n. do elemento, tipo de cdc, valor da cdc]
% Tipos de cdc: 0 : temperatura conhecida
%               1 : fluxo conhecido
cont_el2 = 0;

for l = 1 : length(SEGMENTOS(:,1))
    n_el_lin = MALHA(l,2);
    el_ini = cont_el2 + 1;
    el_fin = cont_el2 + n_el_lin;
    valorCDC=CCSeg(l,3);
    tipoCDC=CCSeg(l,2);
    for el = el_ini : el_fin
        CDC(el,:) = [el,tipoCDC,valorCDC,tipoCDC,valorCDC];
    end;
    cont_el2 = el_fin;
end;
return
