function mostra_elementos(SEGMENTOS,PONTOS,ELEM,NOS)

% Programa de visualização da geometria e malha
nelem=length(ELEM(:,1));
figure
axes('DataAspectRatio',[1,1,1]);
% axis off;
hold on

% Definition of the biggest dimension of the problem
max_dl = 0;
for lin = 1 : length(SEGMENTOS(:,1))
    p1 = SEGMENTOS(lin,2);
    p2 = SEGMENTOS(lin,3);
    xp1 = PONTOS(p1,2);	yp1 = PONTOS(p1,3);
    xp2 = PONTOS(p2,2);	yp2 = PONTOS(p2,3);
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
    if dl > max_dl
        max_dl = dl;
    end;
end;


fat = .02*max_dl;		% Scale factor

% Plot two PONTOS in white color only to define the size of the dimension 
% of the window in which the geometry will be shown

emin = min(PONTOS(:,2))-.1*max_dl;
dmax = max(PONTOS(:,2))+.1*max_dl;
imin = min(PONTOS(:,3))-.1*max_dl;
smax = max(PONTOS(:,3))+.1*max_dl;
plot(emin,imin,'w');	% Botton left point
plot(dmax,smax,'w');	% Top right point

% Plot the PONTOS that defines the geometry as a black x.
for i = 1 : length(PONTOS(:,1))
    plot(PONTOS(i,2),PONTOS(i,3),'kx','markersize',8);
end;

% Compute the continuous shape functions
for el = 1 : nelem
    
    no1 = ELEM(el,2);
    no2 = ELEM(el,3);
    
    x1 = NOS(no1,2);	y1 = NOS(no1,3);
    x2 = NOS(no2,2);	y2 = NOS(no2,3);
    xno=(x1+x2)/2;
    yno=(y1+y2)/2;

    Lx=x2-x1;
    Ly=y2-y1;
    L=sqrt(Lx^2+Ly^2);
    sx=Lx/L; % Component x do vetor tangente
    sy=Ly/L; % Componente y do vetor tangente
    n1=sy;
    n2=-sx;
    
    
    
    % Compute the PONTOS to interpolate the geometry over the element
    x_el = [x1 x2];
    y_el = [y1 y2];
    
    % Plot the element
    plot(x_el,y_el,'k-','markersize',6);	% Plot the node of the elements
    plot(xno,yno,'ro-','markersize',2);	% Plot the node of the elements
    % Plot a line in the beginning and in the end of the element
    plot([x_el(1)+fat*n1 x_el(1)-fat*n1], ...
        [y_el(1)+fat*n2 y_el(1)-fat*n2],'LineWidth',1.2,'color', ...
        'black','LineStyle','-');
    plot([x_el(2)+fat*n1 x_el(2)-fat*n1], ...
        [y_el(2)+fat*n2 y_el(2)-fat*n2],'LineWidth',1.2,'color',...
        'black','LineStyle','-');
end;
return

