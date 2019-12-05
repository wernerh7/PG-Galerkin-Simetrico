function [h,g] = hyp_ar(h,g,xc,yc,elem,conect)
%	HYP equation
%
%	Linear element --- ADJACENT singular integrals
%	REVERSED ORIENTATION
%
eln=elem;
fac = -1.0/(pi+pi);
pi4 = pi/4.0;
ng=30;
neu=size(elem,2);
[xi,wg]=Gauss_Legendre(0,1,ng);

% xi =[0.930568155797026   0.669990521792428   0.330009478207572   0.069431844202974];
% wg = [0.173927422568727   0.326072577431273   0.326072577431273   0.173927422568727];
% ng=4;

%
%	Loop over neumann surface
%
for elp =1:neu
    %
    %	Node numbers for P
    %	element [p1,p2] = P
    %
    pp(1) = eln(1,elp);
    pp(2) = eln(2,elp);
    elq = conect(2,pp(2));
    %
    %	[Q2,Q3] = Q;
    %
    pp(3) = elem(2,elq);
    Q2 = elem(1,elq);
    Q3 = elem(2,elq);
    %
    %	Node coordinates
    %
    xp(1) = xc(pp(1));
    xp(2) = xc(pp(2));
    xp(3) = xc(pp(3));
    yp(1) = yc(pp(1));
    yp(2) = yc(pp(2));
    yp(3) = yc(pp(3));
    %
    %	Jacobians
    %
    t1 = xp(2) - xp(1);
    t2 = yp(2) - yp(1);
    %
    %	P normal
    %
    jnp1 = t2;
    jnp2 = -t1;
    jp2 = t1^2 + t2^2;
    t1 = xp(3) - xp(2);
    t2 = yp(3) - yp(2);
    jnq1 = t2;
    jnq2 = -t1;
    jq2 = t1^2 + t2^2;
    pfac = fac*pi4;
    pnfac = pfac*sqrt(jq2);
    jnjN = jnq1*jnp1 + jnq2*jnp2;
    del1 = abs( yp(3)*(xp(2)-xp(1)) + (xp(1)-xp(3))*yp(2) + ...
        yp(1)*(xp(3)-xp(2)) );
    fint = log( ((xp(3)-xp(1))^2+(yp(3)-yp(1))^2)/jq2 )/2.0;
    %
    %	Integration arrays
    %
    ip = zeros(2,2);
    ipn = zeros(2,2);
    %
    %	Exactly integrated term
    %
    ip(2,1) = fac*fint;
    %
    %	Theta integration loop
    %
    for gp = 1:ng
        %
        %	[0,pi/4] and [pi/4,pi/2]
        %
        for jj = 1:2
            theta = xi(gp)*pi4 + (jj-1)*pi4;
            ct = cos(theta);
            st = sin(theta);
            L = 1.0/ct;
            if jj == 2
                L = 1.0/st;
            end
            L2 = L^2;
            %
            %	Flux integral
            %
            asq = jp2*ct*ct+jq2*st*st+ ...
                2.0*(-xp(2)^2+(xp(1)+xp(3))*xp(2)-yp(2)^2+ ...
                (yp(1)+yp(3))*yp(2)-xp(3)*xp(1)-yp(3)*yp(1))*st*ct;
            kerpn = ( ((xp(2)-xp(1))*ct+(xp(3)-xp(2))*st)*jnp1 + ...
                ((yp(2)-yp(1))*ct+(yp(3)-yp(2))*st)*jnp2 )/asq;
            c11 = -kerpn*( L2*ct*( -st*L/3.0 + 1.0/2.0 ) );
            c12 = -kerpn*( ct*st*L2*L/3.0 );
            c21 = -kerpn*( L*( ct*st*L2/3.0 - (st+ct)*L/2.0 + 1.0 ) );
            c22 = -kerpn*( L2*st*( -ct*L/3.0 + 1.0/2.0 ) );
            factor = wg(gp)*pnfac;
            ipn(1,1) = ipn(1,1) + factor*c11;
            ipn(1,2) = ipn(1,2) + factor*c12;
            ipn(2,1) = ipn(2,1) + factor*c21;
            ipn(2,2) = ipn(2,2) + factor*c22;
            %
            %	Potential integral
            %
            asq = (ct*(xp(1)-xp(2))+st*(xp(2)-xp(3)))^2 + ...
                (ct*(yp(1)-yp(2))+st*(yp(2)-yp(3)))^2;
            br1 = -L*(jnjN+2.0*del1*del1*st*ct/asq)/asq;
            br2 = L*br1/2.0;
            b11 = ct*(br1-st*br2);
            b12 = ct*st*br2;
            b21 = -(st+ct)*br1 + ct*st*br2;
            b22 = st*(br1-ct*br2);
            factor = wg(gp)*pfac;
            ip(1,1) = ip(1,1) + factor*b11;
            ip(1,2) = ip(1,2) + factor*b12;
            ip(2,1) = ip(2,1) + factor*b21;
            ip(2,2) = ip(2,2) + factor*b22;
            %
            %	End theta integration
            %
        end
    end
    g(pp(1),Q2) = g(pp(1),Q2) - ipn(1,1);
    g(pp(1),Q3) = g(pp(1),Q3) - ipn(1,2);
    h(pp(1),Q2) = h(pp(1),Q2) + ip(1,1);
    h(pp(1),Q3) = h(pp(1),Q3) + ip(1,2);
    g(pp(2),Q2) = g(pp(2),Q2) - ipn(2,1);
    g(pp(2),Q3) = g(pp(2),Q3) - ipn(2,2);
    h(pp(2),Q2) = h(pp(2),Q2) + ip(2,1);
    h(pp(2),Q3) = h(pp(2),Q3) + ip(2,2);
end
return
