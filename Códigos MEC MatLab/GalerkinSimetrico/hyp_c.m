function [h,g] = hyp_c(h,g,xc,yc,eln)
%	HYP equation
%
%	Linear element --- Coincident singular integrals
%
fac = -1.0/(pi+pi);
fac2 = fac/2.0;
fac6 = fac*(-pi/6.0);
neu=size(eln,2);

%
%	Loop over coincident element
%
for elp = 1:neu
    %
    %	Node numbers for outer element
    %
    pp(1) = eln(1,elp);
    pp(2) = eln(2,elp);
    %
    %	Node coordinates
    %
    xp(1) = xc(pp(1));
    xp(2) = xc(pp(2));
    yp(1) = yc(pp(1));
    yp(2) = yc(pp(2));
    jp = sqrt( (xp(2)-xp(1))^2 + (yp(2) - yp(1))^2 );
    %
    %	Formulas from mp_hypcf/mp_hypcp codes
    %
    factor = fac6*jp;
    %
    %	From maple code mp_hypcf and mp_hypcp
    %
    g(pp(1),pp(1)) = g(pp(1),pp(1)) + factor + factor;
    g(pp(1),pp(2)) = g(pp(1),pp(2)) + factor;
    g(pp(2),pp(1)) = g(pp(2),pp(1)) + factor;
    g(pp(2),pp(2)) = g(pp(2),pp(2)) + factor + factor;
    h(pp(1),pp(1)) = h(pp(1),pp(1)) - fac2;
    h(pp(1),pp(2)) = h(pp(1),pp(2)) + fac2;
    h(pp(2),pp(1)) = h(pp(2),pp(1)) + fac2;
    h(pp(2),pp(2)) = h(pp(2),pp(2)) - fac2;
end
return
