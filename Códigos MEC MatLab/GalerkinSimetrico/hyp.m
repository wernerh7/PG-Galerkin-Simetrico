function [h,g] = hyp(xc,yc,elem,conect)
%	SYMMETRIC-GALERKIN 2D LAPLACE EQUATION
%
%	HYP equation
%
%	Linear element --- Nonsingular integrals
%
eln=elem;
fac = -1.0/(pi+pi);
ng=30;
nelem=size(elem,2);
neu=nelem;
[xi,wg]=Gauss_Legendre(0,1,ng);

h = zeros(neu,neu); % inicializa a matriz h
g = zeros(neu,neu); % inicializa a matriz g


%
%	Loop over outer P element
%
for elp = 1:neu
    %
    %	Node numbers for outer element
    %
    pp(1) = eln(1,elp);
    pp(2) = eln(2,elp);
    %
    %	Singular elements
    %
    e0 = conect(1,pp(1));
    e1 = conect(2,pp(1));
    e2 = conect(2,pp(2));
    %
    %	Node coordinates
    %
    xp(1) = xc(pp(1));
    xp(2) = xc(pp(2));
    yp(1) = yc(pp(1));
    yp(2) = yc(pp(2));
    %
    %	Normal times jacobian
    %
    jnp1 = yp(2) - yp(1);
    jnp2 = -(xp(2) - xp(1));
    %
    %	Jacobian outer integral
    %
    jp = sqrt(jnp1^2 + jnp2^2);
    %
    %	Integration arrays
    %
    ip = zeros(2,2);
    %
    %	Store gauss/shape functions P integration
    %
    for gp = 1:ng
        x = xi(gp);
        %
        %	Shape functions in elem(elp) at the Gauss points
        %
        mp(1,gp) = 1.0 - x;
        mp(2,gp) = x;
        %
        %	Coordinates of Gauss P
        %
        xpg(gp) = xp(1)*mp(1,gp) + xp(2)*mp(2,gp);
        ypg(gp) = yp(1)*mp(1,gp) + yp(2)*mp(2,gp);
        %
        %	SINGLE INTEGRAL TERM
        %
        factor = wg(gp)*jp;
        ip(1,1) = ip(1,1) + factor*mp(1,gp)*mp(1,gp);
        ip(1,2) = ip(1,2) + factor*mp(1,gp)*mp(2,gp);
        ip(2,2) = ip(2,2) + factor*mp(2,gp)*mp(2,gp);
    end
    %
    %	SINGLE INTEGRAL
    %
    g(pp(1),pp(1)) = g(pp(1),pp(1)) - ip(1,1);
    g(pp(1),pp(2)) = g(pp(1),pp(2)) - ip(1,2);
    g(pp(2),pp(1)) = g(pp(2),pp(1)) - ip(1,2);
    g(pp(2),pp(2)) = g(pp(2),pp(2)) - ip(2,2);
    %
    %	Inner Q integration loop
    %
    for elq = 1:nelem
        if (elq ~= e0) && (elq ~= e1) && (elq ~= e2)
            %
            %	Node numbers for inner element
            %
            pq(1) = elem(1,elq);
            pq(2) = elem(2,elq);
            upper = 1;
            uppr2 = 1;
            %
            %	Inner element nodal coordinates
            %
            xq(1) = xc(pq(1));
            xq(2) = xc(pq(2));
            yq(1) = yc(pq(1));
            yq(2) = yc(pq(2));
            %
            %	Normal times jacobian
            %
            jnq1 = yq(2) - yq(1);
            jnq2 = -(xq(2) - xq(1));
            
            %
            %	Jacobian
            %
            jq = sqrt(jnq1^2 + jnq2^2);
            %
            %	Integration arrays
            %
            ip = zeros(2,2);
            ipn = zeros(2,2);
            %
            %	Q integration
            %
            for gq = 1:ng
                x = xi(gq);
                %
                %	Inner element shape functions
                %
                mq(1) = 1.0 - x;
                mq(2) = x;
                %
                %	Coordinates of Q
                %
                xqg = xq(1)*mq(1) + xq(2)*mq(2);
                yqg = yq(1)*mq(1) + yq(2)*mq(2);
                %
                %	Outer P integration
                %
                for gp = 1:ng
                    r1 = xqg - xpg(gp);
                    r2 = yqg - ypg(gp);
                    rsq = r1^2 + r2^2;
                    jnr = jnq1*r1 + jnq2*r2;
                    jNNr = jnp1*r1 + jnp2*r2;
                    %
                    %	Kernels
                    %
                    G_N = fac*(-jNNr/rsq);
                    G_Nn = fac*( -(jnq1*jnp1+jnq2*jnp2)/rsq + ...
                        2.0*jnr*jNNr/(rsq^2) );
                    factor = wg(gp)*wg(gq);
                    if uppr2 == 1
                        z_N = factor*mp(1,gp)*jq*G_N;
                        ipn(1,1) = ipn(1,1) + z_N*mq(1);
                        ipn(1,2) = ipn(1,2) + z_N*mq(2);
                        z_N = factor*mp(2,gp)*jq*G_N;
                        ipn(2,1) = ipn(2,1) + z_N*mq(1);
                        ipn(2,2) = ipn(2,2) + z_N*mq(2);
                    end
                    if upper == 1
                        z_Nn = factor*mp(1,gp)*G_Nn;
                        ip(1,1) = ip(1,1) + z_Nn*mq(1);
                        ip(1,2) = ip(1,2) + z_Nn*mq(2);
                        z_Nn = factor*mp(2,gp)*G_Nn;
                        ip(2,1) = ip(2,1) + z_Nn*mq(1);
                        ip(2,2) = ip(2,2) + z_Nn*mq(2);
                    end
                    %
                    %	End of P integration
                    %
                end
                %
                %	End of Q integration
                %
            end
            if uppr2 == 1
                g(pp(1),pq(1)) = g(pp(1),pq(1)) - ipn(1,1);
                g(pp(1),pq(2)) = g(pp(1),pq(2)) - ipn(1,2);
                g(pp(2),pq(1)) = g(pp(2),pq(1)) - ipn(2,1);
                g(pp(2),pq(2)) = g(pp(2),pq(2)) - ipn(2,2);
            end
            if upper == 1
                h(pp(1),pq(1)) = h(pp(1),pq(1)) + ip(1,1);
                h(pp(1),pq(2)) = h(pp(1),pq(2)) + ip(1,2);
                h(pp(2),pq(1)) = h(pp(2),pq(1)) + ip(2,1);
                h(pp(2),pq(2)) = h(pp(2),pq(2)) + ip(2,2);
            end
        end
    end
end
return
