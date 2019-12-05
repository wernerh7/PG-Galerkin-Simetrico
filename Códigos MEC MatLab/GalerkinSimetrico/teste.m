TT=XY(:,2);
qq=[zeros(ne,1);-ones(ne,1);zeros(ne,1);ones(ne,1); ...
   zeros(ne,1);-ones(ne,1);zeros(ne,1);ones(ne,1) ];

R1=G*qq-H1*TT
% [G*qq H1*TT]
R2=H2*qq-h*TT
% [H2*qq h*TT]
R3=H2*q'-h*T'

