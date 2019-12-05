function [T,q,q_vet] = Monta_Teq(x,CDC,T_PR)

% Function to reorder de vectors in accordance with boundary conditions


n_el_total = length(CDC(:,1));
n_temp_pr = length(T_PR(:,1));
cdc_val=zeros(1,2*n_el_total);
q=zeros(n_el_total,3);

% Generation of vector of boundady conditions
for el = 1 : n_el_total
    for no = 1 : 2
        cdc_val(2*(el-1)+no) = CDC(el,2*no+1);
    end;
end;
q_vet = cdc_val;
T=x;

for i=1 : n_temp_pr
    n_el=T_PR(i,2);    
    n_no=T_PR(i,1);
    n_no_loc=T_PR(i,3);
    ind_T=n_no;
    ind_q=2*n_el+n_no_loc-2;
    q_vet(ind_q) = x(ind_T);
    valor_x=x(ind_T);
    T(ind_T) = cdc_val(ind_q);
    if(T_PR(i,5)~=0)
        n_el=T_PR(i,4);
        n_no_loc=T_PR(i,5);
        ind_q=2*n_el+n_no_loc-2;
        q_vet(ind_q)=valor_x;
    end;
end;


for el = 1 : n_el_total
    q(el,1) = el;
    for no = 1 : 2
        q(el,no+1) = q_vet(2*(el-1)+no);
    end;
end;

return

