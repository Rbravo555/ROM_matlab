function z = UpdateZMatrix(z,w, current_Dq, Current_basis_index)

number_of_bases = size(w,1);


for m=1:number_of_bases
    k=m+1;
    for p=k:number_of_bases
        %fprintf('%d,%d\n',i,j)
        z(m,p) = z(m,p) + w{Current_basis_index, m,p}'*current_Dq;
        z(p,m) = -z(m,p);
    end
end