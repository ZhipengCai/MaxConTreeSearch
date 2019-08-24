%initialize sedumi model for constrains in socp problem for sx (only need to reset model.obj in BCO function)
%generate model parameters for sedumi (only for A_t,c_t and K, need to generate b_t in BCO)
function modelsx = genModelSocp_BCO_sedumi(A,b,c,d,epsilon,n,d1,d2)
i = 1:n;
bb = b/epsilon;
AA = [A'/epsilon; zeros(n,size(A,1))];
%At
D = [zeros(d2,n); eye(n)];
modelsx.At = zeros(d2+n,n*(2+d1));
modelsx.At(:,1:n) = -D;
bii = [c;eye(n)];
for j = 1:n
    Ati = -[bii(:,j), AA(:,(d1*(j-1)+1):(d1*j))];
    modelsx.At(:,(n+(1+d1)*(j-1)+1):(n+(1+d1)*j)) = Ati;
end
%ct
modelsx.ct = zeros((2+d1)*n,1);
modelsx.ct((i-1)*(d1+1)+n+1,1) = d(i)';
for j = 1:d1
    modelsx.ct((i-1)*(d1+1)+n+1+j,1) = bb(j,i)';
end
modelsx.K.l = n;
modelsx.K.q = (1+d1)*ones(1,n);
end