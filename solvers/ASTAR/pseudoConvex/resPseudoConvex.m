function r = resPseudoConvex(A,b,c,d,x,th,d1)

AxPb = reshape(A*x,d1,numel(d)) + b; % Ax + b
Sqrt_AxPb = sqrt(sum(AxPb.^2,1));                     % sqrt(Ax + b)
if(numel(d) > 0 && numel(c) > 0)
    CxPd = (x'*c + d);
elseif(numel(c)>0)
    CxPd = x'*c;
else
    CxPd = ones(1,numel(d));
end
r = Sqrt_AxPb/th-CxPd;
end