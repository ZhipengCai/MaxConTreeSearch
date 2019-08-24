function [A, b] = genMatrixLinearizeFundamental(x1, x2)

A = [];
b = [];
N = size(x1, 2);
for i=1:N
    x = x1(1,i);
    y = x1(2,i);
    
    xp = x2(1,i);
    yp = x2(2,i);
    
    a = [xp*x xp*y xp yp*x yp*y yp x 1];
    A = [A; a];
    b = [b; -y];    
end


end