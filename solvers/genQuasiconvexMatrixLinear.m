function [A, b, c, d] = genQuasiconvexMatrixLinear(x, y)
    
    N = size(x,1); 
    d = size(x,2);
    
    xp = x'; yp = -y';
    A = [xp; zeros(size(xp))]; A = reshape(A, d, N*2)';  
    b = [yp; zeros(1,length(yp))];
       
    c = repmat(zeros(d-1,1),1, N);
    c = [c; zeros(1, N)];
    
    d = ones(1, N);
    
    
        


end