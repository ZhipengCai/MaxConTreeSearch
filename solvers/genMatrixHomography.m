
function [A, b, c, d] = genMatrixHomography(x1, x2)

    nbpoints = size(x1,2);
    AA = [];
    bb = [];
    cc = [];
    dd = [];    
    for i=1:nbpoints
        x11 = x1(1,i); x12 = x1(2,i);
        xp1 = x2(1,i); xp2 = x2(2,i);
        ai1 = [x11 x12 1 0 0 0 -xp1*x11 -xp1*x12];
        ai2 = [0 0 0 x11 x12 1 -xp2*x11 -xp2*x12];
        AA = [AA; ai1; ai2];
        bb = [bb [-xp1; -xp2]];
        ci = [0 0 0 0 0 0 x11 x12];
        cc = [cc ci'];
        dd = [dd 1];               
    end    
    A = AA;
    b = bb;    
    c = cc;
    d = dd;
end