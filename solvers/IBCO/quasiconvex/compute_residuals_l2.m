function [resn, s, inliers] = compute_residuals_l2(A, b, c, d, x, th, normalizeH)

    if (nargin < 7); normalizeH = 0; end;        
    if (nargin<6);   th = 1;   end;    
    
    if (normalizeH == 1) && (length(x)==9); x = x./x(end); end; 
    if (size(A,2)<length(x))
        x=x(1:size(A,2));
    end
    
    nbrimages = size(b,2);
    d1 = size(A,1)/nbrimages;
    AxPb = reshape(A*x,d1,nbrimages) + b; 
    Sqrt_AxPb = sqrt(sum(AxPb.^2));       
    if(numel(d) > 0 && numel(c) > 0)
        CxPd = (x'*c + d);
    elseif(numel(c)>0)
        CxPd = x'*c;
    else
        CxPd = ones(1,nbrimages);
    end
    resn = zeros(nbrimages, 1);
    id = abs(CxPd)>1e-6;    
    resn(id) = Sqrt_AxPb(id)./CxPd(id);
    
    resn(~id) = 100000000*max(resn);
    inliers = find(resn <= th & resn >= 0); 
    s = max((Sqrt_AxPb./th - CxPd)',0); 
end
