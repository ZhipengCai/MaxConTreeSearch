function idxOut = reshapeIdx(idxIn, d)
    if size(idxIn,1)~=1
        idxIn = idxIn';
    end
    idxOut = reshape(idxIn*d-(d-1:-1:0)',1,d*numel(idxIn));
end