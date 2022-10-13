function Z = Trm_prod_index(C,V)

Z = {};
M = vertcat(C{:});
X = V(1)==M(:,1);
for ii = reshape(find(X),1,[])
    nestfun(M(~X,:),C(ii))
end
    function nestfun(W,A)
    if V(2)==A{end}(2)
        if isempty(Z) || numel(A)<size(Z,2)
            Z = A;
        else
            Z(end+1,:) = A;
        end
    else
        Y = W(:,1)==A{end}(2);
        for jj = reshape(find(Y),1,[])
            nestfun(W(~Y,:),[A,{W(jj,:)}])
        end
    end
    end
end