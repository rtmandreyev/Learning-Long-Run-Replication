function [pr] = ECDF(X,x)

[B,m,n] = size(X);

if m>1 && n>1
    pr = NaN(m,n);
    for i = 1:m
        for j = 1:n
            pr(i,j) = sum(X(:,i,j) <= x(i,j))./B;
        end
    end
else
    pr = sum(X <= x)./B;
end

end