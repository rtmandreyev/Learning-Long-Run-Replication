function bn = factorLoading(rhoQ,n)

K = size(rhoQ,1);
bn = NaN(K,n);
bn(:,1) = zeros(K,1);

for ii = 2:n
    bn(:,ii) = rhoQ'*bn(:,ii-1) - ones(K,1);
end

bn = bn(:,n);