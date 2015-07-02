function ld = logdet(A)

[L,U,P] = lu(A);

for i = 1:size(P,1)
    p(i) = find(P(i,:));
end
s = 1;
for i = 1:size(P,1)
    while p(i) ~= i
        old = p(p(i));
        p(p(i)) = p(i);
        p(i) = old;
        s = s*-1;
    end
end

u = diag(U);
sgn = prod(sign(u))*s;
ld = sum(log(abs(diag(U)))) + log(sgn);
