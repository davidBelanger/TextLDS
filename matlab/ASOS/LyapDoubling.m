function Theta = LyapDoubling( E, M, T )

Phi = E;
Theta = M;

for t = 1:T
    incr = Phi*Theta*Phi';
    Theta = Theta + incr;
    if norm(vec(incr)) < 0.000001*norm(vec(Theta))
        break;
    end
    Phi = Phi^2;
end
%t