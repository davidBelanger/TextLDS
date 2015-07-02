function m = MIL2mat(s)
%s is a struct that represents the inverse of something  implicitly. this
%instantiates that inverse explicitly as a matrix
if(~isstruct(s))
    m = s;
else
    m = MIL2mat(s.firstTerm)  - s.firstFactor*s.secondFactor;
end

