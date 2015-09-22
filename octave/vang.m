## this is an octave file
j = @(s, m) (s .^ (-1 / m) - 1) .^ (1 - m);
sminus = @(s,m,r) ( r ^ (1/(1-m)) * (s .^ (-1/m) - 1) + 1 ) .^ (-m);
jminus = @(y,m) (1+ (y .^ (1/(1-m)))) .^ (-m);
dsminus = @(s,m,r) \
(r ^ (1/(1-m)) *(s .^ (-1/m)-1)+1).^(-1-m) * r ^ (1/(1-m)) .* s .^ (-1/m - 1) 