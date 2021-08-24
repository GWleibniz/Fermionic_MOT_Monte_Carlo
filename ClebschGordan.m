
% ClebschGordan.m by David Terr, Raytheon, 6-17-04
% Modified on 11-9-04

% ClebschGordan(j1,j2,j,m1,m2,m) returns the Clebsch-Gordan coefficient <j1,j2,m1,m2|j1,j2,j,m>. 
% This program requires first downloading Wigner3j.m.

function cg = ClebschGordan(j1,j2,j,m1,m2,m)

% error checking
if ( 2*j1 ~= floor(2*j1) || 2*j2 ~= floor(2*j2) || 2*j ~= floor(2*j) ...
        || 2*m1 ~= floor(2*m1) || 2*m2 ~= floor(2*m2) || 2*m ~= floor(2*m) )
    error('All arguments must be integers or half-integers.');
    return;
end

if m1 + m2 ~= m
    cg = 0;
    return;
end

if ( j1 - m1 ~= floor ( j1 - m1 ) )
    cg = 0;
    return;
end

if ( j2 - m2 ~= floor ( j2 - m2 ) )
    cg = 0;
    return;
end

if ( j - m ~= floor ( j - m ) )
    cg = 0;
    return;
end

if j > j1 + j2 || j < abs(j1 - j2)
    cg = 0;
    return;
end

if abs(m1) > j1
    cg = 0;
    return;
end

if abs(m2) > j2
    cg = 0;
    return;
end

if abs(m) > j
    cg = 0;
    return;
end

cg = (-1)^(j1-j2+m) * sqrt(2*j + 1) * Wigner3j(j1,j2,j,m1,m2,-m);


% Reference: Clebsch-Gordan Coefficient entry of Eric Weinstein's Mathworld: http://mathworld.wolfram.com/Clebsch-GordanCoefficient.html