A = ones(5,6);
b = zeros(1,6);
b(1:3) = 1;
A-b

% C = bsxfun(@minus, A,b)