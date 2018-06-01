function [A] = factorVec(X)
% :(

A = [];
i = 1;
while ~isempty(X)
    a = X(1);
    A(i) = a;
    X = X(X~=a);
    i = i + 1;
end