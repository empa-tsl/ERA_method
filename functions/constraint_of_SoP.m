function [c, ceq] = constraint_of_SoP(x, SoP_min, SoP_max)
% non linear ineuality constraint c(x) <= 0
n = size(SoP_min,1);
c(1:n,1) = 1000 * (SoP_min - x / sum(x));        % SoP_min <= ERA / sum(ERA)
c(n+1:n+n,1) = 1000 * (x / sum(x) - SoP_max);    % SoP_max >= ERA / sum(ERA)

% non linear equality constraint ceq(x) = 0
ceq = [];
end