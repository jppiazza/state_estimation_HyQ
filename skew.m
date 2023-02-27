function [y] = skew(x)

% Convert a vector to its correponding skew symmetric operator.

y = [0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];

end

