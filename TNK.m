function [obj, con] = TNK(x)
% ----------------------------------------------------------------------------
%Deb, K.,Pratap, A.,Agarwal, S., et al. A fast and elitist multiobjective
%genetic algorithm: NSGA-II. IEEE Transactions on Evolutionary Computation
% 2002, 6 (2): 182-197.
% x1 = [0,5]; x2 = [0,3].
% ----------------------------------------------------------------------------
x1 = x(:,1); x2 = x(:, 2);
obj(:, 1) =x1;
obj(:, 2) =x2;
con(:, 1) = -x1.^2 - x2.^2 +1 + 0.1*cos(16*atan(x1./x2));
con(:, 2) = (x1 - 0.5).^2 + (x2 - 0.5).^2 -0.5;


end