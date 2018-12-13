function [obj, con] = BNH(x)
% ----------------------------------------------------------------------------
% the constrains are inactive
% Deb, K. Multi-objective optimization using evolutionary algorithms. 
% John Wiley & Sons, 2001.
% x1 = [0,5]; x2 = [0,3].
% ----------------------------------------------------------------------------
x1 = x(:,1); x2 = x(:, 2);
obj(:, 1) = 4*x1.^2 + 4*x2.^2;
obj(:, 2) =(x1 - 5).^2 + (x2 - 5).^2;
con(:, 1) = (x1 - 5).^2 + x2.^2 -25;
con(:, 2) = 7.7 - (x1 - 8).^2 - (x2 + 3).^2;


end