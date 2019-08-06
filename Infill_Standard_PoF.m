function obj = Infill_Standard_PoF(x, kriging_con)
%------------------------------------
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u = zeros(size(x,1), num_con);
mse = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u(:, ii), mse(:, ii)] = predictor(x, kriging_con{ii});
end
s=sqrt(max(0,mse));
%------------------------------------
% the PoF value
PoF = prod(Gaussian_CDF((0-u)./s), 2);
%-----------------------------------
% the objective is maximized
obj = -PoF;
end
