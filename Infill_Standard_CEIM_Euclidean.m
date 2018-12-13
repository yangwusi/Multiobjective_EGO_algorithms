function obj = Infill_Standard_CEIM_Euclidean(x, kriging_obj, kriging_con, non_dominated_front)
%-----------------------------------------------------
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill 
% Criteria for Expensive Multiobjective Optimization, IEEE Transactions 
% on Evolutionary Computation, 2017, 21 (6): 956-975.
%-----------------------------------------------------
% the input parameters
f = non_dominated_front;
% number of non-dominated points
% number of objectives
[num_pareto,num_obj] = size(f);
% number of input designs
num_x = size(x,1);
y = zeros(num_x,1);
%-----------------------------------------------------
% the kriging prediction and varince
u = zeros(num_x,num_obj);
mse = zeros(num_x,num_obj);
for ii = 1 : num_obj
    [u(:, ii),mse(:, ii)] = predictor(x, kriging_obj{ii});
end
s = sqrt(max(0,mse));
%-----------------------------------------------------
for ii = 1 : num_x
    u_matrix = repmat(u(ii,:),num_pareto,1);
    s_matrix = repmat(s(ii,:),num_pareto,1);
    EIM = (f - u_matrix).*Gaussian_CDF((f - u_matrix)./s_matrix) + s_matrix.*Gaussian_PDF((f - u_matrix)./s_matrix);
   y(ii) = min(sqrt(sum(EIM.^2,2)));
end
%---------------------------------------------------
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u_g = zeros(size(x,1), num_con);
mse_g = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u_g(:, ii), mse_g(:, ii)] = predictor(x, kriging_con{ii});
end
s_g = sqrt(max(0,mse_g));
% the PoF value
PoF = prod(Gaussian_CDF((0-u_g)./s_g), 2);
%-----------------------------------------------------
% the objective is maximized
obj = -y.*PoF;
end
