function obj=Infill_Standard_Hypervolume_EIM(x, kriging_obj, non_dominated_front)
%-----------------------------------------------------
% The criterion is different from the code posted in the paper [1] beacuse
% I modified it to allow multiple design points be evaluated at the same
% time by using matrix opreation. Since the criterion is going to be evaluated 
% a large amount of times,  efficiency is very important.
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill 
% Criteria for Expensive Multiobjective Optimization, IEEE Transactions 
% on Evolutionary Computation, DOI: 10.1109/TEVC.2017.2697503
%-----------------------------------------------------
% the input parameters
f=non_dominated_front;
% number of non-dominated points
% number of objectives
[num_pareto,num_obj] = size(f);
% number of input designs
num_x= size(x,1);
r = 1.1*ones(1, num_obj);
%-----------------------------------------------------
% the kriging prediction and varince
u=zeros(num_x,num_obj);
mse=zeros(num_x,num_obj);
for ii=1:num_obj
    [u(:, ii),mse(:, ii)] = predictor(x, kriging_obj{ii});
end
s=sqrt(max(0,mse));
%-----------------------------------------------------
% the EI matrix (three dimensional matrix)
f_matrix  =  f .* ones(1,1,num_x);
u_matrix = reshape(u', 1, num_obj, num_x).* ones(num_pareto,1);
s_matrix =  reshape(s', 1, num_obj, num_x).* ones(num_pareto,1);
eim_matrix=(f_matrix-u_matrix).*GaussCDF((f_matrix-u_matrix)./s_matrix)+s_matrix.*GaussPDF((f_matrix-u_matrix)./s_matrix);
%-----------------------------------------------------
%  the Euclidean distance-based EI matrix criterion
y = reshape(min(prod(r-f+eim_matrix,2)-prod(r-f,2), [], 1), num_x, 1, 1);
%-----------------------------------------------------
% the objective is maximized
obj=-y;
end
