function obj=Infill_Pseudo_EIM_Maximin(x, kriging_obj, f, point_added)
% the criterion is going to be evaluated a large amount of times, so efficiency is very important
%-----------------------------------------------------
% number of objectives
[num_pareto,num_obj] = size(f);
% number of input designs
[num_x, num_vari] = size(x);
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
eim_matrix=(f_matrix-u_matrix).*Gaussian_CDF((f_matrix-u_matrix)./s_matrix)+s_matrix.*Gaussian_PDF((f_matrix-u_matrix)./s_matrix);
%-----------------------------------------------------
% if this is the first infill point
if ~isempty(point_added)
    num_point_added = size(point_added,1);
    correlation=zeros(num_point_added, num_obj, num_x);
    % the scaling of x is the same for different objectives
    scaled_x = (x - kriging_obj{1}.Ssc(1,:)) ./ kriging_obj{1}.Ssc(2,:);
    scaled_point_added = (point_added - kriging_obj{1}.Ssc(1,:)) ./ kriging_obj{1}.Ssc(2,:);
    % form the three dimensional matrix
    scaled_x_matrix  = reshape(scaled_x', 1, num_vari, num_x).* ones(num_point_added,1);
    scaled_point_added_matrix = repmat(scaled_point_added,[1, 1, num_x]);
    dx = scaled_x_matrix - scaled_point_added_matrix;
    % calculate the correlation
    for ii=1:num_obj
        correlation(:,ii, :)  = feval(kriging_obj{ii}.corr, kriging_obj{ii}.theta, dx);
    end
    % the Pseudo EI matrix
    eim_matrix=eim_matrix.*prod(1-correlation,1);
end
%-----------------------------------------------------
%  the Maximin distance-based EI matrix criterion
y = reshape(min(max(eim_matrix, [], 2), [], 1), num_x, 1, 1);
%-----------------------------------------------------
% the objective is maximized
obj=-y;
end
