function obj = Infill_Pseudo_CEIM_Maximin(x, kriging_obj, kriging_con, f, point_added)
% the criterion is going to be evaluated a large amount of times, so efficiency is very important
%-----------------------------------------------------
% number of objectives
[num_pareto,num_obj] = size(f);
% number of input designs
num_x = size(x,1);
y = zeros(num_x,1);
%-----------------------------------------------------
% the kriging prediction and varince
u = zeros(num_x,num_obj);
mse = zeros(num_x,num_obj);
for ii=1:num_obj
    [u(:, ii),mse(:, ii)] = predictor(x, kriging_obj{ii});
end
s = sqrt(max(0,mse));
%-----------------------------------------------------
for ii = 1 : num_x
    u_matrix = repmat(u(ii,:),num_pareto,1);
    s_matrix = repmat(s(ii,:),num_pareto,1);
    EIM = (f - u_matrix).*Gaussian_CDF((f - u_matrix)./s_matrix) + s_matrix.*Gaussian_PDF((f - u_matrix)./s_matrix);
    % if this is the first infill point
    if ~isempty(point_added)
        num_added = size(point_added,1);
        correlation = zeros(num_added, num_obj);
        % the scaling of x is the same for different objectives
        scaled_x = (x(ii,:) - kriging_obj{1}.Ssc(1,:)) ./ kriging_obj{1}.Ssc(2,:);
        scaled_point_added = (point_added - repmat(kriging_obj{1}.Ssc(1,:),size(point_added,1),1)) ./ repmat(kriging_obj{1}.Ssc(2,:),size(point_added,1),1);
        dx = repmat(scaled_x,num_added,1) - scaled_point_added;
        % calculate the correlation
        for jj = 1:num_obj
            correlation(:,jj)  = feval(kriging_obj{jj}.corr, kriging_obj{jj}.theta, dx);
        end
        % the Pseudo EI matrix
        EIM = EIM.*repmat(prod(1-correlation,1),num_pareto,1);
    end
   y(ii) = min(max(EIM,[],2));
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
