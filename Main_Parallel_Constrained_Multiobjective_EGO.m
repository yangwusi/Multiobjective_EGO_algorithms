% 1. The parallel constrained multiobjective EGO algorithm using PCEIM 
% (pseudo constrained expected improvement matrix) criteria. For detailed 
% description about the EIM criteria, please refer to [1].
% 2. The dace toolbox [2] is used for building the Kriging models in the
%    implementations.
% 3. The non-dominated sorting method by Yi Cao [3] is used to identify the
%    non-dominated fronts from all the design points
% 4. The hypervolume indicators are calculated using the faster algorithm of
%    [4] Nicola Beume et al. (2009).
% 5. The EIM criteria are maximized by DE [5] algorithm.
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill
%      Criteria for Expensive Multiobjective Optimization. IEEE Transactions
%      on Evolutionary Computation, 2017, 21 (6): 956-975.
% [2] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging
%     Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%     Modelling, Technical University of Denmark, 2002. Available at:
%     http://www2.imm.dtu.dk/~hbn/dace/.
% [3] http://www.mathworks.com/matlabcentral/fileexchange/17251-
%      pareto-front.
% [4] N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold,
%     On the Complexity of Computing the Hypervolume Indicator, IEEE
%     Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
% [5] K. Price, R. M. Storn, and J. A. Lampinen, Differential evolution:
%     a practical approach to global optimization: Springer Science & Business Media, 2006.
%     http://www.icsi.berkeley.edu/~storn/code.html
% zhandawei@swjtu{dot}edu{dot}cn
% 2019.07.10 initial creation
% -----------------------------------------------------------------------------------------
clearvars;close all;
% settings of the problem
%  'SRN', 'TNK', 'BNH', 'Welded_Beam'
fun_name ='TNK';
% infill criterion: 'PCEIM_Euclidean','PCEIM_Maximin','PCEIM_Hypervolume'
infill_name= 'PCEIM_Hypervolume';
% number of initial design points
num_initial = 100;
% the maximum allowed evaluations
max_evaluation = 200;
% number of updating points selected in each cycle
num_q = 5;
% get the information about the problem
switch fun_name
    case 'SRN'
        num_obj = 2; num_con = 2; num_vari = 2; design_space=[-20,-20; 20, 20]; ref_point=[250, 50];
    case 'TNK'
        num_obj = 2; num_con = 2; num_vari = 2; design_space=[0,0; pi, pi]; ref_point=[1.2, 1.2];
    case 'BNH'
        num_obj = 2; num_con = 2; num_vari = 2; design_space=[0,0; 5,3]; ref_point=[140, 50];
    case 'Welded_Beam'
        num_obj = 2; num_con = 4; num_vari = 4; design_space=[0.125,0.1,0.1,0.125; 5,10,10,5]; ref_point=[100, 0.08];
end
% the intial design points
sample_x = repmat(design_space(1,:),num_initial,1) +  repmat(design_space(2,:)-design_space(1,:),num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
[sample_y, sample_g] = feval(fun_name, sample_x);
sample_y_scaled = (sample_y - repmat(min(sample_y),size(sample_y,1),1))./repmat(max(sample_y)-min(sample_y),size(sample_y,1),1);
% initialize some parameters
evaluation = size(sample_x,1);
kriging_obj = cell(1,num_obj);
kriging_con = cell(1,num_con);
hypervolume = zeros(max_evaluation - num_initial+1,1);
iteration = 0;
% check is there is at least one feasible solution
index = sum(sample_g <= 0, 2) == num_con;
if sum(index) == 0
    hypervolume(1) = 1E6;
    fprintf(' iteration: %d, evaluation: %d, hypervolume: no feasible solution found\n', 0, evaluation);
else
    feasible_y = sample_y(index, :);
    feasible_y_scaled = sample_y_scaled(index, :);
    index = Paretoset(feasible_y);
    non_dominated_front = feasible_y(index,:);
    non_dominated_front_scaled = feasible_y_scaled(index,:);
    hypervolume(1) = Hypervolume(non_dominated_front,ref_point);
    % plot current non-dominated front points
    if num_obj == 2
        scatter(non_dominated_front(:,1), non_dominated_front(:,2),'ro', 'filled');title(sprintf('%s problem, iteration: %d, evaluations: %d',fun_name, 0,evaluation));drawnow;
    elseif num_obj == 3
        scatter3(non_dominated_front(:,1), non_dominated_front(:,2),non_dominated_front(:,3),'ro', 'filled');title(sprintf('%s problem,iteration: %d, evaluations: %d',fun_name,0,evaluation));drawnow;
    end
    % print the hypervolume information
    fprintf(' iteration: %d, evaluation: %d, hypervolume: %f\n', 0, evaluation, hypervolume(1));
end
%-------------------------------------------------------------------------
% beginning of the iteration
while evaluation < max_evaluation
    % build the initial kriging model for each objective
    for ii = 1:num_obj
        kriging_obj{ii} = dacefit(sample_x,sample_y_scaled(:,ii),'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    % build the Kriging model for each constraint function
    for ii = 1: num_con
        kriging_con{ii} = dacefit(sample_x, sample_g(:, ii),'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    % find q points usign PCEIM criterion
    num_k = min(num_q,max_evaluation - evaluation);
    best_x = zeros(num_k, num_vari);
    point_added = [];
    for ii = 1 : num_k
        if sum(index) == 0
            infill_criterion = @(x)Infill_Standard_PoF(x, kriging_con);
        else
            switch infill_name
                case 'PCEIM_Euclidean'
                    infill_criterion = @(x)Infill_Pseudo_CEIM_Euclidean(x,kriging_obj,kriging_con,non_dominated_front_scaled,point_added);
                case 'PCEIM_Maximin'
                    infill_criterion = @(x)Infill_Pseudo_CEIM_Maximin(x,kriging_obj,kriging_con,non_dominated_front_scaled,point_added);
                case 'PCEIM_Hypervolume'
                    infill_criterion = @(x)Infill_Pseudo_CEIM_Hypervolume(x,kriging_obj,kriging_con,non_dominated_front_scaled,point_added);
                otherwise
                    error('you should select infill_name from PCEIM_Euclidean, PCEIM_Maximin, and PCEIM_Hypervolume');
            end
        end
        best_x(ii,:) = DE(infill_criterion, num_vari, design_space(1,:), design_space(2,:), 50, 200);
        % check if the candidate point is too close to sampled points
        if min(sqrt(sum((repmat(best_x(ii, :),size([sample_x;point_added],1),1)-[sample_x;point_added]).^2,2)))<1E-8
            best_x(ii,:) = DE(@(x)Infill_Maximal_Distance(x, [sample_x;point_added]), num_vari, design_space(1,:), design_space(2,:), 50, 200);
        end
        % update point_added
        point_added =  best_x(1:ii, :);
    end
    % evaluating the candidate with the real function
    [best_y, best_g] = feval(fun_name, best_x);
    % add the new points to the design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y; best_y];
    sample_g = [sample_g; best_g];
    sample_y_scaled =(sample_y - repmat(min(sample_y),size(sample_y,1),1))./repmat(max(sample_y)-min(sample_y),size(sample_y,1),1);
    iteration = iteration + 1;
    evaluation = evaluation + size(best_x,1);
    % calculate the hypervolume values and print them on the screen
    index = sum(sample_g <= 0, 2) == num_con;
    if sum(index) == 0
        hypervolume(1) = 1E6;
        fprintf(' iteration: %d, evaluation: %d, hypervolume: no feasible solution found\n', iteration, evaluation);
    else
        feasible_y = sample_y(index, :);
        feasible_y_scaled = sample_y_scaled(index, :);
        index = Paretoset(feasible_y);
        non_dominated_front = feasible_y(index,:);
        non_dominated_front_scaled = feasible_y_scaled(index,:);
        hypervolume(iteration+1) = Hypervolume(non_dominated_front,ref_point);
        % plot current non-dominated front points
        if num_obj == 2
            scatter(non_dominated_front(:,1), non_dominated_front(:,2),'ro', 'filled');
            title(sprintf('%s problem, iteration: %d, evaluations: %d',fun_name, iteration,evaluation));drawnow;
        elseif num_obj == 3
            scatter3(non_dominated_front(:,1), non_dominated_front(:,2),non_dominated_front(:,3),'ro', 'filled');
            title(sprintf('%s problem,iteration: %d, evaluations: %d',fun_name,iteration,evaluation));drawnow;
        end
        fprintf(' iteration: %d, evaluation: %d, hypervolume: %f\n',  iteration, evaluation, hypervolume(iteration +1));
    end
end



