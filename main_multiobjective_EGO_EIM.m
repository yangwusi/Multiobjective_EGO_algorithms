% -----------------------------------------------------------------------------------------
% The multiobjective EGO algorithm using EIM(expected improvement
% matrix)-based criteria, which is significant cheaper-to-evaluate than the
% state-of-the-art multiobjective EI criteria. For detailed description
% about the EIM criteria, please refer to [1].
% The dace toolbox [2] is used for building the Kriging models in the
% implementations.
% The non-dominated sorting method by Yi Cao [3] is used to identify the
% non-dominated fronts from all the design points
% The hypervolume indicators are calculated using the faster algorithm of
% [4] Nicola Beume et al. (2009).
% -----------------------------------------------------------------------------------------
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
% -----------------------------------------------------------------------------------------
% zhandawei@hust{dot}edu{dot}cn
% 2017.05.03 initial creation
% 2018.03.19 update
% -----------------------------------------------------------------------------------------
clearvars;close all;
addpath('dace', 'Infill_Criterion', 'Test_Problem');
%-------------------------------------------------------------------------
% settings of the problem
% infill criterion: 'EIM_Euclidean','EIM_Maximin','EIM_Hypervolume'
infill_name = 'EIM_Euclidean';
% for ZDT test problems, the number of objectives should be 2
test_name='DTLZ2';
% number of objectives
num_obj = 3;
% number of design variables
num_vari = 6;
% number of initial design points
num_initial_sample = 60;
% the maximum allowed iterations
max_iteration = 20;
%-------------------------------------------------------------------------
% get the information about the problem
[design_space, ref_point]=Test_Function(test_name, num_obj, num_vari);
% the intial design points, points sampled all at once
sample_x = design_space(1,:)+(design_space(2,:)-design_space(1,:)).*lhsdesign(num_initial_sample,num_vari,'criterion','maximin','iterations',1000);
sample_y=feval(test_name, sample_x, num_obj);
% scale the objectives to [0,1]
sample_y_scaled=(sample_y-min(sample_y))./(max(sample_y)-min(sample_y));
%-------------------------------------------------------------------------
% initialize some parameters
evaluation = size(sample_x,1);
kriging_obj=cell(1,num_obj);
hypervolume=zeros(max_iteration+1,1);
% the settings of PSO optimizer
options=optimoptions('particleswarm','SwarmSize',100,'MaxIterations',100,'MaxStallIterations',100,'Display','off', 'UseVectorized', true);
%-------------------------------------------------------------------------
% calculate the initial hypervolume values and print them on the screen
indx=Paretoset(sample_y);
non_dominated_front=sample_y(indx,:);
non_dominated_front_scaled=sample_y_scaled(indx,:);
hypervolume(1)=Hypervolume(non_dominated_front,ref_point);
fprintf('----------------------------------------------------------------\n')
fprintf('Test Function               : %s\n',test_name)
fprintf('Number of objectives        : %d\n',num_obj)
fprintf('Number of variables         : %d\n',num_vari)
fprintf('Infill Criterion            : %s\n','Euclidean distance EIM crtierion')
fprintf('----------------------------------------------------------------\n')
fprintf(' iteration: %d, evaluation: %d, hypervolume: %f \n', 0, evaluation, hypervolume(1));
%-------------------------------------------------------------------------
% beginning of the iteration
for iter= 1 : max_iteration
    %-------------------------------------------------------------------------
    % build the initial kriging model for each objective
    for ii=1:num_obj
        kriging_obj{1,ii}=dacefit(sample_x,sample_y_scaled(:,ii),'regpoly0','corrgauss',1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    %-------------------------------------------------------------------------
    % select updating points using the EIM criteria
    switch infill_name
        case 'EIM_Euclidean'
            infill_criterion = @(x)Infill_Standard_Euclidean_EIM(x, kriging_obj, non_dominated_front_scaled);
        case 'EIM_Maximin'
            infill_criterion = @(x)Infill_Standard_Maximin_EIM(x, kriging_obj, non_dominated_front_scaled);
        case 'EIM_Hypervolume'
            infill_criterion = @(x)Infill_Standard_Hypervolume_EIM(x, kriging_obj, non_dominated_front_scaled);
        otherwise
            error('you should select infill_name from EIM_Euclidean, EIM_Maximin, and EIM_Hypervolume');
    end
    [best_x, best_EIM]=particleswarm(infill_criterion,num_vari,design_space(1,:),design_space(2,:),options);
    %-------------------------------------------------------------------------
    % add the new points to the design set
    sample_x=[sample_x;best_x];
    sample_y=[sample_y; feval(test_name,best_x, num_obj)];
    sample_y_scaled=(sample_y-min(sample_y))./(max(sample_y)-min(sample_y));
    evaluation=evaluation+size(best_x,1);
    %-------------------------------------------------------------------------
    % calculate the hypervolume values and print them on the screen
    indx=Paretoset(sample_y);
    non_dominated_front=sample_y(indx,:);
    non_dominated_front_scaled=sample_y_scaled(indx,:);
    hypervolume(iter+1)=Hypervolume(non_dominated_front,ref_point);
    fprintf(' iteration: %d, evaluation: %d, hypervolume: %f\n', iter, evaluation, hypervolume(iter +1));
end
%-------------------------------------------------------------------------
% plot the final approximated pareto front is the number of objective is 2
% or 3
figure;
if num_obj == 2
    scatter(non_dominated_front(:,1), non_dominated_front(:,2),'ro', 'filled')
elseif num_obj == 3
    scatter3(non_dominated_front(:,1), non_dominated_front(:,2),non_dominated_front(:,3),'ro', 'filled')
end
%-------------------------------------------------------------------------

