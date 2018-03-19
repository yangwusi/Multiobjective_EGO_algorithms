% -----------------------------------------------------------------------------------------
% The parallel multiobjective EGO algorithm using PEIM(Pseudo Expected Improvement
% Matrix) criteria, which is significant cheaper-to-evaluate than the
% state-of-the-art multiobjective EI criteria. For detailed description
% about the PEIM and EIM criteria, please refer to [1].
% The dace toolbox [2] is used for building the Kriging models in the 
% implementations.
% The non-dominated sorting method by Yi Cao [3] is used to identify the
% non-dominated fronts from all the design points
% The hypervolume indicators are calculated using the faster algorithm of
% [4] Nicola Beume et al. (2009).
% -----------------------------------------------------------------------------------------
% [1]  Dawei Zhan, Yuansheng Cheng, Jun Liu, Expected Improvement Matrix-based Infill 
%      Criteria for Expensive Multiobjective Optimization. IEEE Transactions 
%      on Evolutionary Computation, 2017, 21 (6): 956-975.
% [2] Dawei Zhan, Jiachang Qian, Jun Liu, et al. Pseudo Expected Improvement Matrix Criteria for 
%     Parallel Expensive Multi-objective Optimization. In Advances in Structural and Multidisciplinary
%     Optimization: Proceedings of the 12th World Congress of Structural and Multidisciplinary 
%     Optimization (WCSMO12), Schumacher, A.,Vietor, T.,Fiebig, S., et al., Eds. Springer International
%     Publishing: Cham, 2018; 175-190.
% [3] Lophaven SN, Nielsen HB, and Sodergaard J, DACE - A MATLAB Kriging 
%     Toolbox, Technical Report IMM-TR-2002-12, Informatics and Mathematical
%     Modelling, Technical University of Denmark, 2002. Available at:
%      http://www2.imm.dtu.dk/~hbn/dace/.
% [4] http://www.mathworks.com/matlabcentral/fileexchange/17251-
%      pareto-front.
% [5] N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold, 
%     On the Complexity of Computing the Hypervolume Indicator, IEEE 
%     Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
% -----------------------------------------------------------------------------------------
% zhandawei@hust{dot}edu{dot}cn
% 2018.03.19 initial creation
% -----------------------------------------------------------------------------------------
clearvars;close all;
addpath('dace', 'Infill_Criterion', 'Test_Problem');
%-------------------------------------------------------------------------
% settings of the problem
% infill criterion: 'PEIM_Euclidean','PEIM_Maximin','PEIM_Hypervolume'
infill_name = 'PEIM_Hypervolume';
% number of updating points selected in each cycle
num_q = 5;
% for ZDT test problems, the number of objectives should be 2
test_name='DTLZ7';
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
fprintf('Number of updating points   : %d\n',num_q)
fprintf('Infill Criterion            : %s\n','Euclidean distance PEIM crtierion')
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
    % select q updating points use PEIM criterion
    best_x = zeros(num_q, num_vari);
    point_added = [];
    for ii = 1 : num_q
        % find the maximum of the pseudo EI function
        switch infill_name
            case 'PEIM_Euclidean'
                infill_criterion = @(x)Infill_Pseudo_EIM_Euclidean(x, kriging_obj, non_dominated_front_scaled, point_added);
            case 'PEIM_Maximin'
                infill_criterion = @(x)Infill_Pseudo_EIM_Maximin(x, kriging_obj, non_dominated_front_scaled, point_added);
            case 'PEIM_Hypervolume'
                infill_criterion = @(x)Infill_Pseudo_EIM_Hypervolume(x, kriging_obj, non_dominated_front_scaled, point_added);
            otherwise
                error('you should select infill_name from PEIM_Euclidean, PEIM_Maximin, and PEIM_Hypervolume');
        end
        best_x(ii, :) =  particleswarm(infill_criterion, num_vari,design_space(1,:),design_space(2,:),options);
        % check if the candidate point is too close to sampled points
        if min(sqrt(sum((best_x(ii, :)-[sample_x;point_added]).^2,2)))<1E-8
            best_x(ii, :)=particleswarm(@(x)Infill_Maximal_Distance(x, [sample_x;point_added]), num_vari, design_space(1,:), design_space(2,:), options);
        end
        % update point_added
        point_added = [point_added; best_x(ii, :)];
    end
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

