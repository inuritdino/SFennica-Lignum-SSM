%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The sample simulation code to run the test case reported in the paper
% "Data-based stochastic modeling of tree growth and structure formation"
% in Silva Fennica, 2015 (anonymous).
%
% REQUIREMENTS:
% 1. Matlab, The MathWorks.
% 2. LPFG simulator. The simulator is available as a part of L-studio (Windows)
% or VLab (Unix/Mac) simulation environments. Both can be downloaded at 
% http://algorithmicbotany.org/virtual_laboratory/ along with the 
% documentation and detailed installation guides.
%
% CODE:
% 1. Adds paths to Matlab to find the functions.
% 2. Generates the 'DATA' from the model (stochastic Lignum) itself.
% 3. Sets up the trial model.
% 4. Sets up the technical parameters to the genetic algorithm.
% 5. Runs the optimization routine.
% 6. Optionally, plots the u-space of DATA and MODEL.
% 
% As a result, one gets two trees in figure 22 of Matlab: the data-tree
% (origin at [0 0 0]) and the optimized tree (origin at [5 0 0]). Two
% additional figures are shown during the optimization: one with cumulative
% contribution of parameters to the fitness (structural distance) at each step and
% another with the best and average fitness along the generations
% and the average (genetic) distance between population's individuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add function folders to the MATLAB path
curr_dir = pwd;
addpath([curr_dir '/lpfg-lignum-matlab'],genpath([curr_dir '/optim-tool-matlab']));

%% Enter the LPFG-LIGNUM directory
cd lpfg-lignum/;

%% Generate 'DATA' for the test case
% Parameters of the LPFG-LIGNUM, see lgmconst.h file at lpfg-lignum/ for the description
% LR_GAUSS(0.009,0.001), Q_GAUSS(0.2,0.03), YEARS(15), ANGINCR(10),
% ZETASD(5), BRAANG(35)
% Call the model and generate the data u-spaces.
X0 = [0.009 0.001 0.2 0.03 15 10.0 5 35];
PAR = {'LR_GAUSS_M','LR_GAUSS_S','Q_GAUSS_M','Q_GAUSS_S','YEARS','ANGINCR',...
    'ZETASD','BRAANG'};
data_u = optim_lgm_call(X0,'args',...
    {'LR_GAUSS_1_2','Q_GAUSS_3_4','YEARS__5','ANGINCR__6','ZETASD__7','BRAANG__8'},...
    'C',[1 1 500 0.02],'argsConst',{'PERTTUNEN1998__1','RNDSEED__2','GSA__3','VSA__4'},...
    'scat',{'taper','curv'},'order',1);
% Read the structure into the MATLAB representation and then plot the data tree.
data_tree = read_mtg('lignum.mtg');
figure(22); data_tree.draw;

%% Define a trial model to fit to the DATA
% LPFG-LIGNUM call with 8 parameters as input
model_u = @(x)optim_lgm_call(x,'args',...
    {'LR_GAUSS_1_2','Q_GAUSS_3_4','YEARS__5','ANGINCR__6','ZETASD__7','BRAANG__8'},...
    'C',[1 2 500 0 0.02],...
    'argsConst',{'PERTTUNEN1998__1','RNDSEED__2','GSA__3','VERB__4','VSA__5'},...
    'scat',{'taper','curv'},'order',1);

%% Technical parameters to the genetic algorithm
% initial range of the parameters to vary (here coincides with the global
% boundaries)
InitRange = [0.0045 0.0005 0.1 0.015 7  5  2.5 17.5;
             0.0135 0.0015 0.3 0.045 23 15 7.5 52.5];
% Global LOW boundary for the parameter values (-50% of the target values)
LB = [0.0045 0.0005 0.1 0.015 7  5  2.5 17.5];
% Global UP boundary for the parameter values (+50% of the target values)
UB = [0.0135 0.0015 0.3 0.045 23 15 7.5 52.5];

% More specific parameters of the genetic algorithm (see elsewhere, or
% MATLAB documentation for GA)
INTCON = 5;% the 5th parameter is integer
POPSIZE = 40;% population size, 40
ELITE = 4;% elite size, 4
GENS = 50;% number of generations, 50
STALL = 5;% stall generations limit,5

%% Call optimization (genetic algorithm)
DIM = length(LB);% Dimension of the problem (number of parameters to optimize)
% Call the optimization routine resulting in the best solution X, i.e. 
% optimal values of the 8 parameters.
X = optim_call(data_u,model_u,'nvars',DIM,...
    'opts',{'PlotFcns',{@gaplotbestf,@gaplotdistance},'PopulationSize',POPSIZE,...
    'EliteCount',ELITE,'Generations',GENS,'StallGenLimit',STALL,...
    'PopInitRange',InitRange,'OutputFcns',@optim_lgm_output},...
    'lb',LB,'ub',UB,'intcon',INTCON);
% Generate the final u-space set and the tree with the optimal solution X
optim_model_u = model_u(X);
% Read the resulting tree structure
model_tree = read_mtg('lignum.mtg');

% Plot the tree along with the data-tree
figure(22); hold on;
model_tree = model_tree.move_tree([5 0 0]);% move tree in space
model_tree.draw;
hold off;

%% Optional plot of u-spaces
h = figure;
set(h,'Position',[0 0 1000 350]);
subplot(1,3,1);
plot(data_u{1}(1,:),data_u{1}(2,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{1}(1,:),optim_model_u{1}(2,:),'or','MarkerSize',8,'LineWidth',2);
hold off
legend('DATA','MODEL');
title('Tapering function');
xlabel('Branch length, m');
ylabel('Radius, m');
subplot(1,3,2);
plot(data_u{2}(3,:),(180/pi)*data_u{2}(1,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{2}(3,:),(180/pi)*optim_model_u{2}(1,:),'or','MarkerSize',8,'LineWidth',2);
hold off;
title('Spatial curvature');
xlabel('Relative branch length');
ylabel('Horizontal plane angle, deg');
subplot(1,3,3);
plot(data_u{2}(3,:),(180/pi)*data_u{2}(2,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{2}(3,:),(180/pi)*optim_model_u{2}(2,:),'or','MarkerSize',8,'LineWidth',2);
hold off;
title('Spatial curvature');
xlabel('Relative branch length');
ylabel('Vertical plane angle, deg');

%% Error reports
Err = 100.*abs(X0-X)./X0;
fprintf('=== ERRORS ===\n');
fprintf('Parameter,\t Error\n');
for ii = 1:DIM
    fprintf('%s,\t %g %%\n',PAR{ii},Err(ii));
end
fprintf('==============\n');
