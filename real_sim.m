%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The sample simulation code to run the real case reported in the paper
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
% additional figures are shown during the optimization: 1)the cumulative
% contribution of parameters to the fitness (structural distance) at each step and
% 2) the best and average fitness along the generations and the average 
% (genetic) distance between population's individuals.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add function folders to the MATLAB path
curr_dir = pwd;
addpath(genpath(curr_dir));

%% Enter the LPFG-LIGNUM directory
cd lpfg-lignum/;

%% Extract the DATA from the structural data files in the QSM library
lib_fn = 'pine_Ruotsinkyla_tree3_model1.mat';% library file name
% Generate all possible structural features and the data tree
[s, data_trunk, data_tree] = gen_scatter(['../qsm-library/' lib_fn]);
% Extract only the target u-spaces
data_u = extract_scatter(s,{'taper','curv','bra'},1);
% Plot the data tree
figure(22); data_tree = data_tree.move_tree([0 0 0]); data_tree.draw;

%% Define the trial model to fit to the DATA
% Parameters of the LPFG-LIGNUM to vary, see lgmconst.h file at 
% lpfg-lignum/ for the full description (_GAUSS indicates a random parameter
% following gaussing distribution):
% LR_GAUSS, Q_GAUSS, YEARS, SHLEN, SH(30), ENVSHRAD_GAUSS, BRAANG, 
% ANGINCR, ZETASD, MAXANG, GAMMASD, SHEDYEAR, SHEDMU.
model_u = @(x)optim_lgm_call(x,'args',{'LR_GAUSS_1_2','Q_GAUSS_3_4','YEARS__5',...
    'SHLEN__6','SH__7','ENVSHRAD_GAUSS_8_9','BRAANG__10','ANGINCR__11',...
    'ZETASD__12','MAXANG__13','GAMMASD__14','SHEDYEAR__15','SHEDMU__16'},...
    'C',[0 1 1 0.02 500 60.0 2],'argsConst',{'VERB__1','RNDSEED__2',...
    'SIEVANEN2008__3','VSA__4','GSA__5','QRAD__6','SHEDDIST__7'},...
    'scat',{'taper','curv','bra'},'order',1);

%% Technical parameters to the genetic algorithm
% initial range of the parameters to vary (here coincides with the global
% boundaries)
InitRange = [];
% Global LOW boundary for the parameter values
LB = [0.009 0.001 0.1 0    20 0.6 30 0.1 0.05 20 0  0  40 0  5  0.0];
% Global UP boundary for the parameter values (+50% of the target values)
UB = [0.011 0.003 0.4 0.05 40 1.5 80 1.5 0.3  50 15 10 80 10 15 6.0];

% More specific parameters of the genetic algorithm (see elsewhere, or
% MATLAB documentation for GA)
INTCON = [5 15];% the 5th and 15th parameters are integers
POPSIZE = 1;% population size, 50
ELITE = 0;% elite size, 2
GENS = 3;% number of generations, 50
STALL = 1;% stall generations limit,5

%% Call optimization routine (genetic algorithm)
DIM = length(LB);% Dimension of the problem (number of parameters to optimize)
% Call the optimization routine resulting in the best solution X, i.e. 
% optimal values of the 16 parameters.
% OPTIM_CALL here takes additional args 'stat' (the type of 1D statistic to
% compute, see dt_distance()) and 'dirs' (number of random directions to
% generate distance).
X = optim_call(data_u,model_u,'nvars',DIM,...
    'opts',{'PlotFcns',{@gaplotbestf,@gaplotdistance},'PopulationSize',POPSIZE,...
    'EliteCount',ELITE,'Generations',GENS,'StallGenLimit',STALL,...
    'PopInitRange',InitRange,'OutputFcns',@optim_lgm_output},...
    'lb',LB,'ub',UB,'intcon',INTCON,'stat',1,'dirs',1000);
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
subplot(1,4,1);
plot(data_u{1}(1,:),data_u{1}(2,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{1}(1,:),optim_model_u{1}(2,:),'or','MarkerSize',8,'LineWidth',2);
hold off
legend('DATA','MODEL');
title('Tapering function');
xlabel('Branch length, m');
ylabel('Radius, m');
subplot(1,4,2);
plot(data_u{2}(3,:),(180/pi)*data_u{2}(1,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{2}(3,:),(180/pi)*optim_model_u{2}(1,:),'or','MarkerSize',8,'LineWidth',2);
hold off;
title('Spatial curvature');
xlabel('Relative branch length');
ylabel('Horizontal plane angle, deg');
subplot(1,4,3);
plot(data_u{2}(3,:),(180/pi)*data_u{2}(2,:),'ob','MarkerSize',8,'LineWidth',2);
hold on;
plot(optim_model_u{2}(3,:),(180/pi)*optim_model_u{2}(2,:),'or','MarkerSize',8,'LineWidth',2);
hold off;
title('Spatial curvature');
xlabel('Relative branch length');
ylabel('Vertical plane angle, deg');
subplot(2,4,4);
[n1,x1] = hist((180/pi)*data_u{3}); bar(x1,n1,'b');
subplot(2,4,8);
[n2,x2] = hist((180/pi)*optim_model_u{3}); bar(x2,n2,'r');
xlabel('Branching angle, deg');
low = min(min(x1),min(x2)) - max(x1(2)-x1(1),x2(2)-x2(1));
high = max(max(x1),max(x2)) + max(x1(2)-x1(1),x2(2)-x2(1));
subplot(2,4,4);xlim([low high]);
subplot(2,4,8);xlim([low high]);
