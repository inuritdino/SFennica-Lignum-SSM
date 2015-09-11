function lpfg_lgm_1par_run(parname,values,varargin)

%%
argsConst = {''};
ConstValues = [];
tf = strcmp('argsConst',varargin);
if(find(tf))
    argsConst = varargin{find(tf)+1};
end
tf = strcmp('C',varargin);
if(find(tf))
    ConstValues = varargin{find(tf)+1};
end
%% Init and error check


%% Change directory and remember the current one
curr_dir = pwd;
cd /Users/potapov/docs/programs/vlab/lignum;

%% Run the simulations

for ii = 1:length(values)
    optim_lgm_call(values(ii),'args',{[parname '__1']},'visual',...
        'C',ConstValues,'argsConst',argsConst);
end

%% Return back
cd(curr_dir);

end