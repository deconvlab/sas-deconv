%% REWEIGHT_LOOP  Run inner and outer iterations with reweighting
%
%   PARAMETERS:
%   ===========
%
%   solver:  iPALM iterator for the original problem solve.
%
%   reweights:  number of reweighting loops to perform. Can be set to 0.
%
%   lambda2:  sparsity-tradeoff parameter for reweighting.
%
%   eps:  reweighting parameter.
%   
%   update_script:  path to a script for displaying updates. By default
%     text updates are displayed. Setting to empty displays no update.
%
%   maxit:  an array with REWEIGHTS+1 positive integers. Each element
%     provides the maximum iterations to be performed by the solver. The
%     first element corresponds to the initial solve (no reweighting), etc.
%
%   updates: an array of positive integers providing the iteration numbers 
%     in each inner loop to display updates.
%
%
%   VARIABLES CREATED:
%   ==================
%   It is useful to keep track of some of the variables created within the
%   loop for writing UPDATE_SCRIPTs, etc.
%
%   rwgt:  the current reweighting loop. The initial solve (no reweighting)
%     happens for RWGT = 1.
%
%   it:  inner loop iteration.
%
%   solvers:  a cell array containing solvers copied over at each 
%     reweighting loop. Probably a bad idea for large problem sizes.
%
%   costs:  a cell array for each reweighting loop containing costs 
%     generated at inner loop iterations.
%
%   stime:  starting TIC time of the loop.
%
%


%% Script

if ~exist('update_script', 'var')
    update_script = 'text_update.m';
end

solvers = cell(reweights+1,1);
costs = cell(reweights+1,1);  
stime = tic;  
title('Original observation');
for rwgt = 1:reweights+1
    if rwgt > 1
        solver = reweight(solver, lambda2, eps);
    else
        figure(2); clf;
    end
    costs{rwgt} = NaN(maxit(rwgt),1);
    
    for it = 1:maxit(rwgt)
        solver = iterate(solver);
        costs{rwgt}(it) = solver.cost;
        
        if ismember(it, updates) && ~isempty(update_script)
            run(update_script);
        end
    end
    disp(' ');
    solvers{rwgt} = copy(solver);
end