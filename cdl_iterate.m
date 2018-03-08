solvers = cell(reweights+1,1);
costs = cell(reweights+1,1);  
stime = tic;  %profile on;
figure(1);  subplot(3,2,[1 3]);  imagesc(abs(Y-b0));
for r = 1:reweights+1
    if r > 1
        solver = reweight(solver, lambda2, eps);
    else
        figure(2); clf;
    end
    costs{r} = NaN(maxit(r),1);
    
    for i = 1:maxit(r)
        solver = iterate(solver);
        costs{r}(i) = solver.cost;

        if centerfq > 0 && mod(i, centerfq) == 0
            for k = 1:K
                [A, X, ~, A_, X_] = center(...
                    solver.A{k}, solver.X{k}, [],...
                    solver.A_{k}, solver.X_{k});
                solver.A{k} = A;  solver.X{k} = X;
                solver.A_{k} = A;  solver.X_{k} = X;
            end
        end
        
        if ismember(i, updates)
            tmp = ['Iteration ' num2str(solver.it)];
            
            figure(1);
            subplot(3,2,[2 4]); 
                imagesc(abs(synthesize(solver.A, solver.X, {0})));
            subplot(313); 
            for k = 1:r
                c = costs{k};
                if k == 1;  c = c(50:end);  end
                plot(linspace(0,1,numel(c)), c-min(c));  hold on;
            end
            hold off; title(tmp); drawnow;
            
            figure(2);
            for k = 1:K
                subplot(2,2*K,2*K*(r>1) + 2*k-1); 
                imagesc(abs(solver.A{k}));
                subplot(2,2*K,2*K*(r>1) + 2*k);
                imagesc(solver.X{k});
            end
            subplot(2,2*K,2*K*(r>1)+1); title(tmp);
            drawnow;
            
            fprintf(['Iter. %d:%d.  '...
                'Cost %.4e. Elapsed time %.2fs.\n'], ...
                r-1, solver.it, solver.cost, toc(stime));
        end
    end
    disp(' ');
    solvers{r} = copy(solver);
end
%profile off; profile viewer;