%% CDL_UPDATE

%% Centering + plotting + text
K = numel(solver.A);
tmp = ['Iteration ' num2str(solver.it)];

if centerfq > 0 && mod(it, centerfq) == 0
    [solver.A, solver.X, solver.A_, solver.X_] = ...
        center(solver.A, solver.X, solver.A_, solver.X_);
end

figure(1);
subplot(3,2,[2 4]); 
imagesc(abs(synthesize(solver.A, solver.X, {0})));
title('Reconstruction');

subplot(313); 
for k = 1:rwgt
    c = costs{k};
    if k == 1;  c = c(50:end);  end
    plot(linspace(0,1,numel(c)), c-min(c));  hold on;
end
title(tmp); xlabel('Cost'); ylabel('Iteration (fraction)');
hold off; 
drawnow;

figure(2);
for k = 1:K
    subplot(2,2*K,2*K*(rwgt>1) + 2*k-1); 
    imagesc(abs(solver.A{k}));
    xlabel(['Kernel ' num2str(k)]);

    if k == 1
        if rwgt == 1
            ylabel('Before reweighting');
        else
            ylabel('After reweighting');
        end
    end

    subplot(2,2*K,2*K*(rwgt>1) + 2*k);
    imagesc(solver.X{k});
    xlabel(['Activation ' num2str(k)]);
end
subplot(2,2*K,2*K*(rwgt>1)+1); title(tmp);
drawnow;

text_update;