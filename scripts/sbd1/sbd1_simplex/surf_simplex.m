%% Plot phi
surf(gsamps, gsamps, phi_g);  hold on;  
colormap parula;  shading interp;  

if plotcontour
    contour(gsamps, gsamps, phi_g);  colormap parula
end

a = [d C(:,1)+d C(:,2)+d d];    % show the simplex
plot(a(1,:), a(2,:), 'k', 'LineWidth', 1);
plot(a(1,:), a(2,:), 'ro', 'LineWidth', 1.2, 'MarkerSize', 7);

for i = 1:3
    text(a(1,i), a(2,i), sprintf('  S_{%d}[a0]', s(i)));
end

a = C*[1 1]'/3 + d;             % middle point
plot(a(1), a(2), 'bx', 'LineWidth', 1.2, 'MarkerSize', 7);
text(a(1), a(2), '  center');

view(180,-90);  drawnow;
xlabel('g_1');  ylabel('g_2');