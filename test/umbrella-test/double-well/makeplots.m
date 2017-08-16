clear; clc; clf;

xrec_nobias = dlmread('nobias.dat');
[cts, edges] = histcounts(xrec_nobias,20);
p_nobias = cts/sum(cts);
ctrs = edges(1:end-1) + 0.5*(edges(2)-edges(1));
plot(ctrs, p_nobias, 'LineWidth', 2)

% xrec_bias = dlmread('bias.dat');
% [cts, edges] = histcounts(xrec_bias,2000);
% p_bias = cts/sum(cts);
% ctrs = edges(1:end-1) + 0.5*(edges(2)-edges(1));
% plot(ctrs, p_bias, 'LineWidth', 2)
% ub = @(x) -(abs(x)<1).*(-(1/2)*x.^2 + (1/4)*x.^4 + (1/4));
% p_unbiased = cts.*exp(ub(ctrs)) / sum(cts.*exp(ub(ctrs)));
% plot(ctrs, p_unbiased, 'LineWidth', 2)

set(gca, 'Fontsize', 14)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$\log p(x)$', 'Interpreter', 'latex')


%%
hold on
x = 1:0.1:10;
semilogy(x, exp(-x))
% xlim([0 1])
%%

clf

utot = @(x) (-(1/2)*x.^2 + (1/4)*x.^4)-(abs(x)<1).*(-(1/2)*x.^2 + (1/4)*x.^4 + (1/4));
plot(ctrs, utot(ctrs))
