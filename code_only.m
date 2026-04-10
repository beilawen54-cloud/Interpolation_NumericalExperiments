%% Phase 3 - Interpolation
% 1.Sample Points

clear; clc; format long;

alpha = 0.25;
f = @(v)9.8 - 0.5*v.^2.*(1+alpha*log(1+v));

v_nodes = [3 3.5 4 4.5 5];
f_nodes = f(v_nodes);

f_nodes
% 2.Newton Interpolating Polynomial

n = length(v_nodes);
DD = zeros(n,n);
DD(:,1) = f_nodes';

for j = 2:n
    for i = 1:n-j+1
        DD(i,j) = (DD(i+1,j-1) - DD(i,j-1)) / (v_nodes(i+j-1) - v_nodes(i));
    end
end

a = DD(1,:);
a
% 3.Solving p(v)=0

function val = newton_eval(x, nodes, a)
    n = length(a);
    val = a(n);
    for k = n-1:-1:1
        val = val.*(x - nodes(k)) + a(k);
    end
end

p = @(v) newton_eval(v, v_nodes, a);

v_interp = fzero(p, [3,5]);
v_true = fzero(f, [3,5]);

v_interp
v_true
%% Phase 5 - Numerical Experiments
% 1.Bisection model setup

function [v_approx, n] = bisection_terminal(alpha, tol, maxIter)

f = @(v) 9.8 - 0.5*v.^2.*(1 + alpha*log(1+v));
a = 3; b = 5;

if f(a)*f(b) >= 0
    error('Method failed: f(a) and f(b) do not have opposite signs.')
end
converged = false;
for n = 1:maxIter
    c = (a + b)/2;

    if abs(f(c)) < tol || (b-a)/2 < tol
        converged = true;
        break
    end

    if f(a)*f(c) < 0
        b = c;
    else
        a = c;
    end
end

if ~converged
    error('Method failed: max iterations exceeded')
end

v_approx = (a + b)/2;

end
% 2.Different alpha values

clear; clc; format long;

alphas = [0, 0.25, 0.5];
tol = 1e-8;
maxIter = 200;

for j = 1:length(alphas)
    alpha = alphas(j);

    [v_approx, n] = bisection_terminal(alpha, tol, maxIter);

    fprintf('alpha = %.2f\n', alpha);
    fprintf('terminal velocity = %.10f\n', v_approx);
    fprintf('iterations used = %d\n\n', n);
end