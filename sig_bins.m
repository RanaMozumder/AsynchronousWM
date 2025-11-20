function h_fdr = sig_bins(baseline_rate, delay_rate)
for bin = 1:size(delay_rate, 2)
    [~, p(bin)] = ttest(delay_rate(:, bin), baseline_rate);
end

% Apply FDR correction
[h_fdr, crit_p, adj_p] = fdr_bh(p, 0.05, 'pdep', 'yes');
end

function [h, crit_p, adj_p] = fdr_bh(pvals, q, method, report)
% FDR_BH   Benjamini & Hochberg (1995) False Discovery Rate control.
%
% [h, crit_p, adj_p] = fdr_bh(pvals, q, method, report)
%
% Input:
% - pvals: vector of p-values
% - q: FDR threshold (e.g., 0.05)
% - method: 'pdep' (independent) or 'dep' (dependent; more conservative)
% - report: 'yes' or 'no' to print summary
%
% Output:
% - h: 1 if significant after FDR, 0 otherwise
% - crit_p: p-value threshold for significance
% - adj_p: adjusted p-values

if nargin < 2 || isempty(q), q = 0.05; end
if nargin < 3, method = 'pdep'; end
if nargin < 4, report = 'no'; end

pvals = pvals(:);
V = length(pvals);
[sorted_pvals, sort_ids] = sort(pvals);
unsort_ids(sort_ids) = 1:V;

if strcmpi(method, 'pdep')
    thresh = (1:V)'/V*q;
elseif strcmpi(method, 'dep')
    cV = sum(1./(1:V));
    thresh = (1:V)'/(V*cV)*q;
else
    error('Method must be "pdep" or "dep".');
end

w = find(sorted_pvals <= thresh);
if isempty(w)
    crit_p = 0;
    h = zeros(V,1);
else
    crit_p = sorted_pvals(max(w));
    h = pvals <= crit_p;
end

adj_p = zeros(V,1);
min_adj_p = 1;
for i = V:-1:1
    min_adj_p = min(min_adj_p, sorted_pvals(i)*V/i);
    adj_p(sort_ids(i)) = min_adj_p;
end
adj_p = min(adj_p,1);

h = h(:);
adj_p = adj_p(:);

if strcmpi(report, 'yes')
    fprintf('FDR threshold q = %.4f\n', q);
    fprintf('%d/%d tests are significant after FDR correction.\n', sum(h), V);
end
end
