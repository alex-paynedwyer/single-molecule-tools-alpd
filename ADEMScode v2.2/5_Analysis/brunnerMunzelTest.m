%values1 = [comma separated list of values for condition 1];
%values2 = [comma separated list of values for condition 2];

function BrunnerMunzelTest(values1, values2, alpha, numtests)

if min(length(values1),length(values2))<21
    disp("Too few values for a reliable comparison under this implementation of the Brunner-Munzel test.")
    disp("Try again with more data (>20).")
else
if nargin<4
    alpha = 0.05;
    numtests = 5; %the number of truly independent tests for one set of data:
    %e.g. stoichiometry, periodicity, diffusivity, linked v. unlinked. 
    % Five is usually a safe bet for our work.
end

reduced_alpha = alpha/numtests; %apply Bonferroni correction
% to avoid spurious 'significance' across multiple comparisons

[p_value] = brunnerMunzel(values1, values2); %compute the Brunner-Munzel metric from Max Becker's code

if p_value < reduced_alpha
    disp(strcat("significant difference in distributions under BM test: p =", num2str(p_value), " < ", num2str(reduced_alpha)))
else
    disp(strcat("not significant (distributions can't be distinguished) under BM test: p =", num2str(p_value), " > ", num2str(reduced_alpha)))
end
end