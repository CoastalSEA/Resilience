ispos = true;  %logical flag to invert the y values
x = (0:0.01:2);

A = 0;      %left horizontal asymptote
B = 10;     %growth rate
C = 1;      %upper asymptote of A+(K-A)/C^(1/nu) - typically 1
K = 1;      %right horizontal asymptote when C=1
nu = 1;     %>0 - affects proximity to which asymptote maximum growth occurs
M = 0.7;    %start point
Q = 2;      %related to the value y(0).
params = struct('A',A,'B',B,'C',C,'K',K,'nu',nu,'M',M,'Q',Q);

n = [.2,.5,1,2,5];
testparam = 'nu';

figure('Tag','PlotFig');
axes;
hold on
if ispos
    plot([0,0.8,1,2],[0,0.5,0.8,1],'DisplayName','Linear segments')
else
    plot([0,0.8,1,4],[1,0.8,0.2,0],'DisplayName','Linear segments')
end
for i=1:5
    params.(testparam) = n(i);
    y = general_logistic(x,params,~ispos);
    plot(x,y,'DisplayName',num2str(n(i)))
end
legend
title(testparam)
hold off