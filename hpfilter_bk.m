
function [tr, cyc] = hpfilter_bk(y,w,graph)

% Applies the Hodrick-Prescott filter.
% Usage:   [tr, cyc] = hpfilter_bk(y,w,graph)
% Inputs:  y:     Original series
%          w:     Smoothing parameter; default is 1600 (choose 100 for
%                 annual, 1600 for quarterly and 14400 for monthly data)
%          graph: (optional) If graph=1, the results will be plotted.
% Outputs: tr:  Trend component of filtered series
%          cyc: Cyclical component of filtered series
% by Benedikt Kolb (21/01/14), based on code by Wilmer Henao
%
% For a derivation of the formulas used, see e.g.
% http://economics.sas.upenn.edu/~jesusfv/lecturetechnical5_spectral.pdf

if nargin == 1
    warning('hp:no_w','No smoothing parameter w specified. I assume you want to use w=1600 in the following.');
    w = 1600;
end

if size(y,1) < size(y,2)
    y = y';
end
[T N] = size(y);

d = repmat([w -4*w (6*w+1)/2], T, 1);
d(1,2) = -2*w;      d(T-1,2) = -2*w;
d(1,3) = (1+w)/2;   d(T,3) = (1+w)/2;
d(2,3) = (5*w+1)/2; d(T-1,3) = (5*w+1)/2;
B = spdiags(d, -2:0, T, T);    % sparse version of B (due to W. Henao)
B = B+B';
tr  = B\y;
cyc = y - tr;

%% (optional) Plot BP-filtered vs original series
if nargin == 3 && graph == 1
    figure('Name','BP-filtered vs original series')
    for i = 1:N
        subplot(max(floor(N/3),1),ceil(N/max(floor(N/3),1)),i);
        plot(1:T,tr(:,i),'r',1:T,y(:,i),'k--'); grid on; title(['Series #',num2str(i)]);
        if i==1, legend('HP trend','original','Location','Northwest'), end
    end
end
