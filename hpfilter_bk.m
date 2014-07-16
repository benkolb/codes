

%   [s] = hpfilter(y,w,'makeplot')
%   'makeplot' in the input, plots the graphics of the original series
%   against the filtered series, if more than one series is being
%   considered the program will plot all of them in different axes

function fltd = hpfilter_bk(y,w,graph)

% Applies the Hodrick-Prescott filter.
% Usage:   fltd = hpfilter_bk(y,w,graph)
% Inputs:  y:     original series
%          w:     smoothing parameter (choose 100 for annual, 1600 for 
%                 quarterly and 14400 for monthly data)
%          graph: logical; if 1, filtered series will be plotted    
% Outputs: fltd: filtered series
% by B.Kolb (21/01/14), based on code by Wilmer Henao. For a derivation of the formulas used, see e.g. http://economics.sas.upenn.edu/~jesusfv/lecturetechnical5_spectral.pdf


if nargin == 1
    warning('hp:no_w','No smoothing parameter w specified. I assume you want to use w = 1600 in the following.');
end

if size(y,1) < size(y,2)
    y = y';     
end
[t n] = size(y);

d = repmat([w -4*w ((6*w+1)/2)], t, 1);
d(1,2) = -2*w;      d(t-1,2) = -2*w;
d(1,3) = (1+w)/2;   d(t,3) = (1+w)/2;
d(2,3) = (5*w+1)/2; d(t-1,3) = (5*w+1)/2;
B = spdiags(d, -2:0, t, t);    %I use a sparse version of B, because when m is large, B will have many zeros     
B = B+B';
fltd = B\y;

if nargin(3) == 1
    for i = 1:t
        figure(i)
        subplot(floor(n/3),ceil(n/floor(n/3)),i);
        plot(fltd(:,i),'r');   grid on;   hold on;   plot(y(:,i));   title(['Series #',num2str(i)]);
    end
end
