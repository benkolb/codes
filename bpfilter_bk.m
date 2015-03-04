
function [tr, cyc] = bpfilter_bk(y,pl,pu,graph)

% Applies the bandpass filter as in Christiano and Fitzgerald (1999).
%
% Usage:   [tr, cyc] = bpfilter_bk(y,pl,pu,graph)
%
% Inputs:  y:     original series
%          pl:    minimum duration of desired cycle (2 <= pl < pu)
%          pu:    maximum duration of desired cycle (pu > pl)
%          graph: (optional) If graph=1, the results will be plotted. 
% Outputs: tr:    Trend component of filtered series
%          cyc:   Cyclical component of filtered series
%
% Examples: Quarterly data: pl=6, pu=32 returns component with periods 
%                           between 1.5 and 8 years
%           Monthly data:   pl=2, pu=24 returns component with all periods
%                           less than 2 years 
%
% by Benedikt Kolb (21/01/14), largely based on code by Eduard Pelz
%
% For a derivation of the formulas used, see e.g.
% http://economics.sas.upenn.edu/~jesusfv/lecturetechnical5_spectral.pdf


if pu <= pl
   error(' (BP filter): pu must be larger than pl')
end
if pl < 2
   warning('bp:pl_low',' (BP filter): pl less than 2, I reset it to 2')
   pl = 2;
end

if size(y,1) < size(y,2)
    y = y';
end

[T N] = size(y);

%% Remove drift (CF filter requires input w/o drift) 
% drift = (y(T) - y(1)) / (T-1).                     
tim   = 0:T-1;
drift = (y(T,1)-y(1,1)) / (T-1);
yun   = y - tim'*drift*ones(1,N);

%% Create the ideal Bs, then construct the AA matrix
b    = 2*pi/pl;
a    = 2*pi/pu;
bnot = (b-a)/pi;
bhat = bnot/2;

B = ( (sin(tim*b)-sin(tim*a))./(tim*pi) )';
B(1,1) = bnot;

AA = zeros(2*T,2*T);

for i=1:T
   AA(i,i:i+T-1) = B';
   AA(i:i+T-1,i) = B;   
end

AA = AA(1:T,1:T);

AA(1,1) = bhat;
AA(T,T) = bhat;

for i=1:T-1
   AA(i+1,1) = AA(i,1)-B(i,1);
   AA(T-i,T) = AA(i,1)-B(i,1);
end

%% Filtering
cyc = AA*yun;
tr  = y - cyc;


%% (optional) Plot BP-filtered vs original series
if nargin == 4 && graph == 1
    figure('Name','BP-filtered vs original series')
    for i = 1:N
        subplot(max(floor(N/3),1),ceil(N/max(floor(N/3),1)),i);
        plot(1:T,tr(:,i),'r',1:T,y(:,i),'k--'); grid on; ...
            title(['Series #',num2str(i)]);
        if i==1, legend('BP trend','original','Location','Northwest'), end
    end
else
    return
end
    