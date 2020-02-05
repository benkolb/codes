
function MTH = gen_LP_correction_matrix_HerbstJohannsen(H,T)
% generate matrix for correction of small-sample bias of local projection
% regressions à la Herbst and Johannsen (2020) "Bias in Local Projections"
% 
% MTH = gen_LP_correction_matrix_HerbstJohannsen(H,T)
% 
% generate the matrix M_{T,H} as given on page 10 of the Fed WP:
%
% MTH = |          0                (1-1/T)/(T-1)        (1-2/T)/(T-1)     ...     (1-H/T)/(T-1)     |
%       |  (1-1/(T-1))/(T-2)              0            (1-1/(T-1))/(T-2)   ... (1-(H+1)/(T-1))/(T-2) |    ...
%       |         ...                    ...                  ...          ...          ...          |
%       | (1-H/(T-H))/(T-H-1)  (1-(H-1)/(T-H))/(T-H-1)        ...          ...           0           |
%
% by Benedikt Kolb, Deutsche Bundesbank, Feb. 2020


% 0. start with matrix (0 1 2 ... H; 1 0 1 ... H-1; 2 1 0 1 ... etc.) 
A0 = NaN(H+1,H+1);
HH = [H:-1:1,0:H];
for i = 1:H+1
        A0(i,:) = HH(end-H-i+1:end-i+1);
end

% 1. divide by T (first row), T-1 (second row) etc.
T1 = diag(T-(0:H));
A1 = T1\A0;

% 2. subtract this from 1's (except diagonal)
inveye = -(eye(H+1)-1);
A2 = inveye - A1;

% 3. divide by T-1 (first row), T-2 (second row) etc.
T3 = diag((T-(1:H+1)));
A3 = T3\A2;

MTH = A3;

end
