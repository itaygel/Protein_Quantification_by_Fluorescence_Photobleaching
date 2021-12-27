function [fitresult, gof] = two_exp_fit(x,a)
%CREATEFIT(A)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      Y Output: a
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 31-Dec-2019 13:15:37


%% Fit: 'untitled fit 1'.
%[xData, yData] = prepareCurveData( [], a );
[xData, yData] = prepareCurveData( x, a );
% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c*exp(-d*x)+f', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [1000 0 1000 0 -200];
opts.StartPoint = [a(1)/2 -log(a(5)/a(4)) a(1)/2 -log(a(5)/a(4)) 0];
opts.Upper = [Inf 2 Inf 2 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

%Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'a', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% ylabel a
%grid on


