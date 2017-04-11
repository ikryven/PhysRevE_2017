% Exapmle_4

% In this example, a component size distribution is computed for 
% a degree distribution that features a heavy-tail and finite moments. 
% The example demonstrates the emergence of a transient asymptote with
% exponent -2/3, while a switch to a true asymptote (exponent -5 ) takes 
% place only for n > 10^5.
% This example requires around 15 min of cpu time on a single core.
%
% See also Figure 4 in I.Kryven, PhysRevE 2017.

n_max = 1e4;
nn = 0:n_max;


% degree distribution

beta = 6;  % degree exponent 

s = 9.4244;
u   =  s*( beta - 2 ) * nn .^ ( -beta );

u( 1:2 ) = 0;
u(  2 ) = 1 - sum( u );
u = u / sum( u );

[mu1 mu2 mu3] = get_moments( u );

nn = 1:n_max;


cla
 
% asymptote
w_a  = conf_asymptote( nn, beta, s, mu1, mu2, mu3 );
plot( nn, w_a , '--r', 'LineWidth', 1.5 );

hold on
% transient asymptote
C1 = mu1^2 / sqrt(   2 * pi * ( mu1 * mu3 - mu2^2  )  );
w_t  =  C1 * nn.^( - 3 / 2 );

plot( nn, w_t, '--g', 'LineWidth', 1.5 );


refresh;
drawnow;

% Lagrangian inversion + fft convolution,  Equation (8)
disp('Computing Lagrangian inversion up to 10^5, will take around 15 min...')
tic
w_lag = components_lagrange( u, n_max );
toc



%% plotting

plot( nn, w_lag, '-b', 'LineWidth', 1.5 );
plot( nn, w_t, '--g',  'LineWidth', 1.5 );
plot( nn, w_a , '--r', 'LineWidth', 1.5 );

h = legend(' Asymptote',' Transient asymptote', 'Lagrangian inversion + fft convolution' );
set( h, 'FontSize', 16 );

xlabel( 'Component size, n' );
ylabel( 'w(n)' );

set( gca', 'xscale', 'log' );
set( gca', 'yscale', 'log' );