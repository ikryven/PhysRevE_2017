% Exapmle_2

% In this example, an oscillatory component size distribution is computed with 
% Lagrangian inversion/fft algorithm. The component size distribution
% slowly converges to asymptote. The example corresponds to Figure 1 in
% I.Kryven, PhysRevE 2017.

n_max = 1e3;
nn    = 1:n_max;


% degree distribution

u = zeros( 1, n_max + 1 );

u( 2 )  = 0.97;
u( 3 )  = 0.015;
u( 11 ) = 0.015;

u = u / sum( u );


[ mu1, mu2, mu3 ] = get_moments( u );

beta = inf;
s    = 1;


% Lagrangian inversion + fft convolution
w_lag = components_lagrange( u, n_max );

% asymptote
w_a   = conf_asymptote( nn, beta, s, mu1, mu2, mu3 );

%% plotting


cla
stairs( nn, w_lag, '-b', 'LineWidth', 1.5 );
hold on
plot( nn, w_a , '-r', 'LineWidth', 1.5 );


h = legend( 'Lagrangian inversion + fft convolution', ' Asymptote' );
set( h, 'FontSize', 16 );

xlabel( 'Component size, n' );
ylabel( 'w(n)' );

set( gca', 'xscale', 'log' );
set( gca', 'yscale', 'log' );