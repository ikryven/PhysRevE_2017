% Exapmle_1

% In this example, component size distribution is computed with three different methods.
% The example illustrates that the results of simple iterations, 
% Lagrangian inversion/fft algorithm  and point equations all coincide 
% for a sample degree distribution. Source:  I.Kryven, PhysRevE 2017.

n_max = 500;
nn = 0:n_max;


% degree distribution
u = exp( - nn / 1.05 );
u = u / sum( u );


mu1 = get_moments( u );


% methiod 1: simple iterations, Equation (2):
w_num  = components_num( u, n_max );

% methiod 2: Lagrangian inversion + fft convolution,  Equation (8)
w_lag  = components_lagrange( u, n_max );

% method 3: dirrect calculation,  Equation (9)
w_an(1) = u(1);
w_an(2) = 1 / mu1   *  u(2)^2;
w_an(3) = 3 / mu1^2 *  u(2)^2 * u(3);
w_an(4) = 4 / mu1^3 *  u(2)^2 * ( 2 * u(3)^2 + u(2) *u(4));
w_an(5) = 5 / mu1^4 *  u(2)^2 * ( 4 * u(3)^3 + 6* u(2) * u(3) * u(4) + u(2)^2 * u(5) );

%% plotting

w_num = w_num( 1:15 );
w_lag = w_lag( 1:15 );
nn    = 1:15;

cla

stairs( nn, w_num, '-b', 'LineWidth', 1.5   );
hold on
stairs( nn, w_lag, '--r', 'LineWidth',  1.5 );
plot( 1:5,   w_an,  'go', 'MarkerSize', 20  );
plot( 1:5,   w_an,  'g+', 'MarkerSize', 20  );

h = legend( 'Simple iterations', 'Lagrangian inversion + fft convolution', ' dirrect calculation' );
set( h, 'FontSize', 16 );

xlabel( 'Component size, n' );
ylabel( 'w(n)' );

set( gca', 'xscale','lin' );
set( gca', 'yscale','log' );


