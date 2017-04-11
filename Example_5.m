% Exapmle_3

% In this example, a component size distribution is computed for 
% a degree distribution that features a heavy-tail, -2.6.
% In this case, the excess degree distribution has a divergent mean value.
% The example demonstrates presence of transient asymptotic modes in with
% algebraic decay and exponent -1.6;
%
% See also Figure 3 in I.Kryven, PhysRevE 2017.


beta = 2.6;     % degree distribution exponent 
sa   = 8.33e-2; % degree distribution scale for case 1
sb   = 8.33e-5; % degree distribution scale for case 2



n_max = 5e4;
nn = 0:n_max;

%% sub-critical case

ua   =  sa * ( beta - 2 ) * nn .^ ( - beta );

ua( 1:2 ) = 0;
ua(  2 )  = 1 - sum( ua );
ua = ua / sum( ua );

[mu1a, mu2a, mu3a ] = get_moments( ua );


%% critical case

ub   =  sb * ( beta - 2 ) * nn .^ ( - beta );

ub( 1:2 ) = 0;
ub(  2 )  = 1 - sum( ub );
ub = ub / sum( ub );

[ mu1b, mu2b, mu3b ] = get_moments( ub );


%%
nn = 1:n_max;


w_lag_a = components_lagrange( ua, n_max );
w_lag_b = components_lagrange( ub, n_max );

w_a_b   = conf_asymptote( nn, beta, sb, mu1b, mu2b, mu3b );
w_a_a   = conf_asymptote( nn, beta, sa, mu1a, mu2a, mu3a );

%%

% transient asymptotes
 
  w_tb = mu1b * sb * gamma( beta - 1 ) / gamma( beta - 2 ) * nn.^( - beta + 1 );
  w_ta = mu1a * sa * gamma( beta - 1 ) / gamma( beta - 2 ) * nn.^( - beta + 1 );

%% plotting
cla
hold on
plot( nn, w_lag_a, '-b', 'LineWidth', 1.5 );
plot( nn, w_lag_b, '-b', 'LineWidth', 1.5 );

plot(nn,w_ta,'-')
plot(nn,w_tb,'-')


plot( nn, w_a_a , '--r', 'LineWidth', 1.5 );
plot( nn, w_a_b , '--g', 'LineWidth', 1.5 );


h = legend('Lagrangian inversion' , 'Lagrangian inversion' ,'Transient asymptote, \eta=-1.6','transient asymptote, \eta=-1.6 ','True asymptote', 'True asymptote ' );
set( h, 'FontSize', 16 );

xlabel( 'Component size, n' );
ylabel( 'w(n)' );

set( gca', 'xscale', 'log' );
set( gca', 'yscale', 'log' );