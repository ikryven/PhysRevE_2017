function w  = component_num( u, n_max )
%COMPONENT_NUM   component size distribution in configuraiton model
%   W = COMPONENT_NUM( u, n_max) computes the component size
%   distribution in configuration network via simple point iteration.
%
%   U  provides the degree distribution, so that u(1) is the probability 
%   of sampling a node with zero connections, u(2) probability of sampling 
%   a node with 1 connection, etc.
%   N_MAX provides the maximum component size.
%   W is a vector, such that W(n), n = 1,..,n_max gives probability that 
%   a randomply sampled node belongs to a finite component of size n.
%   
%   The notation and comments are in accordance with
%   I.Kryven, PhysRevE 2017: "General expression for 
%   the component size distribution in infinite configuration networks". 
%   
%   Licensed under, CC BY, April, 2017. For attribution refer to the publication.
                                                   


warning off

F_tol = 1e-11;

nn=1:n_max;


%% Generating funcitons

U = u( end:-1:1 );
U = U / sum( U );

U1 = polyder( U );
U1 = U1 ./ sum( U1 );

x = linspace( 0, 2*pi, n_max*4 );
xi = cos(x) + 1i * sin( x ) ;


W1 = ones( 1, length( xi ) );
W1 = W1 / W1( end );

%% Fixed point iteration
nrm = 1;
while nrm > F_tol;
    
    W10 = W1;
    W1 = xi .* polyval( U1, W1 );
    
    nrm = norm( W1 - W10 );
    
end;

W = xi .* polyval( U, W1 );


%% Compute Cauchy integral
for k=1:length(nn)

    w( k )  = real( trapz( xi, W  .* xi.^( -nn( k ) - 1 ) ) / ( 2 * pi * 1i ) );
    
end;

