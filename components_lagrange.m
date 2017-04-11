function [ w ] = components_lagrange( u, n_max )
%COMPONENTS_LAGRANGE  exact component size distribution in
%   configuration network.
%   W  = COMPONENTS_LAGRANGE( U, N_MAX ) computes the component size
%   distribution in configuration network with an exact algorithm.
%   The convolution operation is implemented with FFT.
%
%   U  provides the degree distribution, so that u(1) is the probability 
%   of sampling a node with zero connections, u(2) probability of sampling 
%   a node with 1 connection, etc.
%   N_MAX provides the maximum component size.
%   W is a vector, such that W(n), n=1..n_max gives probability that 
%   a randomly sampled node belongs to a finite component of size n.
%   
%   The notation and comments are in accordance with
%   Equation (8), I.Kryven, PhysRevE 2017: "General expression for 
%   the component size distribution in infinite configuration networks". 
%   
%   Licensed under, CC BY, April, 2017. For attribution refer to the publication.

    %%  initialize

    w  = zeros( 1, n_max );
    u( n_max + 2 : end ) = [ ];
    u( n_max + 1 ) = 0;

    u = u / sum( u );
    nn = 0:length( u ) - 1;

    %% compute excess distribution
    

    u0 = u .* nn;
    u0( 1 ) = [];

    mu1 = sum( u0 );

    u1 = u0 / mu1;
    u1( end + 1 ) = 0;


    %% convolution powers
    
    fft_u1  = fft( u1 );
    fft_u1n = fft_u1;
    
    for n = 2 : n_max - 1; 
        
        fft_u1n = fft_u1n .* fft_u1;
        u1n     = real( ifft( fft_u1n ) );
        
        w( n ) = mu1 / ( n - 1 ) * u1n( n - 1  ) ;
           
    end;
      
    w( 1 ) = u( 1 );  

    
    