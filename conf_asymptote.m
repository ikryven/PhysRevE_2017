function w = conf_asymptote( nn, beta, s, mu1, mu2, mu3 )
%CONF_ASYMPTOTE  analytical asymptote for component size distribution in
%   configuration network.
%   w = conf_asymptote ( NN, BETA, S, MU1, MU2, MU3 ) computes the asymptote
%   for component size distribution in a non-directional configuration 
%   network at points NN.
%   NN   is a vector of points with elements greater than 1;
%   The degree distribution is parametrised with:
%   BETA, the exponent of the heavy-tailed degree distribution,
%         set BETA=inf for light tails;             
%   S,    the scale of the heavy tailed degree distribution;
%   MU1,MU2,MU3, first three moments of the degree distribution.
%   
%   The notation and comments are in accordance with
%   Table II, I.Kryven, PhysRevE 2017: "General expression for
%   the component size distribution in infinite configuration networks". 
%
%   Licensed under, CC BY, April, 2017. For attribution refer to the publication.

%% constants

    theta = mu2 - 2 * mu1;
    alpha = beta - 2;
    
    C1  = mu1^2 / sqrt(   2 * pi * ( mu1 * mu3 - mu2^2  )  );
    C2  = ( mu2 - 2 * mu1 )^2 / ( 2 * mu1 * mu3 - 2 * mu2^2  );
    C1a = mu1 / sqrt(   2 * pi * s );    
    C2a = ( 2 * mu1 - mu2 )^2 / ( 2 * s * mu1^2 );
    C3  = s * mu1^( alpha + 2 ) / ( 2 * mu1 - mu2 )^( alpha + 1 ) * gamma( alpha + 1 ) / gamma( alpha ); 
    C4  = mu1 * pi^( - 1 / alpha - 1 ) * gamma( 1 + 1 / alpha ) * sin( pi / alpha ) * ( 2 / s * ( gamma( alpha ) * sin( pi * alpha / 2 ) ) )^( 1 / alpha );
    C5  = mu1 / sqrt( alpha - 1 ) * (  2^( 2 - alpha ) * ( mu2 / mu1 - 2 )^( 2 - alpha ) * gamma( alpha ) * sin( pi * alpha / 2 ) / ( alpha * pi^alpha * s )     )^( 1 / ( 2 * alpha - 2 ));
    C6  = ( 1 - alpha ) * ( 2 * ( mu2 / mu1 - 2 )^alpha * gamma( alpha ) * sin( ( pi * alpha ) / 2 ) /  ( alpha^( alpha ) * pi * s ) )^( ( 1 / ( alpha - 1 ) ) )        ;
    C7  = mu1 * sqrt( 2 ) / ( pi^( 3 / 2 ) * s );
    C8  = 1 / ( pi * s ) + 1 / 2;
    C9  = exp( - 1 - 2 / ( pi * s ) );
    C10 = mu1 / sqrt( 2 - 2 * alpha )  * ( ( gamma( alpha ) * sqrt( 2 ) * sin( pi * alpha / 2 ) ) / ( s * pi^alpha * alpha ) )^( 1 / ( 2 * alpha - 2 )) ;
    C11 = ( 1 - alpha ) * ( ( ( sqrt( 2 ) * sin( (pi * alpha ) / 2 ) * gamma( alpha ) ) / ( pi * s * alpha^alpha ))^( (1 / ( alpha - 1 )) )) ;

%% cases of beta and theta   

    if isinf( beta )                                               % CASE A

            if     abs( theta  ) < 0.001

                   w  = C1 * nn.^( - 3 / 2 );
                   disp( 'critical regime, light tail' );

            else 

                   w  = C1 * exp( - C2 * nn ).* nn.^( - 3 / 2 );
                   disp( 'sub/sup-critical regime, light tail' );

            end;

    elseif beta > 4.0000001                                        % CASE B

           if     abs( theta  ) < 0.001

                  w  = C1 * nn.^( - 3 / 2 );
                  disp( 'critical regime,  alpha > 2 ' );

           elseif theta < 0 

                  w = C3 * nn.^( - alpha - 1 );
                  disp( 'sub-critical regime, alpha > 2' );

           elseif theta > 0 

               w  = C1 * exp( - C2 * nn ).* nn.^( - 3 / 2 );
               disp( 'supper-critical regime, alpha > 2' );

           end;



    elseif beta == 4                                               % CASE C

           if     abs( theta ) < 0.001

                  w  = C1a * nn.^( - 3 / 2 ) ./ sqrt( log( nn ) );
                  disp( 'critical regime,  alpha = 2 ' );

           elseif theta > 0

                   w  = C1a * nn.^( - 3 / 2 ) ./ sqrt( log( nn ) ) .* exp( - C2a * nn ./ log( nn ) );
                   disp( 'supper-critical regime, alpha = 2' );

           elseif theta <  0

                  w = C3 * nn.^( - alpha - 1 );
                  disp( 'sub-critical regime, alpha = 2' );

           end;


    elseif beta > 3 && alpha < 4                                   % CASE D

           if     abs( theta ) < 0.01

                  w = C4 * nn.^( - 1 / alpha - 1 );
                  disp( 'critical regime,  alpha in ( 1, 2 ) ' );

           elseif theta < 0

                  w = C3 * nn.^( - alpha - 1 );
                  disp( 'sub-critical regime, alpha in ( 1, 2 ) ' );

           elseif theta > 0

                  w = C5 * nn.^( - 3 / 2 ) .* exp( C6 * nn );
                  disp( 'supper-critical regime, alpha in ( 1, 2 )' );

           end;

    elseif beta == 3 

           w = C7 * exp( - C8 - C9 * nn .^( 2 / pi ) ) .* nn .^ ( 1 / pi - 2 );

    elseif beta > 2 && beta < 3                                    % CASE E

           w = C10 * exp( - C11 * nn ) .* nn .^ ( - 3 / 2 );

    elseif beta >1 && beta < 2                         % SEMI-DENSE NETWORK                         

           disp( 'Mean value of the excess degree distribtion is divergent' );
           w = nn * NaN;

    else                                                    % DENSE NETWORK           

           disp( 'No finite components: Mean value of the degree distribtion is divergent' );
           w = nn * 0;

    end;

 