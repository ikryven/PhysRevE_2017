function [mu1, mu2, mu3] = get_moments( u )

nn  = 0:length(u)-1;
mu1 = sum( nn    .* u );
mu2 = sum( nn.^2 .* u );
mu3 = sum( nn.^3 .* u );