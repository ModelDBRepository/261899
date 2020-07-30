%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% bisection for self-consistent steady-state
function [P0,p0,J0,r0,x0] = theory0(mu_in, sigma2, params, xi)   

gx = params(13);

%%% fixed-point iteration
num_x_its = 20;
x0_in = 0;
x0_out = 0;

if gx > 0
    for i=1:num_x_its

        [~,~,~,~,x0_out] = thin_x(params,x0_in,mu_in,sigma2,xi);
        x0_in = x0_out;

    end
end    

%%% self-consistent solution
[P0,p0,J0,r0,x0] = thin_x(params,x0_out,mu_in,sigma2,xi);

end