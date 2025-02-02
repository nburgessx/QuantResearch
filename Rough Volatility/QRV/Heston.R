    
#-----------------------------------------------------------------------------
# Option pricing from characteristic function
# Equation (5.6) of The Volatility Surface

# callOption <- function(phi, k, t){
#     integrand <-  function(u){Re(exp(-1i*u*k)*phi(u - 1i/2, t)/(u^2 + 1/4))};
#     res <- 1 - exp(k/2)/pi*integrate(integrand,lower=0,upper=Inf,rel.tol=0.0000000001,subdivisions=1000)$value;
#     return(res);
# }
# 
# bsvol <- function(phi,k, t){BSImpliedVolCall(1, exp(k), t, 0, callOption(phi,k,t))}

#-----------------------------------------------------------------------------
# Heston characteristic function

#subHeston <- list(lambda = 0.6067,rho = -0.7571,eta = 0.2928,vbar = 0.0707,v = .0654);

# The Heston model
phiHeston <- function(params){
    
    lambda <- params$lambda;
    rho <- params$rho;
    eta <- params$eta;
    vbar <- params$vbar;
    v <- params$v;
    
    function(u, t){

    al <- -u*u/2 - 1i*u/2;
    bet <- lambda - rho*eta*1i*u;
    gam <- eta^2/2;
    d <- sqrt(bet*bet - 4*al*gam);
    rp <- (bet + d)/(2*gam);
    rm <- (bet - d)/(2*gam);
    g <- rm / rp;
    D <- rm * (1 - exp(-d*t))/ (1 - g*exp(-d*t));
    C <- lambda * (rm * t - 2/eta^2 * log( (1 - g*exp(-(d*t)))/(1 - g) ) );
    return(exp(C*vbar + D*v));
}
}


# #-----------------------------------------------------------------------------
# # Function to generate implied vols. corresponding to Heston model
# 
# impvolHeston <- function(params){
#     function(k, t){
#     	bsvol(phiHeston(params),k,t)
#         }
#     }
#     
# #impvolHeston(subHeston)(0,1)
