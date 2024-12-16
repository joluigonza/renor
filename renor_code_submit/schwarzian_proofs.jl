#############################################################
###################Schwarzian derivative 

df = derCheb_even2odd(f0); # f'
d2f = derCheb_odd2even(df); # f''
d3f = derCheb_even2odd(d2f);

f_fun(x)=eval_cheb_even(f0, x);
df_fun(x)= eval_cheb_odd(df, x);
d2f_fun(x)= eval_cheb_even(d2f, x);
d3f_fun(x)= eval_cheb_odd(d3f, x);

###########################################################
nb=1000
### numerics

thisdelta=1.0*10^(-10);
#thisdelta=1.0*10^(-8);

SfK1 = 2*d3f_fun(interval.(range(thisdelta,1,nb))).*df_fun(interval.(range(thisdelta,1,nb))) - 3*(d2f_fun(interval.(range(thisdelta,1,nb)))).^2

SfK1= sup.(real.(SfK1))

SfK_max= maximum(SfK1)
###################

M1K= sup(norm_C0_Enu1(df,rho))
M2K = sup(norm_C0_Enu1(d2f,rho))
M3K = sup(norm_C0_Enu1(d3f,rho))

r1min=sigma_new(interval(1.0),rho)*rmin  #constant in {lem:C0_VS_ell1} with h even with p=1 v=2

thisepsilon= (mid(rho)-1)/3

r2min=(2+3*thisepsilon)/3*sigma_new(interval(1.0),interval(1.0 + thisepsilon))*sigma_new(interval(1.0 + 2*thisepsilon),rho)*rmin
r2min=sup(r2min)

thisepsilon=(mid(rho)-1)/5

r3min=((1+2*thisepsilon)^2+(1+thisepsilon)^2)/((2 + 3*thisepsilon)*thisepsilon)*(2+7*thisepsilon)/thisepsilon* sigma_new(interval(1.0),interval(1.0 + thisepsilon))*sigma_new(interval(1.0 + 2*thisepsilon),interval(1.0 + 3*thisepsilon))*sigma_new(interval(1.0 + 4*thisepsilon),rho)*rmin
r3min=sup(r3min)

Sf_bound = 2*(M3K*r1min + M1K*r3min + r3min*r1min) + 3*r2min*(2*M2K+ r2min)

Sf_max= SfK_max + Sf_bound

sup(Sf_max)

#{lem:C0_VS_ell1}
#{prop:derivative}

##########################################################
foo= 1.0001

SfK2 = 2*d3f_fun(interval.(range(0.0,foo*thisdelta,nb))).*df_fun(interval.(range(0.0,foo*thisdelta,nb))) - 3*(d2f_fun(interval.(range(0.0,foo*thisdelta,nb)))).^2
+ 2*d3f_fun(interval.(range(0.0,foo*thisdelta,nb))).*interval(-sup(r1min),sup(r1min))
+ df_fun(interval.(range(0.0,foo*thisdelta,nb))).*interval(-sup(r2min),sup(r2min))
.+ interval(-sup(r1min),sup(r1min)).*interval(-sup(r2min),sup(r2min))
- 6*d2f_fun(interval.(range(0.0,foo*thisdelta,nb))).*interval(-sup(r2min),sup(r2min))
.- 3*interval(-sup(r2min),sup(r2min))^2;

maximum(sup.(real.(SfK2)))

minimum(inf.(real.(SfK2)))

##########################################
using PyPlot
pygui(true)

plot(range(thisdelta,1,nb),SfK1)

plot(range(0.0,foo*thisdelta,nb), sup.(real.(SfK2)))

