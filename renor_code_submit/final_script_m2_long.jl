
import Pkg;

export dropdims


using LinearAlgebra

using IntervalArithmetic

using GenericLinearAlgebra

using GenericSchur

using FileIO

using JLD2

using Serialization

using TickTock

#using PyPlot

#############################

include("general_renor_functs_final.jl")

################################
pre=2^12
setprecision(pre);
setprecision(Interval, pre);


m=2; # Renormalization order #######################

K=620; # Chebyshev x order


rho = Interval(BigFloat(2.0)); #Bernstein ellipse radius

rstar = BigFloat(10.0)^(-585) # Neighborhood radius in the contraction proofs


##########################################



theta_K, MK, MKinv= cheb_nodes_even(K);
x_K = cos.(theta_K);


#########################################

global jacF


thisdoc=open("thisJacobian"*string(m)) 


thisstring=read(thisdoc,String);
thisJac=replace(thisstring,"\n"=>"");
thisJac=replace(thisJac,";"=>"");

close(thisdoc)



jacF=max_jac_formula0(thisJac,m);

thismiu=get_miu(m)*1000/BigFloat(1000);

thisvals= 1.0 .- thismiu.*x_K.^2;
thiscoeffs=MKinv*thisvals .+ BigFloat.(zeros(K+1,1)).*im;

#####################################################################

eval0 = (-1).^(0:K);
eval0[2:end] = 2*eval0[2:end];
eval0=transpose(eval0);

f=thiscoeffs;

thisf0=eval0*f;


alpha =BigFloat(0.0);
for i=1:m
   global alpha =eval_cheb_even(f,alpha);    
end

alpha=copy(alpha[1]);


alpha0, Rf, DRf, DRalpha = Rm_DRm_even_1D0(f,m,x_K,MKinv,jacF,alpha);
res = [thisf0[1]-1.0; alpha*f-Rf]; #F(f)

Jac = [0 eval0; DRalpha alpha*eye(K+1) - DRf]; #DF(f)

err = sum(abs.(res));
thiserr=Float64(err,RoundNearest)

for it=0:10 # Newton's iterations
    println(it);

    thisd= - Jac\res;
    global    alpha=alpha + thisd[1];
    global f=f + thisd[2:end];


    global thisf0=eval0*f;

    global alpha0, Rf, DRf, DRalpha = Rm_DRm_even_1D0(f,m,x_K,MKinv,jacF,alpha);
    global res = [thisf0[1]-1.0; alpha*f-Rf]; #F(f)

    global Jac = [0 eval0; DRalpha alpha*eye(K+1) - DRf]; #DF(f)

    global err = sum(abs.(res));
    global thiserr=Float64(err,RoundNearest);
    println(thiserr)

end

##############################

for k=1:m
  
    alpha0, Rf=Rm_even_1D(f,k,x_K,MKinv);    
    res = f- Rf ./alpha0;
    err = sum(abs.(res));
    thiserr=Float64(err,RoundNearest);
    println(thiserr)

end



#################################
f0=copy(f);
serialize("f0"*"_"*string(m)*"long", f0)

alpha0=copy(alpha);
serialize("alpha0"*"_"*string(m)*"long", alpha0);


###################################
#f= deserialize("f0")

##############################################################
#############################################################################


tick()

##########################################
J=inv(Jac);


f = Interval.(f);
alpha = Interval(alpha);
J = Interval.(J);


####################################################
##############################  YK

ipi= @interval pi 


theta_K = (K.-(0:K)).*(ipi/BigFloat(2*K)); #In the even case, we only need the nodes in [0,1]
ind_k = 0:2:2*K;
x_K = cos.(theta_K);

# Construction of the matrix to go from values at Chebyshev nodes to
# Chebyshev coefficients
MKinv = transpose(cos.(theta_K*transpose(ind_k)))/K; ######
MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]]/2;
MKinv[K+1,:] = MKinv[K+1,:]/2;

# weights, with a 1 for the alpha component
weightsrho=[1; 1; 2*rho.^ind_k[2:end]];

df = derCheb_even2odd(f); # f'
d2f = derCheb_odd2even(df); # f''

###############################################
f_coeffs = copy(f);
df_coeffs = copy(df);

f_fun(x)=eval_cheb_even(f, x);

f_num=mid.(f);
f_fun_num(x)=eval_cheb_even(f_num, x);

df_fun(x)= eval_cheb_odd(df, x);

df_fun_num(x)=eval_cheb_odd(mid.(df), x);

d2f_fun(x)= eval_cheb_even(d2f, x);

d2f_fun_num(x)= eval_cheb_even(mid.(d2f), x);

### YK
#fprintf('Computing YK...\n')
# 1st estimate


thisf0=eval0*f;

this, Rf, DRf, DRalpha = Rm_DRm_even_1D0i(f,m,x_K,MKinv,jacF,alpha);

Phi = [thisf0[1]-1.0; alpha*f-Rf]; #F(f)
Jac = [0 eval0; DRalpha alpha*eye(K+1) - DRf]; #DF(f)

# 2nd estimate: We compute an explicit bound of |phi(alpha,f)| on an 
# ellipse, in order to get a geometrically decaying bound on each Chebyshev
# coefficient of ipi^K(phi(alpha,f))

alpha_num=mid(alpha);


######################################################
# to guarantee iterations are in a small enough ellipse for eigenvalue problem (after ZKK for fp we compute the check for the fp problem). saving until eigenvalue problem
thist = Interval(maximum(sup.(abs.(fksalphax))));

s_eig = sup(thist + 2*sqrt( thist^2 -1)); 


#################################################
#################################################

cste_YK(rho,alpha) = 1;

#################### 

YK = transpose(weightsrho)*abs.(J*Phi);

YK=sup(YK[1])

###############################################
YK0=copy(YK);

serialize("YK0", YK0)

################################################################
####################### Yinf

#########################################################
#fprintf('\nComputing Yinf...\n')


###########################################################

function cste_Yinf(rho,thisalpha)

    if isa(thisalpha,Interval{BigFloat})
        return 1/abs(alpha)*Upsilon(rho,thisalpha,1,0,2*K);
    else
        return 1/abs(mid(alpha))*Upsilon(mid(rho),thisalpha,1,0,2*K);        
    end

end

############### alpha_Yinf precomputed for smaller K
tol=10^(-2);

alpha_max = 9;


alpha_Yinf = Interval(BigFloat(7.341228274042058)) ;


nb=10^4;
sub = BigFloat.(range(0,2*pi,nb)); #sub==int in original code



ffl_coeffs=MKinv*ffl(x_K);

@time Yinf = cste_Yinf(rho,alpha_Yinf)*norm_C0_Enu1(ffl_coeffs,alpha_Yinf);


Yinf=sup(Yinf)

###################################

Yinf0= copy(Yinf)
serialize("Yinf0", Yinf0)


##
###################################################################
########################## Z1KK
### Creating \tilde{f} for Z1

tr = 0 ± rstar;
dtr = sigma(1,rho,0,1)*tr; #the error bound for f' ###################### 
tl = alpha + tr;

tf(x)= f_fun(x) .+ tr;

tf_num(x)= f_fun_num(x);

dtf(x)= df_fun(x) .+ dtr;

dtf_num(x)= df_fun_num(x);


thisf0=eval0*f;
alpha0, Rf1, DRf1, DRalpha1 = tRm_DRm_even_1D0i(f,m,x_K,MKinv,jacF,alpha);
Phi1 = [thisf0[1]-1.0; alpha*f-Rf1]; #F(f)
Jac1 = [0 eval0; DRalpha1 alpha*eye(K+1) - DRf1]; #DF(f)


Z1KK_mat=Interval.(BigFloat.(eye(K+2))) - J*Jac1;

this=(transpose(weightsrho)*abs.(Z1KK_mat))./transpose(weightsrho);

Z1KK=norm(sup.(norm.(this,Inf)),Inf)

###############################

Z1KK0= copy(Z1KK)

serialize("Z1KK0", Z1KK0)

##########################################
# to guarantee iterations are in a small enough ellipse for fp problem 
thist = Interval(maximum(sup.(abs.(fksalphax))));

thiss = sup(thist + 2*sqrt( thist^2 -1)); 


##################################
###########################################################################
################################ Z1Kinf

#fprintf('\nComputing Z1Kinf...\n')

################################
# 1st estimate
# Constructing bound in formula 14. Notice that df_fksalphax[:,j] are evaluated using tf and the nodes when ZKK is computed. The derivative formulas explicit in section 5.1 and stored in the text file were in practice computed symbolically. The same comment applies to the Z1Kinf bound for the eigenvalue problem. 

global Z1KinfUpsilon

Z1KinfUpsilon=Upsilon(thiss,rho,0,1,2*K);



thisZ1Kinf=jacF;
this=@Name(Z1KinfUpsilon);
for j=1:m 
    global   thisZ1Kinf=replace(thisZ1Kinf, "h_fksalphax[:,:,"*string(j)*"]" => this);
    global   thisZ1Kinf=replace(thisZ1Kinf, "df_fksalphax[:,"*string(j)*"]" => "abs.(df_fksalphax[:,"*string(j)*"])");
    
end

PiKDPhi = [0; abs.(MKinv)*eval(Meta.parse(thisZ1Kinf))];

##############################################################
####################################

func_alpham(z)=func_alpha(z,m,0).*dtf(tl.*z).*z;

# 2nd estimate

thistl= norm(inf.(real.(tl)),Inf);

vartheta1(beta)= ((thistl)*(beta- 1 ./ beta) .+ sqrt.(4 .+ ((thistl)*(beta-1 ./beta)).^2))./2;


cste1_Z1Kinf(rho,beta)= Upsilon(vartheta1(beta),rho,0,1,2*K);

this=real((rho-1/rho)/abs(tl)+sqrt(4+((rho-1/rho)/abs(tl))^2))/2;

beta_max = inf(this);

nb=100;
sub = BigFloat.(range(0,2*pi,nb)); 

term1(z)=func_alpha(z,m,0);  

term1_coeffs=MKinv*term1(x_K);
@time est_PiKDPhi_1 = est_PiKg_cste1([term1],cste1_Z1Kinf,K,rho,beta_max,sub,nb); 


gammas=interval.(BigFloat.(zeros(m,1)));
gammas[1]=beta_max;
thissum=est_PiKDPhi_1;

###########################################


for n=1:m-1
    
    thisfun(z)=tffl(z,n);    
    
    thisvartheta(gamma)=inclusion_gE1(thisfun,gamma,sub,nb);

    this_cste_Z1Kinf(rho,gamma)= Upsilon(thisvartheta(gamma),rho,0,1,2*K);

        
    thisgamma = opt_incl1(thisfun, inf(rho),alpha_max,tol,sub,nb);
    
    
    thisgamma = Interval(thisgamma);

      this1=inclusion_gE1(thisfun,thisgamma,sub,nb);
        
    while sup(this1[1]) > rho
        
        thisgamma = 0.99*thisgamma;    
       
        this1=inclusion_gE1(thisfun,thisgamma,sub,nb);    
    end

    gammas[n+1]=thisgamma;


    k=m-n;
    
    thisterm(z)= func_alpha(z,k,n);  
    
    thisterm_coeffs= MKinv* thisterm.(x_K);

    @time this_est_PiKDPhi = est_PiKg_cste1([thisterm], this_cste_Z1Kinf,K,rho,sup(thisgamma),sub,nb); 


   global thissum=thissum+this_est_PiKDPhi;

end


PiKDPhi[2:end] =  min.(real.(PiKDPhi[2:end]), thissum);

boundDR0f= PiKDPhi[2:end]; # for eigenvalue proof

Z1Kinf = transpose(weightsrho)*abs.(J)*PiKDPhi;

Z1Kinf=sup(real(Z1Kinf))

##########################

Z1Kinf0= copy(Z1Kinf)

serialize("Z1Kinf0", Z1Kinf0)

########################################
##################################################################
############################### Z1inf


nb=100;
sub = BigFloat.(range(0,2*pi,nb)); 

### Z1inf
#fprintf('\nComputing Z1inf...\n')


cste_Z1inf(rho,alpha)= Upsilon(rho,alpha,1,0,2*K);
 

alpha_Z1inf = BigFloat(7.341228274042058);


this=Interval(alpha_Z1inf[1]);

func_alpham_coeffs= MKinv*func_alpham(x_K);

Z1inf_alpha = cste_Z1inf(rho,this)*norm_C0_Enu1(func_alpham_coeffs,this);

Z1inf_alpha

###########################

betas=BigFloat.(zeros(m,1));
betas[1]=beta_max;
thissum=BigFloat(0.0);

for n=0:m-2
    
    thisfun(z)=tffl(z,n+1);
    
    thisbeta=opt_incl1(thisfun, inf(rho),alpha_max,tol,sub,nb);

    betas[n+2]=thisbeta;
    thisbeta=Interval(thisbeta);

    k=m-n;
    
    thisterm(z)= func_alpha(z,k,n);   

    # precomputed
    # this1= opt_para(thisterm, mid(rho), cste_Z1inf, betas[n+1],tol,sub);    
    # this=this1[1];
    foo=Interval(BigFloat(4.000002385341133258989759384247602424431142338487781069929132136601670291309435));
    
    thisterm_coeffs= MKinv*thisterm(x_K);
    thisbound = cste_Z1inf(rho,foo)*norm_C0_Enu1(thisterm_coeffs,foo);


    global thissum=thissum+thisbound

end

#add last one h(f(m)) at the end
thissum=thissum + Upsilon(rho,betas[end],1,0,2*K);

Z1inf_f= thissum; # for eigenvalue proof

Z1inf= 1/abs(tl) * ( max(Z1inf_alpha, Z1inf_f ) + tr);

Z1inf=sup(Z1inf)

##############################################

Z1inf0= copy(Z1inf)

serialize("Z1inf0", Z1inf0)

############################################

################################################
######################################## rstar

### Checking whether we can apply Banch fixed point theorem
Y = YK+Yinf;
Z1 = Z1KK + Z1Kinf + Z1inf;
rmin = sup(Y/(1-Z1));

if Z1>1
   error("negative Z1>1 not allowed")
end

println(rmin)
rmax = rstar;
println(rmax)

if rmin < rmax
    println("Proof successful");
    
else  
    println("No proof");

end

###############################

rmin0= copy(rmin)

serialize("rmin0_"*string(m), rmin0)

##

###################################################
############################ start eigen
#############################################

f0=copy(f);

alpha0=copy(alpha);

##

global jacF_alpha
jacF_alpha= get_JacF(m);


alphaf, Rf0, DRf0, = Rm_DRm_even_1Dv1i(f0,m,x_K,MKinv,jacF_alpha,alpha);

##
###########################################################################
################################ Z1Kinf

#fprintf('\nComputing Z1Kinf...\n')

################################
# 1st estimate

global Z1KinfUpsilon

Z1KinfUpsilon=Upsilon(s_eig,rho,0,1,2*K);

thisZ1Kinf=jacF_alpha;
this=@Name(Z1KinfUpsilon);

thisZ1Kinf=replace(thisZ1Kinf, "+" => " .+");

thisZ1Kinf=replace(thisZ1Kinf, "-" => " .-");

for j=1:(m+1)
    global    thisZ1Kinf=replace(thisZ1Kinf, "h_fksalphax[:,:,"*string(j)*"]" => this);
    global   thisZ1Kinf=replace(thisZ1Kinf, "h_fks0[:,:,"*string(j)*"]" => this); 
    
    
end


PiKDPhi= [1/BigFloat(2.0)^(2*K); abs.(MKinv)*eval(Meta.parse(thisZ1Kinf))];
##############################################


DR0f0=jacF;
this=@Name(Z1KinfUpsilon);


for j=1:(m+1)
    global DR0f0=replace(DR0f0, "h_fksalphax[:,:,"*string(j)*"]" => this);
    global DR0f0=replace(DR0f0, "df_fksalphax[:,"*string(j)*"]" => "abs.(df_fks0["*string(j)*"])");
    
end


boundDR0f0 = eval(Meta.parse(DR0f0));

####################################################

thistl= norm(inf.(real.(alpha)),Inf);

vartheta1(beta)= ((thistl)*(beta- 1 ./ beta) .+ sqrt.(4 .+ ((thistl)*(beta-1 ./beta)).^2))./2;

cste1_Z1Kinf(rho,beta)= Upsilon(vartheta1(beta),rho,0,1,2*K);

this=real((rho-1/rho)/abs(alpha)+sqrt(4+((rho-1/rho)/abs(alpha))^2))/2;

beta_max = inf(this);

term1(z)=func_alpha(z,m,0).*dtf(alpha.*z);  # (f^(m)(alpha * x))'

###############################################
nb=100;
sub = BigFloat.(range(0,2*pi,nb)); 

term1_coeffs= MKinv*term1(x_K)

@time est_dRf = est_PiKg_cste1([term1],cste1_Z1Kinf,K,rho,beta_max,sub,nb); 


@time est_Phi0 = est_PiKg_cste1([ffl],cste1_Z1Kinf,K,rho,beta_max,sub,nb); 


##############################################################

boundDRc0= 1/abs(alpha)^2 .*boundDR0f0.*est_Phi0 .+ 1/abs(alpha)*boundDR0f .+ 1/abs(alpha) .*est_dRf.*boundDR0f0;



PiKDPhi[2:end] =  min.(abs.(real.(PiKDPhi[2:end])), real.(boundDRc0));


Z1Kinf = transpose(weightsrho)*abs.(J)*PiKDPhi;

Z1Kinf=sup(real(Z1Kinf))

##################################################################

serialize("Z1Kinf_eig", Z1Kinf)


##########################################
##################################################################
############################### Z1inf

#fprintf('\nComputing Z1inf...\n')
cste_Z1inf(rho,alpha)= Upsilon(rho,alpha,1,0,2*K);

alpha_Z1inf= Interval(BigFloat(2.373604723386412));

Z1inf_ffl = cste_Z1inf(rho,alpha_Z1inf)*norm_C0_Enu1(ffl_coeffs,alpha_Z1inf);


Z1inf_ffl

this=Interval(BigFloat(2.373604723386412));

Z1inf_term1 = cste_Z1inf(rho,this)*norm_C0_Enu1(term1_coeffs,this);


Z1inf_term1

###########################
#####################################################

alpha=mid(alpha);
x_K=mid.(x_K);

c0=mid.(f0);

alphaf, Rf0, DRf0, = Rm_DRm_even_1Dv1(c0,m,x_K,mid.(MKinv),jacF_alpha,alpha);

D, Q = eigen(DRf0);

lambda=copy(D[end]);
lambda=real(lambda);
thislambda=Float64(lambda,RoundNearest);

M=1;
c1=reshape(Q[:,end], (K+1,1));
C=BigFloat.(zeros(K+1,M+1)) .+ BigFloat.(zeros(K+1,M+1)).*im;
C[:,1]=c0;  #fixed pt coeffs
C[:,2]=c1;  #eigen func coeffs

c1= 1.0/c1[1].*c1;

res = [c1[1]-1.0; DRf0*c1-lambda*c1]; #F(f)

Jac = [0.0 1.0 BigFloat.(zeros(1,K)) .+ BigFloat.(zeros(1,K)).*im;
-c1 DRf0-lambda*eye(K+1)]; #DF(f)

err = sum(abs.(res));
thiserr=Float64(err,RoundNearest)

for it=0:4 # Newton's iterations
    println(it);

    thisd= - Jac\res;
    
    
    global lambda=lambda + thisd[1];
    global c1=c1 + thisd[2:end];

    global res = [c1[1]-1.0; DRf0*c1-lambda*c1]; #F(f)

    
    global Jac = [0.0 1.0 BigFloat.(zeros(1,K)) .+ BigFloat.(zeros(1,K)).*im;
    -c1 DRf0-lambda*eye(K+1)]; #DF(f)
    
    global err = sum(abs.(res));
    global thiserr=Float64(err,RoundNearest)
    println(thiserr)

end


serialize("c1_2_long", c1)
serialize("lambda_2_long", lambda)
############################################################

J=inv(Jac);


f = Interval.(c1);
lambda = Interval(lambda);
J = Interval.(J);

alpha=copy(alpha0);



df = derCheb_even2odd(f); # f'
d2f = derCheb_odd2even(df); # f''

###############################################
f_coeffs = copy(f);
df_coeffs = copy(df);

f_fun(x)=eval_cheb_even(f, x);

f_num=mid.(f);
f_fun_num(x)=eval_cheb_even(f_num, x);

df_fun(x)= eval_cheb_odd(df, x);

df_fun_num(x)=eval_cheb_odd(mid.(df), x);

d2f_fun(x)= eval_cheb_even(d2f, x);

d2f_fun_num(x)= eval_cheb_even(mid.(d2f), x);

#################################
### YK
#fprintf('Computing YK...\n')
# 1st estimate

alphaf, Rf0, DRf0, = Rm_DRm_even_1Dv1i(f0,m,x_K,MKinv,jacF_alpha,alpha);

Phi = [f[1]-1.0; DRf0*f-lambda*f]; #F(f)


Jac = [Interval(0.0) Interval(1.0) Interval.(BigFloat.(zeros(1,K)) .+ BigFloat.(zeros(1,K)).*im); -f DRf0-lambda*eye(K+1)]; #DF(f)
    
########################################################

# 2nd estimate: We compute an explicit bound of |phi(alpha,f)| on an 
# ellipse, in order to get a geometrically decaying bound on each Chebyshev
# coefficient of ipi^K(phi(alpha,f))

alpha_num=mid(alpha);
lambda_num=mid(lambda);
df0f=DRf0*f;
df0f_num= mid.(df0f);

df0df = derCheb_even2odd(df0f); # f'
df0d2f = derCheb_odd2even(df0df); # f''


function phi(z)

    if isa(z[1],Complex{Interval{BigFloat}})
        
        ans= eval_cheb_even(df0f, z)-lambda*f_fun(z)
        return ans
        
    else 
        
        ans= eval_cheb_even(df0f_num, z)-lambda_num*f_fun_num(z)
        return ans
    end
        

end

function dphi(z)

    
    
    ans= eval_cheb_odd(df0df, z) - lambda*df_fun(z)
    return ans

end

function d2phi(z)

        
    ans= eval_cheb_even(df0d2f, z) - lambda*d2f_fun(z)
    return ans

end

###########################################

cste_YK(rho,alpha) = 1;

#################### 

nb=100;
sub = BigFloat.(range(0,2*pi,nb)); 

alpha_max=10*mid(rho);

phi_coeffs= MKinv*phi(x_K);

@time est_Phi = est_PiKg_cste1([phi,dphi], cste_YK,K,rho,alpha_max,sub,nb); 


Phi[2:end] = intersect.(0.0 .± sup.(abs.(real.(Phi[2:end]))), 0.0 .± sup.(est_Phi));

YK = transpose(weightsrho)*abs.(J*Phi);

YK=sup(YK[1]);

########################################


serialize("YK_eig", YK)

###################################
################################################################
####################### Yinf

#fprintf('\nComputing Yinf...\n')

function phi_inf(z)

    if isa(z[1],Complex{Interval{BigFloat}})

        
        ans= eval_cheb_even(df0f, z); 
        return ans
        
    else 
        
        ans= eval_cheb_even(df0f_num, z)   
        return ans
    end
        

end


function dphi_inf(z)        

    ans= eval_cheb_odd(df0df, z);

    return ans

end

function d2phi_inf(z)

    ans= eval_cheb_even(df0d2f, z);
    return ans

end


function cste_Yinf(rho,v)

    if isa(v,Interval{BigFloat})
        return 1/abs(lambda)*Upsilon(rho,v,1,0,2*K);
    else
        return 1/abs(lambda_num)*Upsilon(mid(rho),mid(v),1,0,2*K);        
    end

end

###################################

tol=10^(-2);
alpha_max = 9.5;

nb=1000;
sub = BigFloat.(range(0,2*pi,nb)); 

alpha_Yinf=Interval(BigFloat(8.99953226669109));

##

nb=10^4;
sub = BigFloat.(range(0,2*pi,nb)); 

phi_inf_coeffs=MKinv*phi_inf(x_K);

Yinf = cste_Yinf(rho,alpha_Yinf)*norm_C0_Enu1(phi_inf_coeffs, alpha_Yinf);

Yinf=sup(Yinf)
#######################

serialize("Yinf_eig", Yinf)


################################
##################################################################
########################## Z1KK


Z1KK_mat=Interval.(BigFloat.(eye(K+2))) - J*Jac;

this=(transpose(weightsrho)*abs.(Z1KK_mat))./transpose(weightsrho);


Z1KK=norm(sup.(norm.(this,Inf)),Inf)

##########################################

serialize("Z1KK_eig", Z1KK)

#############################################
Z1inf= 1/abs(lambda)*( 1/abs(alpha)^2*boundDR0f0*Z1inf_ffl + 1/abs(alpha)*Z1inf_f + 1/abs(alpha)*Z1inf_term1*boundDR0f0);

Z1inf=sup(Z1inf)

##########################################

serialize("Z1inf_eig", Z1inf)

#################################################
########### Z2


this=(transpose(weightsrho)*abs.(J))./transpose(weightsrho);

Z2=norm(sup.(norm.(this,Inf)),Inf) + sup(real(1/lambda))


serialize("Z2", Z2)
######################################
################################################
######################################## rstar

## Checking whether we can apply Banch fixed point theorem
Y = YK+Yinf;
Z1 = Z1KK + Z1Kinf + Z1inf;

disc= (Z1-1)^2-4*Z2*Y;

println(disc)

if disc<0
   println("No proof")
else
    println("Proof successful");
end


r1= ((1-Z1)-sqrt(disc))/(2*Z2);
r2= ((1-Z1)+sqrt(disc))/(2*Z2);

println(r1)


serialize("r1_"*string(m), r1)

tock()

println(real(mid(lambda)))


#################################


println("YK0 : $(YK0)")

println("Yinf0 : $(Yinf0)")

println("Z1KK0 : $(Z1KK0)")

println("Z1Kinf0 : $(Z1Kinf0)")

println("Z1inf0 : $(Z1inf0)")

##################################################

println("YK : $(YK)")

println("Yinf : $(Yinf)")

println("Z1KK : $(Z1KK)")

println("Z1Kinf : $(Z1Kinf)")

println("Z1inf : $(Z1inf)")

println("Z2 : $(Z2)")



