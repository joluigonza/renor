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

using PyPlot

#=
#ENV["PYTHON"] = "C:\\Users\\jorge\\anaconda3\\python.exe"

ENV["PYTHON"]= "C:\\Users\\jorge\\.julia\\conda\\3\\python.exe"

using Pkg
Pkg.add("PyCall")
Pkg.build("PyCall")
Pkg.add("PyPlot")
# make sure:  conda install matplotlib (to install matplotlib in python)
using PyCall
PyCall.python

using PyPlot
pygui(true)
plot(1:10, rand(10))

using PyPlot
println(PyPlot.matplotlib.get_backend())

PyPlot.matplotlib.use("TkAgg")  # or "Qt5Agg"
pygui(true)
plot(1:10, rand(10))

using Conda
Conda.add("tk")        # for TkAgg


using PyPlot
pygui(true)
plot(1:10, rand(10))

=#
#############################


include("general_renor_functs_final.jl")
################################

pre=2^8
setprecision(pre);
setprecision(Interval, pre);


m = 10;  # Renormalization order

d = 2; # Degree of the fixed point

ver = 3; # index of fixed point, after m=5 there is more that one

rho = Interval(BigFloat(1.9)); #Bernstein ellipse radius

rstar = BigFloat(10.0)^(-32) # rstar radius in the contraction proofs

global K0 = 40; # Chebyshev x order


###############################################

dK= Integer(maximum([ ceil(0.05*K0), 2]))


global K= copy(K0); 


# weights, with a 1 for the alpha component
ind_k = 0:2:2*K;
weightsrho=[1; 1; 2*rho.^ind_k[2:end]];


nb1=100;
sub1 = BigFloat.(range(0,2*pi,nb1)); #subdivisions lv1

nb2=500;
sub2 = BigFloat.(range(0,2*pi,nb2)); #subdivisions lv2

nb3=1000; 
sub3 = BigFloat.(range(0,2*pi,nb3)); #subdivisions lv3

#########################################

theta_K, MK, MKinv= cheb_nodes_even(K);
x_K = cos.(theta_K);


global jacF


thisdoc=open("thisJacobian"*string(m)) 
thisstring=read(thisdoc,String);
thisJac=replace(thisstring,"\n"=>"");
thisJac=replace(thisJac,";"=>"");
close(thisdoc)
jacF=max_jac_formula0(thisJac,m);

#####################################################################

eval0 = (-1).^(0:K);
eval0[2:end] = 2*eval0[2:end];
eval0=transpose(eval0);

f= deserialize("f0"*"_m"*string(m)*"_d"*string(d)*"_v"*string(ver));


#f = load("f0"*"_m"*string(m)*"_d"*string(d)*".jld2", "f");

# if length(f) > K
#     println("K0 should be at least :",length(f))
# end


f=[f; BigFloat.(zeros(K+1-length(f),1))]

thisf0=eval0*f;


#alpha= deserialize("alpha0"*"_"*string(m)*"_"*string(ver));

alpha = 0;
for k = 1:m
    alpha = eval_cheb_even(f, alpha);
end

alpha0, Rf, DRf, DRalpha = Rm_DRm_even_1D0(f,m,x_K,MKinv,jacF,alpha);
res = [thisf0[1]-1.0; alpha*f-Rf]; #F(f)

Jac = [0 eval0; DRalpha alpha*eye(K+1) - DRf]; #DF(f)

err = sum(abs.(res));
thiserr=Float64(err,RoundNearest)

for it=0:4 # Newton's iterations for fun
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

############################


for k=1:m
  
    alpha0, Rf=Rm_even_1D(f,k,x_K,MKinv);    
    res = alpha0.*f - Rf;
    
    err = norm(res,Inf);
    thiserr=Float64(err,RoundNearest);
    println(thiserr)

end


##############################################################
#############################################################################


tick()

##########################################
J=inv(Jac);


f = Interval.(f);
alpha = Interval(alpha);
J = Interval.(J);

##############################

D, Q = eigen(Real.(Jac));

####################################################
##############################  YK
#interval version of 
# theta_K, MK, MKinv= cheb_nodes_even(K);
# x_K = cos.(theta_K);

ipi= @interval pi 


theta_K = (K.-(0:K)).*(ipi/BigFloat(2*K)); #In the even case, we only need the nodes in [0,1]
ind_k = 0:2:2*K;
x_K = cos.(theta_K);

# Construction of the matrix to go from values at Chebyshev nodes to
# Chebyshev coefficients
MKinv = transpose(cos.(theta_K*transpose(ind_k)))/K; ######
MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]]/2;
MKinv[K+1,:] = MKinv[K+1,:]/2;


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

############################# YK
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

thist = Interval(maximum([maximum(sup.(abs.(fksalphax))),1.0]));

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

######################################################



function phi(z)

    if isa(z[1],Complex{Interval{BigFloat}})

        
        this=alpha.*z;
        
        for i=1:m
            
            this=f_fun(this);

        end    
        
        ans=alpha*f_fun(z) - this;
        return ans
        
    else 
        this=mid(alpha).*z;
        
        for i=1:m
            
            this=f_fun_num(this);
        end

        ans=mid(alpha)*f_fun_num(z)-this;
        return ans
    end
        

end

function dphi(z)

    this=alpha*z;
    that=df_fun(this);

    for i=2:m
        this=f_fun(this);
        that=df_fun(this).*that;        
    end

    ans=alpha*(df_fun(z) - that);

    return ans

end


thisdoc=open("d2phi"*string(m)) 

thisstring=read(thisdoc,String);
close(thisdoc)


this_d2phi=replace(thisstring,"\n"=>"");
this_d2phi=replace(this_d2phi,";"=>"");


this_d2phi=replace(this_d2phi,"pderivop(f,1)" =>"df_fun");
this_d2phi=replace(this_d2phi,"pderivop(f,2)" =>"d2f_fun");
this_d2phi=replace(this_d2phi,"c" =>"alpha");
this_d2phi=replace(this_d2phi,"f(" =>"f_fun(");
this_d2phi=replace(this_d2phi,"z" =>"thisz");
this_d2phi=replace(this_d2phi,"*" =>" .*");
this_d2phi=replace(this_d2phi,"^" =>".^");


macro Name(arg)
    string(arg);
end

function d2phi(z,this_d2phi) 
    
    global this
    this=z;
    thatz=@Name(this);
    foo=replace(this_d2phi,"thisz" =>"$thatz");

    ans= alpha*d2f_fun(this) - eval(Meta.parse(foo));

    return ans
end

d2phiz(z)=d2phi(z,this_d2phi);

############### alpha_Yinf precomputed for smaller K
tol=10^(-2);

alpha_max = 9;

###############

@time alpha_Yinf= opt_para(ffl, mid(rho), cste_Yinf, alpha_max, tol,sub3); 

alpha_Yinf = Interval(BigFloat(alpha_Yinf));


#################################################


@time temp = compute_sup_E([ffl,d_ffl],alpha_Yinf,sub3,nb3)

sup.(temp) ./ (inf.(temp).+ 1.0)

thisconst= cste_Yinf(rho,alpha_Yinf)


global Yinf = thisconst*temp

Y_inf1=0.0;
foo=0.0
K1= copy(K0) 



if sup(Yinf)>0.5*rstar
    
    while foo < rstar*0.1
    Y_inf1=copy(foo)
    global K1 = K1 + dK
    theta_K1, MKinv1, weightsrho1 = cheb_nodes_even_i(K1);
    x_K1= cos.(theta_K1);

    phi_coeffs0= Interval.(BigFloat.(zeros(K1+1,1)).+ BigFloat.(zeros(K1+1,1)).*im);
    phi_coeffs0[1:K0+1]=MKinv*phi(x_K);


    phi_coeffs1=MKinv1*phi(x_K1);

    foo= 1/abs(alpha)*transpose(weightsrho1[2:end])*abs.(phi_coeffs1-phi_coeffs0);
    foo=sup(foo[1])
    end

    ##########################
    global K= copy(K1)-dK; # Chebyshev x order        

    thisconst= cste_Yinf(rho,alpha_Yinf)

    subtol= 0.5*10.0^(ceil(log(10,sup(rstar*10.0^(-1)/(thisconst)))))

    temp= compute_sup_E1_adap([ffl,d_ffl],alpha_Yinf,subtol)
    
    global Yinf = thisconst*temp

    global K=copy(K0);
end



###############################################

sup(Yinf)

sup(Y_inf1[1])

Yinf=sup(Yinf) + sup(Y_inf1[1])



if Yinf>rstar
    error("No proof: Yinf is greater than rstar")
end
###################################

Yinf0= copy(Yinf)
serialize("Yinf0", Yinf0)

###################################################################
########################## Z1KK
### Creating \tilde{f} for Z1
#rstar = BigFloat(10.0)^(-32)



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

thist = Interval(maximum([maximum(sup.(abs.(fksalphax))),1.0]));

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

PiKDPhi = [0; abs.((MKinv))*eval(Meta.parse(thisZ1Kinf))]



##############################################################
####################################

func_alpham(z)=func_alpha(z,m,0).*dtf(tl.*z).*z;

# 2nd estimate

thistl= norm(inf.(real.(tl)),Inf);

vartheta1(beta)= ((thistl)*(beta- 1 ./ beta) .+ sqrt.(4 .+ ((thistl)*(beta-1 ./beta)).^2))./2;


cste1_Z1Kinf(rho,beta)= Upsilon(vartheta1(beta),rho,0,1,2*K);

this=real((rho-1/rho)/abs(tl)+sqrt(4+((rho-1/rho)/abs(tl))^2))/2

beta_max = inf(this)
beta_max= minimum([beta_max, alpha_max])

term1(z)=func_alpha(z,m,0);  

###########################


@time est_PiKDPhi_1, temp = est_PiKg_cste1([term1],cste1_Z1Kinf,K,rho,beta_max,sub3,nb3) 


subtol1 = (sup.(temp) ./ (inf.(temp).+ 1.0))


subtol1 = tail_func.(subtol1,10.0,10.0,maximum([0.05*maximum(subtol1), 10.0]))


@time est_PiKDPhi_1 = est_PiKg_cste_adap([term1],cste1_Z1Kinf,K,rho,beta_max,subtol1)


gammas=interval.(BigFloat.(zeros(m,1)));
gammas[1]=beta_max;
thissum=est_PiKDPhi_1

###########################################


PiKDPhi[2:end] =  min.(real.(PiKDPhi[2:end]), thissum);

boundDR0f= PiKDPhi[2:end]; # save for eigenvalue proof


Z1Kinf = transpose(weightsrho)*abs.(J)*PiKDPhi

#########################################


Z1Kinf=sup(real(Z1Kinf))

##########################

Z1Kinf0= copy(Z1Kinf)

serialize("Z1Kinf0", Z1Kinf0)

########################################
##################################################################
############################### Z1inf


### Z1inf
#fprintf('\nComputing Z1inf...\n')


cste_Z1inf(rho,alpha)= Upsilon(rho,alpha,1,0,2*K);
 
##############################
tol=10^(-2);

###############################################

nb=500; 
sub = BigFloat.(range(0,2*pi,nb)); 

@time alpha_Z1inf, Z1inf_alpha =opt_para_glob([func_alpham], mid(rho), cste_Z1inf, 4.0,tol, sub,nb)

#############################################
if sup(Z1inf_alpha) > 0.5

    thisconst= cste_Z1inf(rho,this)

    subtol= 10.0^(ceil(log(10,sup(10.0^(-1)/(thisconst)))))

    #compute_sup_E([func_alpham],this,sub3,nb3) #keep here 

    @time temp= compute_sup_E1_adap([func_alpham],this,subtol)

    Z1inf_alpha = thisconst*temp

end
###########################


betas=BigFloat.(zeros(m,1));
betas[1]=minimum([beta_max, alpha_max]);

thissum=BigFloat(0.0);

for n=0:m-2
    println(n)
    thisfun(z)=tffl(z,n+1);
    
    thisbeta=opt_incl(thisfun, inf(rho),alpha_max,tol,sub3,nb3)


    betas[n+2]=thisbeta;
    thisbeta=Interval(thisbeta);

    k=m-n;
    
    thisterm(z)= func_alpha(z,k,n);   


    foo1= opt_para(thisterm, mid(rho), cste_Z1inf, betas[n+1],tol,sub3);    
    foo=Interval(foo1[1])


    thisconst= cste_Z1inf(rho,foo)

    #@time temp= compute_sup_E([thisterm],foo,sub3,nb3) #keep here
    
    subtol= 0.5*10.0^(ceil(log(10,sup(10.0^(-1)/(thisconst)))))
        
    @time temp= compute_sup_E1_adap([thisterm],foo,subtol)
    
    thisbound = thisconst*temp


    global thissum=thissum+thisbound

end

#add last one h(f(m)) at the end
thissum=thissum + Upsilon(rho,Interval(betas[end]),1,0,2*K);

Z1inf_f= thissum; # save for eigenvalue proof

Z1inf= 1/abs(tl) * ( max(Z1inf_alpha, Z1inf_f ) + tr);

Z1inf=sup(Z1inf)

##############################################

Z1inf0= copy(Z1inf)

serialize("Z1inf0", Z1inf0)

############################################
################################################
######################################## rstar

### Checking whether we can apply Banach fixed point theorem
Y = YK+Yinf;
Z1 = Z1KK + Z1Kinf + Z1inf;
rmin = sup(Y/(1-Z1));

if Z1>1
   error("negative Z1>1 not allowed")
end

if rmin < rstar
    println("Proof successful");
    
else  
    println("No proof");

end

println(rmin)
println(rstar)

###############################

rmin0= copy(rmin)

serialize("rmin0_"*string(m), rmin0)

#########################################################################
############################ start eigen ###############################
############################################################################

f0=copy(f);
alpha0=copy(alpha);

global jacF_alpha

jacF_alpha= "";


alphaf, Rf0, DRf0, = Rm_DRm_even_1Dv1i(f0,m,x_K,MKinv,jacF_alpha,alpha);


###########################################################################
################################ Z1Kinf

#fprintf('\nComputing Z1Kinf...\n')

######################################################
global Z1KinfUpsilon

Z1KinfUpsilon=Upsilon(s_eig,rho,0,1,2*K);

thisZ1Kinf0_1D0=num_Z1Kinf0_1D0(m,Z1KinfUpsilon);

thisZ1Kinf_1D0f= num_Z1Kinf_1D0(m,Z1KinfUpsilon);

~, forterm1 = der_f_term(m,0);

thisterm1=forterm1 .* df_fksalphax[:,1];


PiKDPhi= [1/BigFloat(mid(rho))^(2*K); abs.(MKinv)*num_DRf_even_1Dv1(m, thisZ1Kinf0_1D0, thisZ1Kinf_1D0f, thisterm1)];

##############################################
############################################################


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

beta_max = minimum([inf(this), alpha_max])


#####################################################

term1(z)=func_alpha(z,m,0).*dtf(alpha.*z);  # (f^(m)(alpha * x))'

###############################################


@time est_dRf, temp = est_PiKg_cste1([term1],cste1_Z1Kinf,K,rho,beta_max,sub3,nb3)

subtol1 = (sup.(temp) ./ (inf.(temp).+ 1.0))

subtol1 = tail_func.(subtol1,10.0,10.0,maximum([0.05*maximum(subtol1), 10.0]))

@time est_dRf = est_PiKg_cste_adap([term1],cste1_Z1Kinf,K,rho,beta_max,subtol1)

##################


@time est_Phi0, temp = est_PiKg_cste1([ffl],cste1_Z1Kinf,K,rho,beta_max,sub3,nb3)

subtol1 = (sup.(temp) ./ (inf.(temp).+ 1.0)) 

subtol1 = tail_func.(subtol1,10.0,10.0,maximum([0.05*maximum(subtol1), 10.0]))


@time est_Phi0 = est_PiKg_cste_adap([ffl],cste1_Z1Kinf,K,rho,beta_max,subtol1)



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

##########################################


nb=500; 
sub = BigFloat.(range(0,2*pi,nb)); 

@time alpha_Z1inf, Z1inf_ffl =opt_para_glob([ffl], mid(rho), cste_Z1inf, 4.0,tol, sub,nb)

this=Interval(alpha_Z1inf[1])

###############################################

if sup(Z1inf_ffl) > 0.5

    thisconst= cste_Z1inf(rho,this)

    #compute_sup_E([ffl],this,sub3,nb3)  #keep here 


    subtol= 10.0^(ceil(log(10,sup(10.0^(-1)/(thisconst)))))


    @time temp= compute_sup_E1_adap([ffl],this,subtol)


    Z1inf_ffl = thisconst*temp


end

###############################


nb=500; 
sub = BigFloat.(range(0,2*pi,nb)); 

@time alpha_Z1inf, Z1inf_term1 =opt_para_glob([term1], mid(rho), cste_Z1inf, 4.0,tol, sub,nb)

this=Interval(alpha_Z1inf[1])

######################################
if sup(Z1inf_term1)>0.5

    thisconst= cste_Z1inf(rho,this)

    subtol= 10.0^(ceil(log(10,sup(10.0^(-1)/(thisconst)))))

    #compute_sup_E([term1],this,sub3,nb3)  #keep here 

    @time temp= compute_sup_E1_adap([term1],this,subtol)

    Z1inf_term1 = thisconst*temp

end
###########################
#####################################################

alpha=mid(alpha);
x_K=mid.(x_K);

c0=mid.(f0);

alphaf, Rf0, DRf0, = Rm_DRm_even_1Dv1(c0,m,x_K,mid.(MKinv),jacF_alpha,alpha);

#c1=deserialize("c1"*"_"*string(m)*"_"*string(ver));
#c1=[c1; BigFloat.(zeros(K+1-length(c1),1))]

D, Q = eigen(DRf0);

lambda=copy(D[end]);
lambda=real(lambda);
#thislambda=Float64(lambda,RoundNearest);

c1=reshape(Q[:,end], (K+1,1));


#lambda=deserialize("lambda"*"_"*string(m)*"_"*string(ver));

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


alpha_max=9.0;


@time est_Phi, temp = est_PiKg_cste1([phi,dphi], cste_YK,K,rho,alpha_max,sub3,nb3) 


Phi[2:end] = intersect.(0.0 .± sup.(abs.(real.(Phi[2:end]))), 0.0 .± sup.(est_Phi));

YK = transpose(weightsrho)*abs.(J*Phi);

YK=sup(YK[1])

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



########################################

@time alpha_Yinf = opt_para(phi_inf, mid(rho), cste_Yinf, alpha_max, tol,sub3);

alpha_Yinf = Interval(BigFloat(alpha_Yinf))

###################################################

# nb=500; 
# sub = BigFloat.(range(0,2*pi,nb)); 

# @time alpha_Yinf, Yinf =opt_para_glob([phi_inf], mid(rho), cste_Yinf, 4.0,tol, sub,nb)


################################################


@time temp=compute_sup_E([phi_inf,dphi_inf],alpha_Yinf,sub3,nb3)

(sup.(temp) ./ (inf.(temp).+ 1.0)) 
##################################################

thisconst= cste_Yinf(rho,alpha_Yinf)


global Yinf = thisconst*temp

Y_inf1=0.0;
foo=0.0;
K1= copy(K0)

###################################################################

if sup(Yinf)>0.5*rstar
    
    while foo < rstar*0.1
    Y_inf1=copy(foo)
    global K1 = K1 + dK
    theta_K1, MKinv1, weightsrho1 = cheb_nodes_even_i(K1);
    x_K1= cos.(theta_K1);

    phi_coeffs0= Interval.(BigFloat.(zeros(K1+1,1)).+ BigFloat.(zeros(K1+1,1)).*im);
    phi_coeffs0[1:K0+1]=MKinv*phi(x_K);


    phi_coeffs1=MKinv1*phi(x_K1);

    foo= 1/abs(alpha)*transpose(weightsrho1[2:end])*abs.(phi_coeffs1-phi_coeffs0)
    foo=sup(foo[1])

    end 
    ###################################################
    global K= copy(K1) - dK; # Chebyshev x order

    thisconst = cste_Yinf(rho,alpha_Yinf)
    
    #compute_sup_E([phi_inf,dphi_inf],alpha_Yinf, sub3,nb3)

    subtol= 0.5*10.0^(ceil(log(10,sup(rstar*10.0^(-1)/(thisconst)))))

    @time temp = compute_sup_E1_adap([phi_inf,dphi_inf],alpha_Yinf,subtol);
    
    
    global Yinf = thisconst*temp

    global K=copy(K0);
end


############################################

sup(Yinf)

sup(Y_inf1[1])

Yinf=sup(Yinf) + sup(Y_inf1[1])

if Yinf>rstar
error("No proof: Yinf is greater than rstar")
end
#######################

serialize("Yinf_eig", Yinf)

######################### Z1KK

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
   error("No proof")
else
    println("Proof successful");
end



r1= ((1-Z1)-sqrt(disc))/(2*Z2);
r2= ((1-Z1)+sqrt(disc))/(2*Z2);

if r2<0
    error("r2 negative not allowed")
end

println(r1)


serialize("r1_"*string(m), r1)

tock()

println(real(mid(lambda)))

###################################

println(real.(mid.(f0)))


println("YK0 : $(YK0)")

println("Yinf0 : $(Yinf0)")

println("Z1KK0 : $(Z1KK0)")

println("Z1Kinf0 : $(Z1Kinf0)")

println("Z1inf0 : $(Z1inf0)")

println("rmin0 : $(rmin0)")

##################################################

println(real(mid(lambda)))

println("YK : $(YK)")

println("Yinf : $(Yinf)")

println("Z1KK : $(Z1KK)")

println("Z1Kinf : $(Z1Kinf)")

println("Z1inf : $(Z1inf)")

println("Z2 : $(Z2)")

println("r1 : $(r1)")

