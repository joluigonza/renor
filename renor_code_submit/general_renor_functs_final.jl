

#############################################
function cheb_nodes_even(K)
    #pre=precision(f[1])
    #setprecision(pre);

    pre_fac=pi/BigFloat(2*K);
    theta_K = (K.-(0:K)).*pre_fac;
    x_K = cos.(theta_K);

    # Chebyshev nodes x_0,...,x_K
    # By symmetry, it is enough to only keep the nodes that are in [0,pi])

    # Construction of the matrix to go from Chebyshev coefficients to values at
    # Chebyshev nodes
    MK = cos.(theta_K*transpose(0:2:2*K));
    MK[:,2:end] = 2*MK[:,2:end];

    # Construction of the matrix to go from values at Chebyshev nodes to
    # Chebyshev coefficients
    MKinv = transpose(MK) / (2*K);
    MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]] / 2;
    MKinv[K+1,:] = MKinv[K+1,:] / 2;
    MKinv[1,:] = MKinv[1,:] * 2;

    return theta_K, MK, MKinv

end
#######################################################################

function cheb_nodes_even_i(K)
        
    theta_K = (K.-(0:K)).*(ipi/BigFloat(2*K)); #In the even case, we only need the nodes in [0,1]
    ind_k = 0:2:2*K;
    #x_K = cos.(theta_K);

    # Construction of the matrix to go from values at Chebyshev nodes to
    # Chebyshev coefficients
    MKinv = transpose(cos.(theta_K*transpose(ind_k)))/K; ######
    MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]]/2;
    MKinv[K+1,:] = MKinv[K+1,:]/2;

    # weights, with a 1 for the alpha component
    weightsrho=[1; 1; 2*rho.^ind_k[2:end]];
    return theta_K, MKinv, weightsrho

end

#####################################################################
function cheb_nodes(K)
    #pre=precision(f[1])
    #setprecision(pre);

    #K=K-1;
    pre_fac=pi/BigFloat(K);
    theta_K = (K.-(0:K)).*pre_fac;
    x_K = cos.(theta_K);

    # Chebyshev nodes x_0,...,x_K
    # By symmetry, it is enough to only keep the nodes that are in [0,pi])

    # Construction of the matrix to go from Chebyshev coefficients to values at
    # Chebyshev nodes
    MK = cos.(theta_K*transpose(0:K));
    #MK[:,2:end] = 2*MK[:,2:end];

    # Construction of the matrix to go from values at Chebyshev nodes to
    # Chebyshev coefficients
    MKinv = transpose(MK) / (K);
    MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]] / 2;
    MKinv[K+1,:] = MKinv[K+1,:] / 2;
    #MKinv[1,:] = MKinv[1,:] * 2;

    return theta_K, MK, MKinv

end

###############################################
##################################

function eval_cheb_even(f,x)

    ff=copy(f);
    thisx=copy(x);
    # Computation of f_0 + 2*sum_{k=1}^K f_k T_k(x) (default convention), where
    # f is a vector containing the coefficients f_k. If x is a vector, this 
    # formula is applied component-wise to x. If type = 'even' (resp. type = 
    # 'odd'), we assume that f contains only even (resp) odd coefficients, that
    # is we instead compute f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x) (resp. 
    # 2*sum_{k=0}^K f_{2k+1} T_{2k+1}(x)).
    #
    # Here, T_k(x) is computed using the Clenshaw algorithm
    
    K = length(f)-1;
    if size(x,1) == 1
        thisx = transpose(thisx);
        
    end
    
    # if nargin == 2
    #     type = 'none';
    # end
                    
       
    #ff = copy(thisf);
    ff[2:end] = 2*ff[2:end];
    #thisf = zeros(2*K+1,1);
    thisf=0 .*[f;f[1:end-1]]
    
    thisf[1:2:end] = ff;
    K = 2*K;
        
        
    ## Clenshaw algorithm
    if K == 0
        u = thisf*ones(size(x));
    else
        uuu = 0;
        uu = thisf[end];
        for k = K:-1:2
            u = 2*x.*uu .- uuu .+ thisf[k];
            uuu = uu;
            uu = u;
        end
        u = x.*uu .- uuu .+ thisf[1];
    end

    return u
end

###############################################

function eval_cheb_odd(f,x)

    ff=copy(f);
    thisx=copy(x);
    # Computation of f_0 + 2*sum_{k=1}^K f_k T_k(x) (default convention), where
    # f is a vector containing the coefficients f_k. If x is a vector, this 
    # formula is applied component-wise to x. If type = 'even' (resp. type = 
    # 'odd'), we assume that f contains only even (resp) odd coefficients, that
    # is we instead compute f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x) (resp. 
    # 2*sum_{k=0}^K f_{2k+1} T_{2k+1}(x)).
    #
    # Here, T_k(x) is computed using the Clenshaw algorithm
    
    K = length(f)-1;
    if size(x,1) == 1
        thisx = transpose(thisx);
        
    end
    
    # if nargin == 2
    #     type = 'none';
    # end
                    
       
    #ff = copy(thisf);
    ff = 2*ff;
    thisf = 0 .*[f;f];
    
    thisf[2:2:end] = ff;
    K = 2*K+1;
        
        
    ## Clenshaw algorithm
    if K == 0
        u = thisf*ones(size(x));
    else
        uuu = 0;
        uu = thisf[end];
        for k = K:-1:2
            u = 2*x.*uu .- uuu .+ thisf[k];
            uuu = uu;
            uu = u;
        end
        u = x.*uu .- uuu .+ thisf[1];
    end

    return u
end

##############################################

function eval_cheb(f,x)

    thisf=copy(f);
    thisx=copy(x);
    # Computation of f_0 + 2*sum_{k=1}^K f_k T_k(x) (default convention), where
    # f is a vector containing the coefficients f_k. If x is a vector, this 
    # formula is applied component-wise to x. If type = 'even' (resp. type = 
    # 'odd'), we assume that f contains only even (resp) odd coefficients, that
    # is we instead compute f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x) (resp. 
    # 2*sum_{k=0}^K f_{2k+1} T_{2k+1}(x)).
    #
    # Here, T_k(x) is computed using the Clenshaw algorithm
    
    K = length(f)-1;
    if size(x,1) == 1
        thisx = transpose(thisx);
        
    end
    
    # if nargin == 2
    #     type = 'none';
    # end
                    
       
    #ff = copy(thisf);
    thisf[2:end] = 2*thisf[2:end];
            
        
    ## Clenshaw algorithm
    if K == 0
        u = thisf*ones(size(x));
    else
        uuu = 0;
        uu = thisf[end];
        for k = K:-1:2
            u = 2*x.*uu .- uuu .+ thisf[k];
            uuu = uu;
            uu = u;
        end
        u = x.*uu .- uuu .+ thisf[1];
    end

    return u
end

##################################

function derCheb_even2odd(u)
    # pre=precision(u[1])
    # setprecision(pre);
    # The input u is supposed to be a column. The convention is that
    # u = u_0 + 2*\sum u_{2k} T_{2k}. The output v contains the Chebyshev
    # coefficients of u' (only those with odd indices, since the other are 0)


    K=length(u)-1;
    v=reverse((2:2:2*K).*u[2:end]);

    v = 2*cumsum(v);
    v=reverse(v);
    return v
end
#######################################

function derCheb_odd2even(u)
    # pre=precision(u[1])
    # setprecision(pre);
    # The input u is supposed to be a column. The convention is that
    # u =  2*\sum u_{2k+1} T_{2k+1}. The output v contains the Chebyshev
    # coefficients of u' (only those with even indices, since the other are 0)


    K=length(u)-1;
    v=reverse((1:2:2*K+1).*u[1:end]);

    v = 2*cumsum(v);
    v=reverse(v);
    return v
end
#####################################

function derCheb(u)
    # pre=precision(u[1])
    # setprecision(pre);
    # The input u is supposed to be a column. The convention is that
    # u =  2*\sum u_{2k+1} T_{2k+1}. The output v contains the Chebyshev
    # coefficients of u' (only those with even indices, since the other are 0)


    K=length(u)-1;

#############################
   

#################

    derMat=BigFloat.(zeros(K+1,K+1));
    derMat[1,2:2:end]=2*(1:2:K);
    for k=2:K
    derMat[k,k+1:2:end]=2*(k:2:K)
    end
v=derMat*u;
v=v[1:end-1];

    return v
end

#########################################
function R2_DR2_even_1D(f)

    # f is a vector of Chebyshev coefficients f_0,f_2,...,f_{2K}, with the
    # normalization f(x) = f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x).

    # Outputs the Chebyshev coefficients of R(f), and the jacobian matrix DR(f)
    # for m = 2, i.e. for R(f)(x) = 1/lambda * f(f( lambda*x )), where
    # lambda = [f(f(0))]
    # pre=precision(f[1]);
    # setprecision(pre);

    K = length(f)-1;


    pre_fac=pi/BigFloat(2*K);
    theta_K = (K.-(0:K)).*pre_fac;
    x_K = cos.(theta_K);
    MK = cos.(theta_K*transpose(0:2:2*K));

    # Chebyshev nodes x_0,...,x_K
    # By symmetry, it is enough to only keep the nodes that are in [0,pi])

    # Construction of the matrix to go from Chebyshev coefficients to values at
    # Chebyshev nodes

    MK[:,2:end] = 2*MK[:,2:end];

    # Construction of the matrix to go from values at Chebyshev nodes to
    # Chebyshev coefficients
    MKinv = transpose(MK) / (2*K);
    MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]] / 2;
    MKinv[K+1,:] = MKinv[K+1,:] / 2;
    MKinv[1,:] = MKinv[1,:] * 2;


    f0 = eval_cheb_even(f,BigFloat(0.0)); # f(0)
    
    lambda = eval_cheb_even(f, f0); # lambda=f(f(0))
    lambdax = lambda.*x_K; # f(f(0))*x_k, k=0,...,K
    flambdax = eval_cheb_even(f, lambdax); # f(lambda*x_k), k=0,...,K
    fflambdax = eval_cheb_even(f, flambdax); # f(f(lambda*x_k)), k=0,...,K


    # R
   
    R = MKinv*((1 ./ lambda).*fflambdax);

    # DR

    df = derCheb_even2odd(f); # f'
    dff0 = eval_cheb_odd(df, f0); # df(f(0))
    dflambdax = eval_cheb_odd(df, lambdax); # f'(lambda*x_k), k=0,...,K
    dfflambdax = eval_cheb_odd(df, flambdax); # f'(f(lambda*x_k)), k=0,...,K

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the value h(0).
    eval0 = (-1).^transpose(0:K);
    eval0[2:end] = 2*eval0[2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the value h(f(0)).
    evalf0 = cos.(myacos.(f0)*transpose(0:2:2*K));
    evalf0[2:end] = 2*evalf0[2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the values h(lambda*x_k).
    MKlambda = cos.(myacos.(lambdax)*transpose(0:2:2*K));
    MKlambda[:,2:end] = 2*MKlambda[:,2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the values h(f(lambda*x_k)).
    MKflambda = cos.(myacos.(flambdax)*transpose(0:2:2*K));
    MKflambda[:,2:end] = 2*MKflambda[:,2:end];

    DR= MKinv*((-1 ./ lambda^2).*fflambdax*(dff0*eval0+evalf0)
                +(1 ./ lambda).*(dfflambdax.*dflambdax.*x_K)*(evalf0+dff0*eval0)
                +(1 ./ lambda).*Diagonal(dfflambdax)*MKlambda
                +(1 ./ lambda).*MKflambda);

    

    return R, DR, x_K, MK, MKinv
end


############################

function RmP_even(c,m,x_K,MKinv)


    alpha=BigFloat(0.0);

    for i=1:m

        alpha=eval_cheb_even(c,alpha);

    end

    Rf=alpha.*x_K;

    for i=1:m

        Rf=eval_cheb_even(c, Rf);

    end

    Rf=MKinv*((1 ./ alpha).*Rf);

    return Rf

end

##########################################################################

function R0mP_even(c,m,x_K,MKinv)


    alpha=BigFloat(0.0);

    for i=1:m

        alpha=eval_cheb_even(c,alpha);

    end

    Rf=alpha.*x_K;

    for i=1:m

        Rf=eval_cheb_even(c, Rf);

    end

    Rf=MKinv*((1 ./ alpha).*Rf);

    return Rf,alpha

end

###########################################################################


function RmP_alpha_even(c,m,x_K,MKinv,alpha)

    
    Rf=alpha.*x_K;

    for i=1:m

        Rf=eval_cheb_even(c, Rf);

    end

    Rf=MKinv*((1 ./ alpha).*Rf);

    return Rf

end


############################################################

function R3_DR3_even_1D(f)


    # f is a vector of Chebyshev coefficients f_0,f_2...,f_{2K}, with the
    # normalization f(x) = f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x).

    # Outputs the Chebyshev coefficients of R(f), and the jacobian matrix DR(f)
    # for m = 3, i.e. for R(f)(x) = 1/f(f(f(0))) * f(f(f( f(f(f(0)))x )))

    K = length(f)-1;

    pre_fac=pi/BigFloat(2*K);
    theta_K = (K.-(0:K)).*pre_fac;
    x_K = cos.(theta_K);
    MK = cos.(theta_K*transpose(0:2:2*K));

    # Construction of the matrix to go from Chebyshev coefficients to values at
    # Chebyshev nodes

    MK[:,2:end] = 2*MK[:,2:end];

    # Construction of the matrix to go from values at Chebyshev nodes to
    # Chebyshev coefficients
    MKinv = transpose(MK) / (2*K);
    MKinv[:,[1,K+1]] = MKinv[:,[1,K+1]] / 2;
    MKinv[K+1,:] = MKinv[K+1,:] / 2;
    MKinv[1,:] = MKinv[1,:] * 2;


    f0 = eval_cheb_even(f, BigFloat(0.0)); # f(0)
    ff0 = eval_cheb_even(f, f0); # ff(0)
    lambda = eval_cheb_even(f, ff0); # lambda=f(f(f(0)))
    lambdax = lambda.*x_K; # f(f(1))*x_k, k=0,...,K
    flambdax = eval_cheb_even(f, lambdax); # f(lambda*x_k), k=0,...,K
    fflambdax = eval_cheb_even(f, flambdax); # f(f(lambda*x_k)), k=0,...,K
    ffflambdax = eval_cheb_even(f, fflambdax); # f(f(f(lambda*x_k))), k=0,...,K


    R = MKinv*((1 ./ lambda).*ffflambdax);
    
    df = derCheb_even2odd(f); # f'
    dff0 = eval_cheb_odd(df, f0); # df(f(0))
    dfff0 = eval_cheb_odd(df, ff0); # df(f(f(0)))
    dflambdax = eval_cheb_odd(df, lambdax); # f'(lambda*x_k), k=0,...,K
    dfflambdax = eval_cheb_odd(df, flambdax); # f'(f(lambda*x_k)), k=0,...,K
    dffflambdax = eval_cheb_odd(df, fflambdax); # f'(f(f(lambda*x_k))), k=0,...,K

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the value h(0).
    eval0 = (-1).^transpose(0:K);
    eval0[2:end] = 2*eval0[2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the value h(f(0)).
    evalf0 = cos.(myacos.(f0)*transpose(0:2:2*K));
    evalf0[2:end] = 2*evalf0[2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the value h(f(f(0))).
    evalff0 = cos.(myacos.(ff0)*transpose(0:2:2*K));
    evalff0[2:end] = 2*evalff0[2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the values h(lambda*x_k).
    MKlambda = cos.(myacos.(lambdax)*transpose(0:2:2*K));
    MKlambda[:,2:end] = 2*MKlambda[:,2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the values h(f(lambda*x_k)).
    MKflambda = cos.(myacos.(flambdax)*transpose(0:2:2*K));
    MKflambda[:,2:end] = 2*MKflambda[:,2:end];

    # Construction of the matrix to go from the Chebyshev coefficients h_{2k}
    # to the values h(f(f(lambda*x_k))).
    MKfflambda = cos.(myacos.(fflambdax)*transpose(0:2:2*K));
    MKfflambda[:,2:end] = 2*MKfflambda[:,2:end];

    DR= MKinv*((-1 ./ lambda^2).*ffflambdax*(dfff0.*dff0.*eval0+dfff0.*evalf0+evalff0)
                +(1 ./ lambda).*(dffflambdax.*dfflambdax.*dflambdax.*x_K)*(dfff0.*dff0.*eval0+dfff0.*evalf0+evalff0)
                +(1 ./ lambda).*Diagonal(dffflambdax)*(Diagonal(dfflambdax)*MKlambda+MKflambda)
                +(1 ./ lambda).*MKfflambda);

    return R, DR, x_K, MK, MKinv
    
end



#########################################################

function max_jac_formula(thisJac,m)
    
    ###################################
    gks0_string=[];
    gen="0";
    gks0_string=cat(gks0_string,gen,dims=1);
    for k=1:m
       gen="g0("*gen*")";
       gks0_string=cat(gks0_string,gen,dims=1);
    end

    ##################################

    thisJac=replace(thisJac, gen => "alpha");
    thisJac=replace(thisJac,"alpha*x" => "alphax");
    thisJac=replace(thisJac,"pderivop(g0,1)" =>"dg");
    thisJac=replace(thisJac,"x*" =>"x_K*");
    thisJac=replace(thisJac,"/"=>"*1/");
    thisJac=replace(thisJac,"*"=>".*");

    h_gks0_string="g1(".*gks0_string.*")";
    dg_gks0_string="dg(".*gks0_string.*")";

    ####################################
    gksalphax_string=[];
    gen="alphax";
    gksalphax_string=cat(gksalphax_string,gen,dims=1);

    for k=1:m
       gen="g0("*gen*")";
       gksalphax_string=cat(gksalphax_string,gen,dims=1);
    end

    h_gksalphax_string="g1(".*gksalphax_string.*")";
    dg_gksalphax_string="dg(".*gksalphax_string.*")";
    brunos = [h_gks0_string, h_gksalphax_string, dg_gks0_string, dg_gksalphax_string];


    thisJac=replace(thisJac,gksalphax_string[end] =>"fksalphax[:,"*string(m+1)*"]");

    subvars=["h_fks0[:,:,", "h_fksalphax[:,:,", "df_fks0[",  "df_fksalphax[:,"];


    for i=1:4
        thisbruno=brunos[i];
        for k=1:m+1
        j=m+1-k+1;
        
        thisJac=replace(thisJac,thisbruno[j] => subvars[i]*string(j)*"]");
        end
    end

    return thisJac

end

################################################

function max_jac_formula0(thisJac,m)
    

    thisJac=replace(thisJac,"(x)" => "(alphax)");


    thisJac=replace(thisJac,"pderivop(g0,1)" =>"dg");
    thisJac=replace(thisJac,"x*" =>"x_K*");
    thisJac=replace(thisJac,"/"=>"*1/");
    thisJac=replace(thisJac,"*"=>".*");
    thisJac=replace(thisJac,"+"=>".+");


    ####################################
    gksalphax_string=[];
    gen="alphax";
    gksalphax_string=cat(gksalphax_string,gen,dims=1);

    for k=1:m
       gen="g0("*gen*")";
       gksalphax_string=cat(gksalphax_string,gen,dims=1);
    end

    h_gksalphax_string="g1(".*gksalphax_string.*")";
    dg_gksalphax_string="dg(".*gksalphax_string.*")";
    
    brunos = [h_gksalphax_string, dg_gksalphax_string];


    
    subvars=["h_fksalphax[:,:,", "df_fksalphax[:,"];


    for i=1:length(brunos)
        thisbruno=brunos[i];
        for k=1:m+1
        j=m+1-k+1;
        
        thisJac=replace(thisJac,thisbruno[j] => subvars[i]*string(j)*"]");
        end
    end

return thisJac

end



################################################


function Rm_DRm_even_1D(f,m,x_K,MKinv,jacF)


    K = length(f)-1;


    global fks0
    fks0=BigFloat.(zeros(1,m+1));

    for i=1:m
        
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alpha=fks0[end];

    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1));
    fksalphax[:,1]=alpha.*x_K;

    for i=1:m
        
        fksalphax[:,i+1]=eval_cheb_even(f, fksalphax[:,i]);
    end

    Rf=MKinv*((1 ./ alpha).*fksalphax[:,end]);

    df = derCheb_even2odd(f); # f'

    global df_fks0
    df_fks0=eval_cheb_odd(df, fks0);

    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1));
    for i=1:m+1
    df_fksalphax[:,i]=eval_cheb_odd(df, fksalphax[:,i]);
    end

    ##################################
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1));


    for i=1:m+1
        evalf0 = cos.(myacos.(fks0[i].*ones(K+1,1))*transpose(0:2:2*K));
        evalf0[:,2:end] = 2*evalf0[:,2:end];
        h_fks0[:,:,i]=evalf0;
    end

    ########
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1));
    for i=1:m+1
        MKflambda = cos.(myacos.(fksalphax[:,i])*transpose(0:2:2*K));
        MKflambda[:,2:end] = 2*MKflambda[:,2:end];
        h_fksalphax[:,:,i]=MKflambda;
    end

    DR=MKinv*eval(Meta.parse(jacF));


    return Rf, DR

end
##############################################################

function Rm_DRm_even_1Dv1(f,m,x_K,MKinv,jacF_alpha,alpha)


    K = length(f)-1;


    global fks0
    fks0=BigFloat.(zeros(1,m+1)) .+ BigFloat.(zeros(1,m+1)).*im;

    for i=1:m
    
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alphaf=fks0[end];

    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    fksalphax[:,1]=alpha.*x_K;

    for i=1:m
        
        fksalphax[:,i+1]=eval_cheb_even(f, fksalphax[:,i]);
    end

    Rf=MKinv*((1 ./ alpha).*fksalphax[:,end]);

    df = derCheb_even2odd(f); # f'

    global df_fks0
    df_fks0=eval_cheb_odd(df, fks0);

    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    for i=1:m+1
    df_fksalphax[:,i]=eval_cheb_odd(df, fksalphax[:,i]);
    end

    ##################################
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;


    for i=1:m+1
        
        evalf0= BigFloat.(zeros(K+1,K+1)) .+ BigFloat.(zeros(K+1,K+1)).*im
        evalf0[:,1]=ones(K+1);
        evalf0[:,2]=2*fks0[i]^2 .*ones(K+1,1) .- 1;

        for k=3:K+1
            evalf0[:,k]=2*(evalf0[:,2]).*evalf0[:,k-1] .- evalf0[:,k-2]
        end

        evalf0[:,2:end] = 2*evalf0[:,2:end];
        h_fks0[:,:,i]=evalf0;
    end

    ########
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;
    for i=1:m+1
        

        MKflambda = BigFloat.(zeros(K+1,K+1)) .+ BigFloat.(zeros(K+1,K+1)).*im
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end

    ########################################################

    thisDR0f0=num_DRf0_even_1D0(m);

    thisDR0f= num_DRf_even_1D0(m);

    ~, forterm1 = der_f_term(m,0);

    thisterm1=forterm1 .* df_fksalphax[:,1];

    DR= MKinv*num_DRf_even_1Dv1(m, thisDR0f0, thisDR0f, thisterm1);
    ###########################################


    return alphaf, Rf, DR

end

#############################################

function Rm_DRm_even_1Dv1i(f,m,x_K,MKinv,jacF_alpha,alpha)

    K = length(f)-1;

    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;
    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;

    if isa(f[1],Complex{Interval{BigFloat}})
         fks0=Interval.(fks0); 
         fksalphax=Interval.(fksalphax);
         df_fksalphax=Interval.(df_fksalphax);
         h_fks0=Interval.(h_fks0);
         h_fksalphax=Interval.(h_fksalphax);
    end
    
    thismat= copy(h_fksalphax[:,:,1]);

    for i=1:m        
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alphaf=fks0[end];

    fksalphax[:,1]=alpha.*x_K;

    for i=1:m
        
        fksalphax[:,i+1]=eval_cheb_even(f, fksalphax[:,i]);
    end

    Rf=MKinv*((1 ./ alpha).*fksalphax[:,end]);
   
    df = derCheb_even2odd(f); # f'

    global df_fks0
    df_fks0=eval_cheb_odd(df, fks0);

    for i=1:m+1
    df_fksalphax[:,i]=eval_cheb_odd(df, fksalphax[:,i]);
    end

    ##################################
    

    for i=1:m+1
        

        evalf0= thismat;
        evalf0[:,1]=ones(K+1);
        evalf0[:,2]=2*fks0[i]^2 .*ones(K+1,1) .- 1;

        for k=3:K+1
            evalf0[:,k]=2*(evalf0[:,2]).*evalf0[:,k-1] .- evalf0[:,k-2]
        end

        evalf0[:,2:end] = 2*evalf0[:,2:end];

        h_fks0[:,:,i]=evalf0;
    end

    ########
    

    for i=1:m+1
        

        MKflambda = thismat;
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end
    
    ########################################################

    thisDR0f0=num_DRf0_even_1D0(m);

    thisDR0f= num_DRf_even_1D0(m);
    
    ~, forterm1 = der_f_term(m,0);
    
    thisterm1=forterm1 .* df_fksalphax[:,1];

    DR= MKinv*num_DRf_even_1Dv1(m, thisDR0f0, thisDR0f, thisterm1);

    ###########################################


    return alphaf, Rf, DR

end

#####################################################################################################


function tRm_DRm_even_1Dv1i(f,m,x_K,MKinv,jacF_alpha,alpha)

    K = length(f)-1;

    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;
    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;

    if isa(f[1],Complex{Interval{BigFloat}})
         fks0=Interval.(fks0); 
         fksalphax=Interval.(fksalphax);
         df_fksalphax=Interval.(df_fksalphax);
         h_fks0=Interval.(h_fks0);
         h_fksalphax=Interval.(h_fksalphax);
    end
    
    thismat= copy(h_fksalphax[:,:,1]);

    for i=1:m        
        
        thisval=tf(fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alphaf=fks0[end];

    fksalphax[:,1]=tl.*x_K;

    for i=1:m
        
        fksalphax[:,i+1]=tf(fksalphax[:,i]);
    end

    Rf=MKinv*((1 ./ tl).*fksalphax[:,end]);
    
    df = derCheb_even2odd(f); # f'

    global df_fks0
    df_fks0=dtf(fks0);

    for i=1:m+1
    
    df_fksalphax[:,i]=dtf(fksalphax[:,i]);
    end


    for i=1:m+1
        
        evalf0= thismat;
        evalf0[:,1]=ones(K+1);
        evalf0[:,2]=2*fks0[i]^2 .*ones(K+1,1) .- 1;

        for k=3:K+1
            evalf0[:,k]=2*(evalf0[:,2]).*evalf0[:,k-1] .- evalf0[:,k-2]
        end

        evalf0[:,2:end] = 2*evalf0[:,2:end];

        h_fks0[:,:,i]=evalf0;
    end


    for i=1:m+1
        
        
        MKflambda = thismat;
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end

    DR=MKinv*eval(Meta.parse(jacF_alpha));


    return alphaf, Rf, DR

end

###########################################################


function Rm_even_1D(f,m,x_K,MKinv)


    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;
    
    for i=1:m
        
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end
    
    global alpha0
    alpha0=fks0[end];
    
    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    fksalphax[:,1]=alpha0.*x_K;
    
    for i=1:m
       
        fksalphax[:,i+1]=eval_cheb_even(f, fksalphax[:,i]);
    end
    
    Rf=MKinv*(fksalphax[:,end]);
    
    
    return alpha0, Rf
    
end
    


#################################################


function get_miu(m)
    if m==2
    thismiu=1.392 #m=2
    elseif m==3
    thismiu=1.799 #m=3
    elseif m==4
    thismiu=1.944 #m=4
    elseif m==5
    #thismiu=0.536 #m=5
    thismiu=1.632
    elseif m==6
    thismiu=1.484 #m=6
    elseif m==7
    #thismiu=1.573#m=7
    thismiu=1.675
    elseif m==8
    thismiu=1.523 #m=8
    elseif m==9
    thismiu=1.596; #m=9
    elseif m==10
    thismiu=1.448; #m=10
    elseif m==11
    thismiu=1.56368; #m=11
    end
    return thismiu
end


get_miu_array=[0, [1.392], [1.799], [1.944], [1.632, 1.862, 1.9855], [1.484, 1.782, 1.9075], 
[1.675, 1.577, 1.8327, 1.8851, 1.9272], [1.523, 1.7114, 1.9424], [1.596, 1.6563], [1.448, 1.5019, 1.6310]
]
#####################################################################

function get_JacF(m)
    if m==2
    JacF="h_fksalphax[:,:,2].*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2]).*1/alpha+(h_fks0[:,:,2].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2]).*1/alpha-(h_fks0[:,:,2].*fksalphax[:,3]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*fksalphax[:,3]).*1/alpha^2"
    elseif m==3
    JacF="h_fksalphax[:,:,3].*1/alpha+(h_fksalphax[:,:,2].*df_fksalphax[:,3]).*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2].*df_fksalphax[:,3]).*1/alpha+(h_fks0[:,:,3].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3]).*1/alpha+(h_fks0[:,:,2].*df_fks0[3].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3]).*1/alpha-(h_fks0[:,:,3].*fksalphax[:,4]).*1/alpha^2-(h_fks0[:,:,2].*df_fks0[3].*fksalphax[:,4]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*fksalphax[:,4]).*1/alpha^2"
    elseif m==4
    JacF="h_fksalphax[:,:,4].*1/alpha+(h_fksalphax[:,:,3].*df_fksalphax[:,4]).*1/alpha+(h_fksalphax[:,:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha+(h_fks0[:,:,4].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha+(h_fks0[:,:,3].*df_fks0[4].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha+(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4]).*1/alpha-(h_fks0[:,:,4].*fksalphax[:,5]).*1/alpha^2-(h_fks0[:,:,3].*df_fks0[4].*fksalphax[:,5]).*1/alpha^2-(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*fksalphax[:,5]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*fksalphax[:,5]).*1/alpha^2"
    elseif m==5
    JacF="h_fksalphax[:,:,5].*1/alpha+(h_fksalphax[:,:,4].*df_fksalphax[:,5]).*1/alpha+(h_fksalphax[:,:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fksalphax[:,:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fks0[:,:,5].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fks0[:,:,4].*df_fks0[5].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5]).*1/alpha-(h_fks0[:,:,5].*fksalphax[:,6]).*1/alpha^2-(h_fks0[:,:,4].*df_fks0[5].*fksalphax[:,6]).*1/alpha^2-(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*fksalphax[:,6]).*1/alpha^2-(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*fksalphax[:,6]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*fksalphax[:,6]).*1/alpha^2"
    elseif m==6
    JacF="h_fksalphax[:,:,6].*1/alpha+(h_fksalphax[:,:,5].*df_fksalphax[:,6]).*1/alpha+(h_fksalphax[:,:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fksalphax[:,:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fksalphax[:,:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,5].*df_fks0[6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,4].*df_fks0[5].*df_fks0[6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6]).*1/alpha-(h_fks0[:,:,6].*fksalphax[:,7]).*1/alpha^2-(h_fks0[:,:,5].*df_fks0[6].*fksalphax[:,7]).*1/alpha^2-(h_fks0[:,:,4].*df_fks0[5].*df_fks0[6].*fksalphax[:,7]).*1/alpha^2-(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*fksalphax[:,7]).*1/alpha^2-(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*fksalphax[:,7]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*fksalphax[:,7]).*1/alpha^2"
    elseif m==7
    JacF="h_fksalphax[:,:,7].*1/alpha+(h_fksalphax[:,:,6].*df_fksalphax[:,7]).*1/alpha+(h_fksalphax[:,:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fksalphax[:,:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fksalphax[:,:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fksalphax[:,:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fksalphax[:,:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,5].*df_fks0[6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha+(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*x_K.*df_fksalphax[:,1].*df_fksalphax[:,2].*df_fksalphax[:,3].*df_fksalphax[:,4].*df_fksalphax[:,5].*df_fksalphax[:,6].*df_fksalphax[:,7]).*1/alpha-(h_fks0[:,:,7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,5].*df_fks0[6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2-(h_fks0[:,:,1].*df_fks0[2].*df_fks0[3].*df_fks0[4].*df_fks0[5].*df_fks0[6].*df_fks0[7].*fksalphax[:,8]).*1/alpha^2"
    # elseif m==8
    # thismiu=BigFloat(1.523) #m=8
    end
    return JacF
end

############################################################################################################


function der_f_term(n,j)   
    # also used for general terms in Jacobian
    # n is the order of the term (where term first apeared), j is starting number of iterations alphax inside of h, m=n+j
    # n>=2

    
    that=1;
    
    for k=1:n-1
        

        that= df_fksalphax[:,j+1+k].*that
        
    end

    forterm1 = copy(that)
    that= h_fksalphax[:,:,j+1].*that


    return that, forterm1          

end


function num_DRf_even_1D0(m)

    thissum= h_fksalphax[:,:,m] # this is m-1 iteration , matrix goes up to m iterations, last iteration does not appear in the jacobian 
   
    for n=0:m-2
        
        
        k=m-n;
        
        thisterm, ~ = der_f_term(k,n)  
        

        thissum=thissum .+ thisterm;
       

    end


    return thissum

end

####################################


function der_f0_term(n,j)   
    # also used for general terms in Jacobian
    # n is the order of the term (where term first apeared), j is starting number of iterations alphax inside of h, m=n+j
    # n>=2

    that=1;
    
    for k=1:n-1
        

        that= df_fks0[j+1+k].*that
       
    end

    
    that= h_fks0[:,:,j+1].*that
    

    return that

end


function num_DRf0_even_1D0(m)

    thissum= h_fks0[:,:,m] # this is m-1 iteration , matrix goes up to m iterations, last iteration does not appear in the jacobian 
    

    for n=0:m-2
        
        
        k=m-n;
         
        thisterm= der_f0_term(k,n)  
        

        thissum=thissum .+ thisterm;
       
    end


    return thissum

end

############################################


function num_DRf_even_1Dv1(m, DR0f0, DR0f, term1)

    ans = -(1.0 / alpha^2) .* DR0f0 .* fksalphax[:,m+1] .+ (1.0 / alpha) .* DR0f .+ (1.0 / alpha) .* term1 .* DR0f0 .* x_K; 


    return ans

end

############################################################################

function Z1Kinf_term(n,j, Z1KinfUpsilon)   
    # also used for general terms in Jacobian
    # n is the order of the term (where term first apeared), j is starting number of iterations alphax inside of h, m=n+j
    # n>=2

    that=1;
    
    for k=1:n-1
        
        that= df_fksalphax[:,j+1+k].*that
        
    end

    that= Z1KinfUpsilon.*that


    return that      

end


function num_Z1Kinf_1D0(m, Z1KinfUpsilon)

    #thissum= h_fksalphax[:,:,m] # this is m-1 iteration , matrix goes up to m iterations, last iteration does not appear in the jacobian 
    
    thissum = Z1KinfUpsilon

    for n=0:m-2
        
        
        k=m-n;
        
        thisterm = Z1Kinf_term(k,n, Z1KinfUpsilon)

        thissum=thissum .+ thisterm;
       

    end


    return thissum

end

####################################


function Z1Kinf0_term(n,j,Z1KinfUpsilon)   
    # also used for general terms in Jacobian
    # n is the order of the term (where term first apeared), j is starting number of iterations alphax inside of h, m=n+j
    # n>=2

   
    that=1;
    
    for k=1:n-1
        

        that= df_fks0[j+1+k].*that
        
    end

    
    that= Z1KinfUpsilon.*that
    

    return that

end


function num_Z1Kinf0_1D0(m,Z1KinfUpsilon)

    thissum= Z1KinfUpsilon

    for n=0:m-2
        
        
        k=m-n;
        
        thisterm= Z1Kinf0_term(k,n,Z1KinfUpsilon); 
        

        thissum=thissum .+ thisterm;
       

    end


    return thissum

end


#############################################################################
#############################################################################
############################### general purpose functs



function eye(n)
    #Id= 1*Matrix(I, n, n)
    Id= diagm(ones(n));
    return Id
end

function string_as_varname_function(s::AbstractString, v::Any)
s = Symbol(s)
@eval (($s) = ($v))
end

function heaviside(t)
   0.5 * (sign(t) + 1)
end

function doubleprimesum(i)
if (i==0.0 || i==1.0)
        return 0.5
else
        return 1
end
end


function my_ceil(x)
   
    y = Interval(ceil(inf(x)),ceil(sup(x)));

    return y
end


########################################
#############################################################
########################################### from part 2

function myacos(x)
    

    if typeof(x)==Sym
        ans=acos(x)
        return ans


    else


        try
            ans=acos.(x)

            return ans

        catch
            
            try
            size(x)[1]

            newx=copy(x)            

            thisin=0 .∈ imag(newx);
            thispt=newx[thisin];

            thisim= imag(thispt);

            
            this1=interval.(inf.(thisim),zeros(length(thispt),1)) .- eps(0.0);
            this2=interval.(zeros(length(thispt),1), sup.(thisim)) .+ eps(0.0);
            

            thispt1=real(thispt) + this1*im;
            thispt2=real(thispt) + this2*im;

            that1=acos.(thispt1);
            that2=acos.(thispt2);
           
            that= that1 .∪ that2;

            newx[thisin].=0.0;
            ans=acos.(newx);
            ans[thisin]=that;

            return ans            

            catch


            this1=inf(imag(x))..0.0;
            this2=0.0..sup(imag(x)) .+ eps(0.0);

            thispt1=real(x) + this1*im;
            thispt2=real(x) + this2*im;

            that1=acos.(thispt1);
            that2=acos.(thispt2);
            that= that1 .∪ that2;

            
            return that
            end
            
            

        end
       

    end
end

################################
##########################################################


function Rm_DRm_even_1D0(f,m,x_K,MKinv,jacF,alpha)


    K = length(f)-1;


    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;

    for i=1:m
        
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alpha0=fks0[end];

    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    fksalphax[:,1]=alpha.*x_K;

    for i=1:m
    
        fksalphax[:,i+1]=eval_cheb_even(f, fksalphax[:,i]);
    end

    Rf=MKinv*(fksalphax[:,end]);


    df = derCheb_even2odd(f); # f'

    global df_fks0
    df_fks0=eval_cheb_odd(df, fks0);

    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    for i=1:m+1
    df_fksalphax[:,i]=eval_cheb_odd(df, fksalphax[:,i]);
    end



    ########
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;

    for i=1:m+1
    
        MKflambda = BigFloat.(zeros(K+1,K+1)) .+ BigFloat.(zeros(K+1,K+1)).*im
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end

    DR=MKinv*eval(Meta.parse(jacF));


    DRalpha=df_fksalphax[:,1].*x_K;

    for i=2:m
        DRalpha=df_fksalphax[:,i].*DRalpha;
    end

    fx = eval_cheb_even(f, x_K);

    DRalpha=MKinv*(fx - DRalpha); 

    return alpha0, Rf, DR, DRalpha

end

##########################################################
####################################################################
function Rm_DRm_even_1D0i(f,m,x_K,MKinv,jacF,alpha)


    K = length(f)-1;

    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;
    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;



    if isa(f[1],Complex{Interval{BigFloat}})
         fks0=Interval.(fks0); 
         fksalphax=Interval.(fksalphax);
         df_fksalphax=Interval.(df_fksalphax);
         h_fks0=Interval.(h_fks0);
         h_fksalphax=Interval.(h_fksalphax);
    end

    thismat= copy(h_fksalphax[:,:,1]);

    for i=1:m
        
        thisval=eval_cheb_even(f,fks0[i]);
        fks0[i+1]=thisval[1];
    end

    global alpha
    alpha0=fks0[end];

    this=((alpha))*(x_K);
    fksalphax[:,1]=this;

    for i=1:m
        
        this=(eval_cheb_even((f), this));
       
        fksalphax[:,i+1]=this;
    end

    Rf=(MKinv)*(fksalphax[:,end]);


    df = derCheb_even2odd(f); # f'

    for i=1:m+1
    df_fksalphax[:,i]=eval_cheb_odd(df, fksalphax[:,i]);
    end


    for i=1:m+1
       
        MKflambda = thismat;
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end

    DR=MKinv*eval(Meta.parse(jacF));


    DRalpha=df_fksalphax[:,1].*x_K;

    for i=2:m
        DRalpha=df_fksalphax[:,i].*DRalpha;
    end

    fx = eval_cheb_even(f, x_K);

    DRalpha=MKinv*(fx - DRalpha); 

    return alpha0, Rf, DR, DRalpha

end

#####################################################


function est_PiKg_cste(gg,cste,K,rho,alpha_max,sub,nb)

    # The output coeffs contains upper bounds on (the absolute value of) the 
    # K+1 first Chebyshev coefficients of any analytic map f such that, for all
    # alpha >= rho, sup_{E_\alpha} |f| <= cste(rho,alpha) * sup_{E_\alpha} |g|.
    
   
    nb0=10;
    

    tab_coeffs = zeros(K+1,nb0);
    x=0:0.001:1;
    z=exp.((0+2im)*pi*x);
    

    ind_k = (0:2:2*K);
    
    
    tab_alpha = range(mid(rho),alpha_max,nb0);
    
    #ind = 1;
    g = gg[1];
    for ind=1:nb0
        println(ind)
        thisalpha = tab_alpha[ind]
        coeffs_decay = [1+2/(thisalpha^(2*ind_k[end])-1);
                        1 ./ thisalpha.^ind_k[2:end-1].*(thisalpha^(2*ind_k[end]) .+ thisalpha.^(2*ind_k[2:end-1]))/(thisalpha^(2*ind_k[end])-1);
                        1/thisalpha^ind_k[end]*thisalpha^(2*ind_k[end])/(thisalpha^(2*ind_k[end])-1)];

        Ealpha = (thisalpha*z+(thisalpha*z).^(-1))/2;
        
        this=abs.(g(Ealpha));

        g_C0_alpha =  norm(this,Inf);
       
        that=cste(mid(rho),thisalpha);
        try
        this= real.(that);
    
        tab_coeffs[:,ind] = g_C0_alpha*this[1]*coeffs_decay;
        catch
            
            continue

        end
        
    end


    thisI=Int.(zeros(size(tab_coeffs)));
    
    for k=1:size(tab_coeffs)[1]
    
        thisI[k,:]=sortperm(tab_coeffs[k,:]);

    end


    coeffs_test=tab_coeffs[:,thisI[:,1]];
    coeffs_test=diag(coeffs_test);


    tab_alpha = tab_alpha[thisI[:,1]];
    
    
    ################################## Rigorous computation 
    
    tab_g_C0_alpha = interval.(BigFloat.(zeros(K+1,1)));


    thisI = sortperm(tab_alpha);

    tab_alpha = interval.(tab_alpha);

    alpha_old = -Inf;
    g_C0_alpha_old = Inf;


    for k = 1:K+1
        thisalpha = tab_alpha[k];
        if thisalpha != alpha_old
            g_C0_alpha = compute_sup_E(gg,thisalpha,sub,nb);
        else
            g_C0_alpha = g_C0_alpha_old;
        end
        tab_g_C0_alpha[thisI[k]] = g_C0_alpha;
        alpha_old = thisalpha;
        g_C0_alpha_old = g_C0_alpha;
    end
    
    coeffs_decay = [1+2/(tab_alpha[1]^(2*ind_k[end])-1);
                    1 ./ tab_alpha[2:end-1].^ind_k[2:end-1].*(tab_alpha[2:end-1].^(2*ind_k[end])+tab_alpha[2:end-1].^(2*ind_k[2:end-1]))./(tab_alpha[2:end-1].^(2*ind_k[end]).-1);
                    1/tab_alpha[end]^ind_k[end]*tab_alpha[end]^(2*ind_k[end])/(tab_alpha[end]^(2*ind_k[end])-1)];

    coeffs = tab_g_C0_alpha.*cste(rho,tab_alpha).*coeffs_decay; 
    
    return coeffs

end

###########################################################################


function est_PiKg_cste1(gg,cste,K,rho,alpha_max,sub,nb)

    # The output coeffs contains upper bounds on (the absolute value of) the 
    # K+1 first Chebyshev coefficients of any analytic map f such that, for all
    # alpha >= rho, sup_{E_\alpha} |f| <= cste(rho,alpha) * sup_{E_\alpha} |g|.
    
   
    nb0=10;
    

    tab_coeffs = zeros(K+1,nb0);
    x=0:0.001:1;
    z=exp.((0+2im)*pi*x);
    

    ind_k = (0:2:2*K);
    
    
    tab_alpha = range(mid(rho),alpha_max,nb0);
    
    #ind = 1;
    g = gg[1];
    for ind=1:nb0
        println(ind)
        thisalpha = tab_alpha[ind]
        coeffs_decay = [1+2/(thisalpha^(2*ind_k[end])-1);
                        1 ./ thisalpha.^ind_k[2:end-1].*(thisalpha^(2*ind_k[end]) .+ thisalpha.^(2*ind_k[2:end-1]))/(thisalpha^(2*ind_k[end])-1);
                        1/thisalpha^ind_k[end]*thisalpha^(2*ind_k[end])/(thisalpha^(2*ind_k[end])-1)];

        #global Ealpha
        Ealpha = (thisalpha*z+(thisalpha*z).^(-1))/2;
        
        this=abs.(g(Ealpha));

        g_C0_alpha =  norm(this,Inf);
       
        that=cste(mid(rho),thisalpha);
        try
        this= real.(that);
    
        tab_coeffs[:,ind] = g_C0_alpha*this[1]*coeffs_decay;
        catch
            
            continue

        end
        
    end


    thisI=Int.(zeros(size(tab_coeffs)));
    
    for k=1:size(tab_coeffs)[1]
    
        thisI[k,:]=sortperm(tab_coeffs[k,:]);

    end


    coeffs_test=tab_coeffs[:,thisI[:,1]];
    coeffs_test=diag(coeffs_test);


    tab_alpha = tab_alpha[thisI[:,1]];
    
    
    ################################## Rigorous computation 
    
    tab_g_C0_alpha = interval.(BigFloat.(zeros(K+1,1)));


    thisI = sortperm(tab_alpha);

    tab_alpha = interval.(tab_alpha);

    alpha_old = -Inf;
    g_C0_alpha_old = Inf;


    for k = 1:K+1
        thisalpha = tab_alpha[k];
        if thisalpha != alpha_old
            g_C0_alpha = compute_sup_E(gg,thisalpha,sub,nb);
        else
            g_C0_alpha = g_C0_alpha_old;
        end
        tab_g_C0_alpha[thisI[k]] = g_C0_alpha;
        alpha_old = thisalpha;
        g_C0_alpha_old = g_C0_alpha;
    end
    
    coeffs_decay = [1+2/(tab_alpha[1]^(2*ind_k[end])-1);
                    1 ./ tab_alpha[2:end-1].^ind_k[2:end-1].*(tab_alpha[2:end-1].^(2*ind_k[end])+tab_alpha[2:end-1].^(2*ind_k[2:end-1]))./(tab_alpha[2:end-1].^(2*ind_k[end]).-1);
                    1/tab_alpha[end]^ind_k[end]*tab_alpha[end]^(2*ind_k[end])/(tab_alpha[end]^(2*ind_k[end])-1)];

    coeffs = tab_g_C0_alpha.*cste(rho,tab_alpha).*coeffs_decay; 
    
    return coeffs, tab_g_C0_alpha

end


#############################################################################

function est_PiKg_cste_adap(gg,cste,K,rho,alpha_max,subtol)

    # The output coeffs contains upper bounds on (the absolute value of) the 
    # K+1 first Chebyshev coefficients of any analytic map f such that, for all
    # alpha >= rho, sup_{E_\alpha} |f| <= cste(rho,alpha) * sup_{E_\alpha} |g|.
    
   
    nb0=10;
    

    tab_coeffs = zeros(K+1,nb0);
    x=0:0.001:1;
    z=exp.((0+2im)*pi*x);
    

    ind_k = (0:2:2*K);
    
    
    tab_alpha = range(mid(rho),alpha_max,nb0);
    
    #ind = 1;
    g = gg[1];
    for ind=1:nb0
        println(ind)
        thisalpha = tab_alpha[ind]
        coeffs_decay = [1+2/(thisalpha^(2*ind_k[end])-1);
                        1 ./ thisalpha.^ind_k[2:end-1].*(thisalpha^(2*ind_k[end]) .+ thisalpha.^(2*ind_k[2:end-1]))/(thisalpha^(2*ind_k[end])-1);
                        1/thisalpha^ind_k[end]*thisalpha^(2*ind_k[end])/(thisalpha^(2*ind_k[end])-1)];

        #global Ealpha
        Ealpha = (thisalpha*z+(thisalpha*z).^(-1))/2;
        
        this=abs.(g(Ealpha));

        g_C0_alpha =  norm(this,Inf);
       
        that=cste(mid(rho),thisalpha);
        try
        this= real.(that);
    
        tab_coeffs[:,ind] = g_C0_alpha*this[1]*coeffs_decay;
        catch
            
            continue

        end
        
    end


    thisI=Int.(zeros(size(tab_coeffs)));
    
    for k=1:size(tab_coeffs)[1]
    
        thisI[k,:]=sortperm(tab_coeffs[k,:]);

    end


    coeffs_test=tab_coeffs[:,thisI[:,1]];
    coeffs_test=diag(coeffs_test);


    tab_alpha = tab_alpha[thisI[:,1]];
    
    
    ################################## Rigorous computation 
    
    tab_g_C0_alpha = interval.(BigFloat.(zeros(K+1,1)));


    thisI = sortperm(tab_alpha);

    tab_alpha = interval.(tab_alpha);

    alpha_old = -Inf;
    g_C0_alpha_old = Inf;


    for k = 1:K+1
        thisalpha = tab_alpha[k];
        if thisalpha != alpha_old
            #g_C0_alpha = compute_sup_E(gg,thisalpha,sub,nb);
            g_C0_alpha = compute_sup_E1_adap(gg,thisalpha,subtol[k])

        else
            g_C0_alpha = g_C0_alpha_old;
        end
        tab_g_C0_alpha[thisI[k]] = g_C0_alpha;
        alpha_old = thisalpha;
        g_C0_alpha_old = g_C0_alpha;
    end
    
    coeffs_decay = [1+2/(tab_alpha[1]^(2*ind_k[end])-1);
                    1 ./ tab_alpha[2:end-1].^ind_k[2:end-1].*(tab_alpha[2:end-1].^(2*ind_k[end])+tab_alpha[2:end-1].^(2*ind_k[2:end-1]))./(tab_alpha[2:end-1].^(2*ind_k[end]).-1);
                    1/tab_alpha[end]^ind_k[end]*tab_alpha[end]^(2*ind_k[end])/(tab_alpha[end]^(2*ind_k[end])-1)];

    coeffs = tab_g_C0_alpha.*cste(rho,tab_alpha).*coeffs_decay; 
    
    return coeffs

end


#################################################################################
function est_PiKg_cste_fft(gg,gg_coeffs,cste,K,rho,alpha_max,sub,nb)

    # The output coeffs contains upper bounds on (the absolute value of) the 
    # K+1 first Chebyshev coefficients of any analytic map f such that, for all
    # alpha >= rho, sup_{E_\alpha} |f| <= cste(rho,alpha) * sup_{E_\alpha} |g|.
    
    
    nb0=10;
   

    tab_coeffs = zeros(K+1,nb0);
    x=0:0.001:1;
    z=exp.((0+2im)*pi*x);
    

    ind_k = (0:2:2*K);
        
    tab_alpha = range(mid(rho),alpha_max,nb0);
    
    #ind = 1;
    g = gg[1];
    for ind=1:nb0
        println(ind)
        thisalpha = tab_alpha[ind]
        coeffs_decay = [1+2/(thisalpha^(2*ind_k[end])-1);
                        1 ./ thisalpha.^ind_k[2:end-1].*(thisalpha^(2*ind_k[end]) .+ thisalpha.^(2*ind_k[2:end-1]))/(thisalpha^(2*ind_k[end])-1);
                        1/thisalpha^ind_k[end]*thisalpha^(2*ind_k[end])/(thisalpha^(2*ind_k[end])-1)];

        #global Ealpha
        Ealpha = (thisalpha*z+(thisalpha*z).^(-1))/2;
        
        this=abs.(g(Ealpha));

        g_C0_alpha =  norm(this,Inf);
        
        
        that=cste(mid(rho),thisalpha);
        try
        this= real.(that);
    
        tab_coeffs[:,ind] = g_C0_alpha*this[1]*coeffs_decay;
        catch
            
            continue

        end
        
    end


    thisI=Int.(zeros(size(tab_coeffs)));
    
    for k=1:size(tab_coeffs)[1]
    
        thisI[k,:]=sortperm(tab_coeffs[k,:]);

    end

    coeffs_test=tab_coeffs[:,thisI[:,1]];
    coeffs_test=diag(coeffs_test);


    tab_alpha = tab_alpha[thisI[:,1]];
    
    
    ################################## Rigorous computation 
    
    tab_g_C0_alpha = interval.(BigFloat.(zeros(K+1,1)));


    thisI = sortperm(tab_alpha);

    tab_alpha = interval.(tab_alpha);

    alpha_old = -Inf;
    g_C0_alpha_old = Inf;

    for k = 1:K+1
        thisalpha = tab_alpha[k];
        if thisalpha != alpha_old
            
            g_C0_alpha = norm_C0_Enu1(gg_coeffs,thisalpha);
            

        else
            g_C0_alpha = g_C0_alpha_old;
        end
        tab_g_C0_alpha[thisI[k]] = g_C0_alpha;
        alpha_old = thisalpha;
        g_C0_alpha_old = g_C0_alpha;
    end
    

    coeffs_decay = [1+2/(tab_alpha[1]^(2*ind_k[end])-1);
    1 ./ tab_alpha[2:end-1].^ind_k[2:end-1].*(interval(1.0) .+ tab_alpha[2:end-1].^(-2*(ind_k[end].-ind_k[2:end-1])))./(interval(1.0) .- tab_alpha[2:end-1].^(-2*ind_k[end]));
    1/tab_alpha[end]^ind_k[end]*(interval(1.0)/(interval(1.0)-tab_alpha[end]^(-2*ind_k[end])))];

    coeffs = tab_g_C0_alpha.*cste(rho,tab_alpha).*coeffs_decay; 
    
    return coeffs

end


#####################################################################


function eval_cheb_even_num(f,x)
    

    thisf=copy(f);
    thisx=copy(x);
    # Computation of f_0 + 2*sum_{k=1}^K f_{2k} T_{2k}(x), where f is a vector
    # containing the coefficients f_{2k}. If x is a vector, this formula is
    # applied component-wise to x.


    K = length(thisf)-1;
    if size(thisx,1) == 1
        thisx = transpose(thisx);
    end


    Mat = cos.(acos.(thisx)*transpose(0:2:2*K));
    thisf[2:end] = 2*thisf[2:end];
    y = Mat * thisf;
    thisf= nothing;
    thisx= nothing;
    return y

end

#######################################################


function compute_sup_E(gg,thisalpha,sub,nb)

        

        #(ff,alpha,int,nb)
    # Rigorously encloses sup_{E_\alpha} |f|. By the maximum modulus principle, 
    # we only need to compute the supremum of |f| on the boundary of E_\alpha.
    #
    # ff must be a cell array, with ff{1} = f. Its is possible to also provide
    # derivatives of f (i.e f' in ff{2}, f'' in ff{3}, etc). If derivates are
    # given, they are used to provide a more accurate answer.
    #
    # The boundary of of E_\alpha is parametrized as 
    # 1/2 * ( alpha*e^{i*theta} + 1/(alpha*e^{i*theta}) ), for theta in [0,2pi]
    # A different range for theta can be inputed (symmetry considerations may
    # be used to reduce this range), via the input "int"
    #
    # nb describes the number of subintervals into which the range of values of
    # theta is split for the interval arithmetic computations.



    if isa(thisalpha[1],Interval{})

        theta_p = interval.(sub);
        
        dtheta = norm(sup.(theta_p[2:end]-theta_p[1:end-1]),Inf)/2;


        theta_i = sub .± dtheta;

        dz = 1/2 * dtheta * (thisalpha+1/thisalpha)/2;

        N = length(gg);
    
        g = interval.(BigFloat.(zeros(nb,1)));

    else

        theta_p = sub;
        
        dtheta = norm(theta_p[2:end]-theta_p[1:end-1],Inf)/2;
        
        
        theta_i = sub;
        
        dz = 1/2 * dtheta * (thisalpha+1/thisalpha)/2;
        
        N = length(gg);
        
        g = BigFloat.(zeros(nb,1));


    end


    for n = 1:N
        
        dng(theta) = abs.( gg[n]( 1/2*(thisalpha*exp.(1im*theta) + 1/thisalpha*exp.(-1im*theta)) ) );

        if n<N 
            arg = theta_p;
        else
            arg = theta_i;
        end
    
        g = g .+ (dng(arg)).*(dz^(n-1)/factorial(n-1));


    end

    val = norm(g,Inf);

    return val

end


################################################################################



function compute_sup_E1(gg,thisalpha,sub)

    nb= length(sub)

        #(ff,alpha,int,nb)
    # Rigorously encloses sup_{E_\alpha} |f|. By the maximum modulus principle, 
    # we only need to compute the supremum of |f| on the boundary of E_\alpha.
    #
    # ff must be a cell array, with ff{1} = f. Its is possible to also provide
    # derivatives of f (i.e f' in ff{2}, f'' in ff{3}, etc). If derivates are
    # given, they are used to provide a more accurate answer.
    #
    # The boundary of of E_\alpha is parametrized as 
    # 1/2 * ( alpha*e^{i*theta} + 1/(alpha*e^{i*theta}) ), for theta in [0,2pi]
    # A different range for theta can be inputed (symmetry considerations may
    # be used to reduce this range), via the input "int"
    #
    # nb describes the number of subintervals into which the range of values of
    # theta is split for the interval arithmetic computations.


    theta_i = copy(sub);
    theta_p = mid.(sub);

    #dtheta = norm(sup.(theta_p[2:end]-theta_p[1:end-1]),Inf)/2;
    #dtheta= sup(theta_p[end])/(2*nb);

    dtheta =(sup(theta_i[1]) - inf(theta_i[1]))/2


    dz = 1/2 * dtheta * (thisalpha+1/thisalpha)/2;

    N = length(gg);

    g = interval.(BigFloat.(zeros(nb,1)));


    for n = 1:N
        
        dng(theta) = abs.( gg[n]( 1/2*(thisalpha*exp.(1im*theta) + 1/thisalpha*exp.(-1im*theta)) ) );

        if n<N 
            arg = theta_p;
        else
            arg = theta_i;
        end

        g = g .+ (dng(arg)).*(dz^(n-1)/factorial(n-1));


    end

    #val = norm(g,Inf);

    return g

end

###################################

function compute_sup_E1_adap(gg, thisalpha,subtol)

        nb=1000;
        sub = BigFloat.(range(0,2*pi,nb)); 

        nall= copy(nb);
        
        theta_p = interval.(sub);
        
        #dtheta = norm(sup.(theta_p[2:end]-theta_p[1:end-1]),Inf)/2;
        dtheta= pi/(nb);    

        theta_i = sub .± dtheta;

        temp= compute_sup_E1(gg,thisalpha,theta_i);
        
        temp1= sup.(temp) ./ (inf.(temp).+ 1.0);

        temp_good= temp1 .< subtol;

        temp_bad= .!(temp_good[:,1]);
        
        thatmax= norm(temp[temp_good[:,1]],Inf)        

        nbad= sum(temp_bad)  


        thissub= copy(theta_i);
        

        while  nbad > 0
                            
            nsubs= maximum([2, Int(ceil(nbad/nall*5))]);
            
            dtheta= dtheta/(nsubs)

            thatsub= inf.(thissub[temp_bad]);

            
            nall=nsubs*nbad;


            thissub= interval.(BigFloat.(zeros(nall,1)));

            thissub[1:(nsubs):nall] = interval.(thatsub,thatsub .+ dtheta);

            Threads.@threads for i=2:nsubs
            
            thissub[i:(nsubs):nall] =  interval.(thatsub .+ (i-1)*dtheta,thatsub .+ i*dtheta)
            end

            temp= compute_sup_E1(gg,thisalpha,thissub);
            
            temp1= sup.(temp) ./ (inf.(temp).+ 1.0);
            
            temp_good= temp1 .< subtol;
            
            temp_bad= .!(temp_good[:,1]);
                
            thismax= norm(temp[temp_good[:,1]],Inf)
             

            thatmax= norm([thismax; thatmax], Inf);
            

            nbad= sum(temp_bad);

            
        end


    return thatmax

end

################################################################

function Upsilon_old(rho,alpha,p,q,K)
    
        # The constants called \Upsilon^{p,q}_{rho,alpha} in the paper
        
        this= norm(inf.(real.(rho)),Inf);
        
    
        if this > sup(alpha)
           
            val = Inf;
            return
        end
        
        if q == 1 
            if p == 1
                val = (1+1 ./ rho.^2).*(rho./alpha).^(K+1);
            elseif p == 0
                val = 1/2*( (1 .+ 1 ./ rho.^2).*(rho./alpha).^(K+1) + (1 .+ rho.^2).*1 ./ (rho.*alpha).^(K+1) );
            else
                error("Invalid value for the input 'p'")
            end
        elseif q == 0
            if rho == alpha
                val = Inf;
                return
            end
            if p == 1
                val = 2 ./ (alpha.^(2*K)-1).*(rho./(alpha-rho).*(1-(rho./alpha).^K)+((alpha.*rho).^K-1)./(alpha.*rho-1)) + 2*(rho./alpha).^K.*rho./(alpha-rho);  
            elseif p == 0
                val = 1 ./ (alpha.^(2*K)-1) .* ( rho./(alpha-rho).*(1-(rho./alpha).^K) + 1 ./(alpha*rho-1).*(1-(1 ./ (alpha.*rho)).^K) + ((alpha.*rho).^K-1)./(alpha.*rho-1) + ((alpha./rho).^K-1)./(alpha./rho-1) ) + (rho ./ alpha).^K.*rho ./ (alpha-rho) + (1 ./ (alpha.*rho)).^K.*1 ./ (alpha.*rho-1);
            else
                error("Invalid value for the input 'p'")
            end
        else
            error("Invalid value for the input 'q'")
        end
    
        return val
    end

########################################################################

function Upsilon(rho,nu,p,q,K)
    
    # The constant called \Upsilon^{p,q,even}_{rho,nu,K} in the paper
        
    if mod(K,2) == 1 
        error("Only the even case is implemented")
    end

    rho_max = maximum(sup.(real.(rho))); 
    
    if p == 0 && q == 1
        if rho_max > inf(nu)
            val = Inf;
        else
            val = 1/2 * (rho/nu).^(K+2) .* (1 .+ rho.^(-4) .+ rho.^(-2*K) .+ rho.^(-2*(K+2)) );
        end
    
    elseif p == 1 && q == 0 
        if rho_max >= inf(nu)
            val = Inf;
        else
            val = 2/(nu^(2*K)-1) * ( rho.^2 ./(nu^2-rho.^2).*(1-(rho/nu).^K) .+ ((nu*rho).^K-1)./((nu*rho).^2-1) ) .+ 2*rho.^2 ./(nu^2-rho.^2).*(rho/nu).^K
        end
    
    else
        error("This combination of 'p' and 'q' is not implemented, you should not need it for this work")
    end

    return val
end

########################################################################


function opt_para(g, rho, cste, alpha_max, tol,sub)
     
    
        # Finds alpha which minimizes sup_{E_\alpha} |g| * cste (rho,alpha). Not
        # rigorous.
        
        
        z = exp.(2im*pi.*sub);
    
    
        alpha_min = rho + tol;
        Ealpha_min = (alpha_min*z+(alpha_min*z).^(-1))/2;
        Ealpha_max = (alpha_max*z+(alpha_max*z).^(-1))/2;
    
    
        val_min = norm(abs.(g(Ealpha_min)),Inf) * cste(rho,alpha_min);
        val_max = norm(abs.(g(Ealpha_max)),Inf) * cste(rho,alpha_max);
        
        
        while alpha_max-alpha_min>tol
            alpha_1 = (2*alpha_min+alpha_max)/3;
            Ealpha_1 = (alpha_1*z+(alpha_1*z).^(-1))/2;
            val_1 = norm(abs.(g(Ealpha_1)),Inf) * cste(rho,alpha_1);
            if val_min < val_1
                alpha_max = alpha_1;
            else
                alpha_2 = (alpha_min+2*alpha_max)/3;
                Ealpha_2 = (alpha_2*z+(alpha_2*z).^(-1))/2;
                val_2 = norm(abs.(g(Ealpha_2)),Inf) * cste(rho,alpha_2);
                if val_1 < val_2
                    alpha_max = alpha_2;
                elseif val_2 < val_max
                    alpha_min = alpha_1;
                else
                    alpha_min = alpha_2;
                end
            end
        end
    
        alpha = (alpha_min + alpha_max)/2;
        
    
        return alpha
 end

 ########################################################################################

 

function opt_para1(g, rho, cste, alpha_max, tol,sub)

     
    
    # Finds alpha which minimizes sup_{E_\alpha} |g| * cste (rho,alpha). Not
    # rigorous.
    
    
    z = exp.(2im*pi.*sub);


    alpha_min = rho + tol;
    Ealpha_min = (alpha_min*z+(alpha_min*z).^(-1))/2;
    Ealpha_max = (alpha_max*z+(alpha_max*z).^(-1))/2;


    val_min = norm(abs.(g(Ealpha_min)),Inf) * cste(rho,alpha_min);
    val_max = norm(abs.(g(Ealpha_max)),Inf) * cste(rho,alpha_max);
    
    
    while alpha_max-alpha_min>tol
        alpha_1 = (2*alpha_min+alpha_max)/3;
        Ealpha_1 = (alpha_1*z+(alpha_1*z).^(-1))/2;
        val_1 = norm(abs.(g(Ealpha_1)),Inf) * cste(rho,alpha_1);
        if val_min < val_1
            alpha_max = alpha_1;
            val_max = val_1;
        else
            alpha_2 = (alpha_min+2*alpha_max)/3;
            Ealpha_2 = (alpha_2*z+(alpha_2*z).^(-1))/2;
            val_2 = norm(abs.(g(Ealpha_2)),Inf) * cste(rho,alpha_2);
            if val_1 < val_2
                alpha_max = alpha_2;
                val_max = val_2
            elseif val_2 < val_max
                alpha_min = alpha_1;
                val_min = val_1
            else
                alpha_min = alpha_2;
                val_min = val_2
            end
        end
    end


    alpha = (alpha_min + alpha_max)/2;
    
    val = (val_min + val_max)/2;    

    return alpha, val
end


##################################################################################



function sigma(rho,alpha,p,q)

    # The constants called \sigma^{p,q}_{rho,alpha} in the paper
    
    if rho > alpha
        val = Inf;
        return
    end
    
    if q == 1 
        if p == 1
            error("Estimate not implemented")
        elseif p == 0
            kalpha = max(my_ceil(1/log(alpha))-1,1);
            if rho == 1
                val1 = exp(-1)/log(alpha)*(1+2*(kalpha-1));
            else
                val1 = exp(-1)/log(alpha)*(1+(rho^kalpha-1)/(rho-1)+(1-rho^(-kalpha))/(1-rho^(-1)));
            end
            val2 = (rho/alpha)^kalpha*((alpha-rho)*kalpha+alpha)/(alpha-rho)^2;
            val3 = rho/(rho*alpha)^kalpha*((alpha*rho-1)*kalpha+alpha*rho)/(alpha*rho-1)^2;
            val = 2*( val1 + val2 + val3 );
        else
            error("Invalid value for the input 'p'")
        end
    elseif q == 0
        if p == 1
            val = 2*( 2*alpha/(alpha^2-1)^2*(alpha+rho)/(alpha-rho) + alpha/(alpha^2-1)*((alpha+rho)/(alpha-rho) + 2*alpha*rho/(alpha-rho)^2) );
        elseif p == 0
            error("Estimate not implemented")    
        else
            error("Invalid value for the input 'p'")
        end
    else
        error("Invalid value for the input 'q'")
    end
    
    return val
end
    
#####################################################

function sigma_new(rho,v)
    a=(rho/v)^2;
    b= (1/v)^2; 
    


    n0= ceil(1/(2*(log(v)-log(rho))));
    n0= sup(n0);

    if rho==1
        thismax=1;
        thissigmas= 4*(1:n0).^2 .*b.^(1:n0);
        this=maximum(thissigmas);
    else 
        thismax= 2*n0/(v^(2*n0))*(rho^(2*n0)*rho/(rho^2-1)+rho^(-1)/(1-rho^(-2)));
        
        thissigmas= (1:n0).*(a.^(1:n0)-b.^(1:n0));
        thissigma=maximum(thissigmas);
        thissigma=2*rho/(rho^2-1)*(thissigma);
        this=maximum([sup(thissigma),sup(thismax)]);

    end

    return this
end 


################################################################################

function tRm_DRm_even_1D0i(f,m,x_K,MKinv,jacF,alpha)


    K = length(f)-1;

    global fks0
    fks0=BigFloat.(zeros(1,m+1)).+ BigFloat.(zeros(1,m+1)).*im;
    global fksalphax
    fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global df_fksalphax
    df_fksalphax=BigFloat.(zeros(K+1,m+1)) .+ BigFloat.(zeros(K+1,m+1)).*im;
    global h_fks0
    h_fks0=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;
    global h_fksalphax
    h_fksalphax=BigFloat.(zeros(K+1,K+1,m+1)) .+ BigFloat.(zeros(K+1,K+1,m+1)).*im;


    if isa(f[1],Complex{Interval{BigFloat}})
         fks0=Interval.(fks0); 
         fksalphax=Interval.(fksalphax);
         df_fksalphax=Interval.(df_fksalphax);
         h_fks0=Interval.(h_fks0);
         h_fksalphax=Interval.(h_fksalphax);
    end

    
    thismat= copy(h_fksalphax[:,:,1]);


    
    for i=1:m
        
        thisval=tf(fks0[i])
        fks0[i+1]=thisval[1];
    end

    global alpha
    alpha0=fks0[end];

    fksalphax[:,1]=tl.*x_K;

    for i=1:m
        
        fksalphax[:,i+1]=tf(fksalphax[:,i]);
    end

    Rf=MKinv*(fksalphax[:,end]);

    for i=1:m+1
   
    df_fksalphax[:,i]=dtf(fksalphax[:,i]);
    
    end


    for i=1:m+1
        
        MKflambda = thismat;
        MKflambda[:,1]=ones(K+1);
        MKflambda[:,2]=2*fksalphax[:,i].^2 .- 1;

        for k=3:K+1
            MKflambda[:,k]=2*(MKflambda[:,2]).*MKflambda[:,k-1] .- MKflambda[:,k-2]
        end

        MKflambda[:,2:end] = 2*MKflambda[:,2:end];

        h_fksalphax[:,:,i]=MKflambda;
    end

    DR=MKinv*eval(Meta.parse(jacF));


    DRalpha=df_fksalphax[:,1].*x_K;

    for i=2:m
        DRalpha=df_fksalphax[:,i].*DRalpha;
    end

    fx=tf(x_K);


    DRalpha=MKinv*(fx - DRalpha); 

    return alpha0, Rf, DR, DRalpha

end

#######################################################################
###################################

function inclusion_gE(g,thisalpha,sub,nb)

    # Finds the smallest beta such that g(E_alpha) is included in E_beta.
    
    thisN = length(thisalpha);
    eta = Interval.(BigFloat.(zeros(thisN,1)));
    
    foo(z)=abs.(g(z).- 1.0) .+ abs.(g(z) .+ 1.0);
    
    for j = 1:thisN
        
        eta[j] = compute_sup_E([foo], thisalpha[j],sub,nb);
    end
    
    
    beta = ( eta .+ sqrt.(eta.^2 .- 4) ) ./ 2;
   
    return beta
end

##########################################################################

function inclusion_gE2(g,thisalpha,subtol)

    # Finds the smallest beta such that g(E_alpha) is included in E_beta.
    
    thisN = length(thisalpha);
    eta = Interval.(BigFloat.(zeros(thisN,1)));
    
    foo(z)=abs.(g(z).- 1.0) .+ abs.(g(z) .+ 1.0);
    
    for j = 1:thisN
        
        #eta[j] = compute_sup_E([foo], thisalpha[j],sub,nb);
        eta[j] = compute_sup_E1_adap([foo], thisalpha[j],subtol)
    end
    
    
    beta = ( eta .+ sqrt.(eta.^2 .- 4) ) ./ 2;
   
    return beta
end


########################################################################
function inclusion_gE1(g,thisalpha,sub,nb)

    # Finds the smallest beta such that g(E_thisalpha) is included in E_beta.
        
    thisN = length(thisalpha);
    eta = Interval.(BigFloat.(zeros(thisN,1)));
    
    foo(z)=abs.(g(z).- 1.0) .+ abs.(g(z) .+ 1.0);

    foo_coeffs = MKinv*foo(x_K);
    
    for j = 1:thisN
        
        
        eta[j] = norm_C0_Enu1(foo_coeffs,thisalpha[j]);


    end
    
    
    beta = ( eta .+ sqrt.(eta.^2 .- 4) ) ./ 2;
   
    return beta
end

###########################################################################


function opt_incl(g, beta,alpha_max,tol,sub,nb)

    # Finds the largest alpha such that g(E_alpha) included in E_beta
    
    
    thisalpha = 1.0;
    gamma = inclusion_gE(g,thisalpha,sub,nb);

    this= norm(inf.(real.(gamma)),Inf);
    
    if this > beta
        
        error("Impossible to find an inclusion")
    else
        it = 0;
        while (alpha_max-thisalpha)>tol && it < 1e5
            alpha_old = thisalpha;
            thisalpha = (thisalpha+alpha_max)/2;
            gamma = inclusion_gE(g,thisalpha,sub,nb);
            
            this= norm(inf.(real.(gamma)),Inf);
    
            if this > beta
               
                alpha_max = thisalpha;
                thisalpha = alpha_old;
            end
            it = it+1;
        end
        if (alpha_max-thisalpha)>tol
            warning("Prescribed tolerance not reached")
        end
    end

    return thisalpha
end
                
##################################################################

function opt_incl1(g, beta,alpha_max,tol,sub,nb)

    # Finds the largest alpha such that g(E_alpha) included in E_beta
    # uses norm_C0_Enu1

    thisalpha = 1.0;
    gamma = inclusion_gE1(g,thisalpha,sub,nb);

    this= norm(inf.(real.(gamma)),Inf);
    
    if this > beta
       
        error("Impossible to find an inclusion")
    else
        it = 0;
        while (alpha_max-thisalpha)>tol && it < 1e5
            alpha_old = thisalpha;
            thisalpha = (thisalpha+alpha_max)/2;
            gamma = inclusion_gE1(g,thisalpha,sub,nb);
            
            this= norm(inf.(real.(gamma)),Inf);
    
            if this > beta
                #gamma > beta
                alpha_max = thisalpha;
                thisalpha = alpha_old;
            end
            it = it+1;
        end
        if (alpha_max-thisalpha)>tol
            warning("Prescribed tolerance not reached")
        end
    end

    return thisalpha
end
    
 ##############################################################################

function my_fft_plus(u)

    if length(u)==1
       return u
    else
        N=length(u)
        u_0= my_fft_plus(u[2:2:end]);
        u_e= my_fft_plus(u[1:2:end]);
        thisfac= exp.(2im * pi *(1:(N/2 ))/N );
        fftu= [u_e .+ thisfac .* u_0; u_e .- thisfac .* u_0];
    end

    return fftu

end

################################################################


function my_fft_plus_1(u,N0)

    if length(u)==1
       return u
    else
        N=length(u)
        u_0= my_fft_plus(u[2:2:end]);
        u_e= my_fft_plus(u[1:2:end]);
        thisfac= exp.(2im * pi *(N0:(N/2 ))/N );
        fftu= [u_e .+ thisfac .* u_0; u_e .- thisfac .* u_0];
    end

    return fftu

end

##################################################################

function my_fft_neg(u)

    if length(u)==1
       return u
    else
        N=length(u)
        u_0= my_fft_neg(u[2:2:end]);
        u_e= my_fft_neg(u[1:2:end]);
        thisfac= exp.(-2im * pi *(1:(N/2 ))/N );
        fftu= [u_e .+ thisfac .* u_0; u_e .- thisfac .* u_0];
    end

    return fftu

end

###############################################################

function norm_C0_Enu1(u,nu)

    # Inputs: a function g, a vector of Chebyshev coefficients u, a real nu>=1,
    # an integer NFFT (optional).
    # Assumptions: g is analytic on a domain containing the range of (the
    # function associated to) u over the Bernstein ellipse E_nu. 
    # Output: sup_{z in E_nu} |g(u(z))|. If intvals are used, the output is a 
    # guaranteed enclosure of the supremum. 
    
    N = length(u) - 1;

    

    NFFT = Int(2^ceil(log2(N+1)));
    

    wplus = Interval.(BigFloat.(zeros(NFFT,1)) + 1im*BigFloat.(zeros(NFFT,1)));
    
    wneg = Interval.(BigFloat.(zeros(NFFT,1)) + 1im*BigFloat.(zeros(NFFT,1)));


    delta = interval(-1,1) * sup(@interval pi/NFFT) + 1im*Interval(BigFloat(0.0));
    
    
    wplus[1:N] = u[2:end] .* (nu).^(1:N) .* exp.(1im*(1:N).*delta);
    wneg[1:N] = u[2:end] .* (nu) .^(-(1:N)) .* exp.(-1im*(1:N).*delta);

    u_grid = my_fft_plus(wplus) .+ my_fft_neg(wneg) .+ u[1];
    
    guC0rho = norm(abs.(u_grid),Inf);

    return guC0rho
end

#################################################################

function fmalpha(z) # f^(m)(alpha*z)

    if isa(z[1],Complex{Interval{BigFloat}})

        this=alpha.*z;
        
        for i=1:m
            
            this=f_fun(this);

        end            
        
        return this
        
    else 
        this=mid(alpha).*z;
        
        for i=1:m
            
            this=f_fun_num(this);
        end
        
        return this
    end
        

end


#############

function tfnalpha(z,n) # tf^(n)(alpha*z)

    if isa(z[1],Complex{Interval{BigFloat}})

        this=tl.*z;
        
        for i=1:n
            
            this=tf(this);

        end
           
        return this
        
    else 
        this=mid(tl).*z;
        
        for i=1:n
            
            this=tf_num(this);
        end
        
        return this
    end
        

end

#####################################


function func_alpha(z,n,j)   
    # also used for general terms in Jacobian
    # n is the order of the term (where term first apeared), j is starting number of iterations alphax inside of h, m=n+j
    # n>2
    if n==1
        return 1
    end
    
    if isa(z[1],Complex{Interval{BigFloat}})

        this=tl.*z;

        i=0;
        while i<j
            this=tf(this);
            i=i+1;
        end


        that=1;
        for k=1:n-1
            
            this=tf(this);
            that= dtf(this) .* that;
            
        end


        return that;

    else

        this=mid(tl).*z;

        i=0;
        while i<j
            this=tf_num(this);
            i=i+1;
        end

        that=1;
        for k=1:n-1
            
            this=tf_num(this);
            that= dtf_num(this) .* that;
            
        end
        return that;

    end

end

###########################################


function ffl(z)

    if isa(z[1],Complex{Interval{BigFloat}})

        this=alpha.*z;
        
        for i=1:m
            
            this=f_fun(this);

        end            
        
        return this
        
    else 
        this=mid(alpha).*z;
        
        for i=1:m
            
            this=f_fun_num(this);
        end
        
        return this
    end
        

end


function d_ffl(z)        # derivative of f^(m)(alpha*z) wrt to z

    this=alpha*z;
    that=df_fun(this);

    for i=2:m
        this=f_fun(this);
        that=df_fun(this).*that;        
    end

    ans=alpha*(df_fun(z));

    return ans

end


macro Name(arg)
    string(arg);
end

function d2_ffl(z,this_d2phi)  # second derivative of f^(m)(alpha*z) wrt to z
    global this
    this=z;
    thatz=@Name(this);
    this_d2phi=replace(this_d2phi,"thisz" =>"$thatz");
    ans= eval(Meta.parse(this_d2phi));
    return ans
end


function d2_fflz(z)

    this=alpha*z;
    fthis=f_fun(this);

    that= -alpha^2*df_fun(this).^2 .*d2f_fun(fthis)-alpha^2*df_fun(fthis).*d2f_fun(this)+alpha.*d2f_fun(z);

    return that

end




#########################################################


function tffl(z,n)

    if isa(z[1],Complex{Interval{BigFloat}})

        this=tl.*z;
        
        for i=1:n
            
            this=tf(this);

        end
           
        return this
        
    else 
        this=mid(tl).*z;
        
        for i=1:n
            
            this=tf_num(this);
        end
        
        return this
    end
        

end
############################################

function tail_func(x,N,tol1,tol2)

    if x<N
        return tol1
    else
        return tol2
    end

end

##################################


function opt_para_glob(g, rho, cste, alpha_max,tol, sub,nb)

    alphas = range(rho + tol,alpha_max, 20)

    thatalpha = interval(alphas[1])
    thisconst= cste(rho,thatalpha)
    thatval = thisconst*compute_sup_E(g,thatalpha,sub,nb)

    for k=2:20

        thisalpha = interval(alphas[k])
        thisconst= cste(rho,thisalpha)
        thisval = thisconst*compute_sup_E(g,thisalpha,sub,nb)
        if sup(thisval) < sup(thatval)
            thatval= copy(thisval)
            thatalpha=copy(thisalpha)
        end    

    end

    return thatalpha, thatval

end

#################################################################