using DifferentialEquations
using Plots
using Latexify

# -------- Dynamical Properties ---------
# 1. EQ model
w = 10.0
f_eq = @ode_def EQ_model begin
    du = -K_on* u* w + K_off* c
    dc = K_on* u* w -K_off* c
    dm = f*c -m/τ_m
end K_on K_off f τ_m


u0 = [1.0, 0.0, 0.0]
tspan = (0.0,1.0)
p = [10.0, 8.0, 4.0, 2.0    ]
prob = ODEProblem(f_eq,u0,tspan,p)
sol = solve(prob, progress = true)
plot(sol)


# 2. FM model
w = 10.0
f_FM = @ode_def FM_model begin
    dA = w*(1- A -R)/(τ_I*(k+w)) - A/τ_A
    dR = A/τ_A - R/τ_R
    dm = f*A -m/τ_m
end τ_I τ_A τ_R f τ_m k

u0 = [1.0, 0.0, 0.0]
tspan = (0.0,1e2)
p_FM = [1.0, 15.0, 40.0, 2.0, 2.0, 0.8 ]
prob_FM = ODEProblem(f_FM,u0,tspan,p_FM)

sol_FM = solve(prob_FM, progress = true)
plot(sol_FM,vars=[:m], label = "FM_model mRNA level") # just ploting the mRNA level

# 3. DM model
w = 10.0
f_DM = @ode_def DM_model begin
    dA = (1- A -R)/τ_I - (k/(k+w))*A/τ_A
    dR = (k/(k+w))*A/τ_A - R/τ_R
    dm = f*A -m/τ_m
end τ_I τ_A τ_R f τ_m k

u0 = [1.0, 0.0, 0.0]
tspan = (0.0,1e2)
p_DM = [1.0, 15.0, 40.0, 2.0, 2.0, 0.8 ]
prob_DM = ODEProblem(f_DM,u0,tspan,p_DM)

sol_DM = solve(prob_DM, progress = true)
plot!(sol_DM,vars=[:m],label = "DM_model mRNA level") # just ploting the mRNA level
yaxis!("mRNA Level")
title!("Model Comparison : Dynamical Response when [TF] = 10")









##############################
####   ----- No ϵ ------   ###
##############################

# ------- Paper models
# Two states telegraph model
τ_A, τ_I, k = 10, 8, 10
V_FM(w) = τ_A*(1-k/(k+w))/(τ_I + τ_A*(1-k/(k+w)))
V_DM(w) = τ_A/(τ_I + τ_A*(1-k/(k+w)))
plot(V_FM,0:50)
plot!(V_DM,0:50)

# Refractory cycling model
τ_A, τ_I, τ_R, k = 15, 1, 40, 100
τ_total = τ_A + τ_I + τ_R
R_FM(w) = (τ_A/τ_total)*(w/(k*τ_I/τ_total + w))

τ_A, τ_I, τ_R, k = 1, 15, 40, 100
τ_total = τ_A + τ_I + τ_R
R_DM(w) = (τ_A/τ_total)*((k+w)/(k + w*τ_A/τ_total))

# normalization of mRNA level
x = 1:1e5
R_FM_norm = (R_FM.(x) - minimum(R_FM.(x)))/(maximum(R_FM.(x)) - minimum(R_FM.(x)))
R_DM_norm = (R_DM.(x) - minimum(R_DM.(x)))/(maximum(R_DM.(x)) - minimum(R_DM.(x)))
plot(log10(x), R_FM_norm,  label  = "Refractory FM model",line=(4,:black))
plot!(log10(x), R_DM_norm, label  = "Refractory DM model",line=(4,:dash))
title!("Model Comparison : Dose Response")
xaxis!("Log[TF]")
yaxis!("Normalized mRNA Level")
savefig("paper_model.png")

# ------ -------- My modol ---------- ----------
# -------3 states simplified models without ϵ

# 1. With Cooperativity | changing Hill coef
# Refractory cycling model
τ_A, τ_I, τ_R, k = 15, 1, 40, 100
τ_total = τ_A + τ_I + τ_R
R_FM(w) = (τ_A/τ_total)*(w^2/(k^2*τ_I/τ_total + w^2))

τ_A, τ_I, τ_R, k = 1, 15, 40, 100
τ_total = τ_A + τ_I + τ_R
R_DM(w) = (τ_A/τ_total)*((k^2+w^2)/(k^2 + w^2*τ_A/τ_total))

# normalization of mRNA level
x = 1:1e4
R_FM_norm = (R_FM.(x) - minimum(R_FM.(x)))/(maximum(R_FM.(x)) - minimum(R_FM.(x)))
R_DM_norm = (R_DM.(x) - minimum(R_DM.(x)))/(maximum(R_DM.(x)) - minimum(R_DM.(x)))
plot(log10(x), R_FM_norm,  label  = "Refractory FM model",line=(4,:black))
plot!(log10(x), R_DM_norm, label  = "Refractory DM model",line=(4,:dash))
title!("Model Comparison : Dose Response (Hill Coef = 2)")
xaxis!("Log[TF]")
yaxis!("Normalized mRNA Level")


# Animation of FM and DM model Dose Responses with varied Hill Coef
anim = @animate for i=1:0.5:20
    τ_A, τ_I, τ_R, k = 15, 1, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_FM(w) = (τ_A/τ_total)*(w^i/(k^i*τ_I/τ_total + w^i))

    τ_A, τ_I, τ_R, k = 1, 15, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_DM(w) = (τ_A/τ_total)*((k^i+w^i)/(k^i + w^i*τ_A/τ_total))

    x = 1:1e4
    R_FM_norm = (R_FM.(x) - minimum(R_FM.(x)))/(maximum(R_FM.(x)) - minimum(R_FM.(x)))
    R_DM_norm = (R_DM.(x) - minimum(R_DM.(x)))/(maximum(R_DM.(x)) - minimum(R_DM.(x)))
    plot(log10(x), R_FM_norm,  label  = "Refractory FM model",line=(4,:black))
    plot!(log10(x), R_DM_norm, label  = "Refractory DM model",line=(4,:dash))
    title!("Model Comparison : Dose Response (Hill Coef = 1 ~ 20)")
    xaxis!("Log[TF]")
    yaxis!("Normalized mRNA Level")
end
gif(anim, "/tmp/paper_model_varied_Hill_coef.gif", fps = 30)




##############################
####  ----- With ϵ ------  ###
##############################

                            # --------- Paper -----------


# 1.
# 3 states FM DM models with modifed ϵ  |  without Cooperativity
# --------- Paper models (Refractory cycling model)  ----♌ ϵ_FM and ϵ_DM have to be different.


ϵ_FM, τ_A, τ_I, τ_R, k = 0.5, 15, 1, 40, 100
τ_total = τ_A + τ_I + τ_R
R_FM_ϵ(w) = (τ_A/τ_total)*((ϵ_FM*k + w)/(k*(τ_I + ϵ_FM*(τ_A + τ_R))/τ_total + w))

ϵ_DM, τ_A, τ_I, τ_R, k = 0.5, 1, 15, 40, 100
τ_total = τ_A + τ_I + τ_R
R_DM_ϵ(w) = (τ_A/τ_total)*((k+w)/(k*(τ_A + (τ_I + τ_R)/ϵ_DM)/τ_total + w))

# normalization of mRNA level
x = 1:1e4
R_FM_ϵ_norm = (R_FM_ϵ.(x) - minimum(R_FM_ϵ.(x)))/(maximum(R_FM_ϵ.(x)) - minimum(R_FM_ϵ.(x)))
R_DM_ϵ_norm = (R_DM_ϵ.(x) - minimum(R_DM_ϵ.(x)))/(maximum(R_DM_ϵ.(x)) - minimum(R_DM_ϵ.(x)))
plot(log10(x), R_FM_ϵ_norm,  label  = "Refractory FM model",line=(4,:black))
plot!(log10(x), R_DM_ϵ_norm, label  = "Refractory DM model",line=(4,:dash))
title!("Model Comparison : Dose Response")
xaxis!("Log[TF]")
yaxis!("Normalized mRNA Level")
savefig("paper_model_ϵ.png")

# 2.
# 3 states FM DM models with modifed ϵ  |  with Cooperativity (changing Hill coef)
# --------- Modified Paper models (Refractory cycling model)


anim = @animate for i=1:0.5:20
    ϵ_FM, τ_A, τ_I, τ_R, k = 0.5, 15, 1, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_FM_ϵ(w) = (τ_A/τ_total)*((ϵ_FM*k^i + w^i)/(k^i*(τ_I + ϵ_FM*(τ_A + τ_R))/τ_total + w^i))

    ϵ_DM, τ_A, τ_I, τ_R, k = 2, 1, 15, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_DM_ϵ(w) = (τ_A/τ_total)*((k^i+w^i)/(k^i*(τ_A + (τ_I + τ_R)/ϵ_FM)/τ_total + w^i))

    # normalization of mRNA level
    x = 1:1e4
    R_FM_ϵ_norm = (R_FM_ϵ.(x) - minimum(R_FM_ϵ.(x)))/(maximum(R_FM_ϵ.(x)) - minimum(R_FM_ϵ.(x)))
    R_DM_ϵ_norm = (R_DM_ϵ.(x) - minimum(R_DM_ϵ.(x)))/(maximum(R_DM_ϵ.(x)) - minimum(R_DM_ϵ.(x)))
    plot(log10(x), R_FM_ϵ_norm,  label  = "Refractory FM model",line=(4,:black))
    plot!(log10(x), R_DM_ϵ_norm, label  = "Refractory DM model",line=(4,:dash))
    title!("Model Comparison : Dose Response")
    xaxis!("Log[TF]")
    yaxis!("Normalized mRNA Level")
end
gif(anim, "/tmp/ϵ_cooperativity.gif", fps = 30)






                                # --------- My modifed -----------

# 3 states FM DM models with modifed ϵ  |  without Cooperativity
# --------- My Modified models (Refractory cycling model)    |   τ_IB = ϵ_FM_m * τ_I,  τ_IU = τ_I
ϵ_FM_m = 0.5

τ_A, τ_I, τ_R, k = 15, 1, 40, 100
τ_total = τ_A + τ_I + τ_R
R_FM_ϵ_m(w) = (τ_A/τ_total)*((k + w/ϵ_FM_m)/(k + w*(τ_I + (τ_A + τ_R)/ϵ_FM_m)/τ_total ))

τ_A, τ_I, τ_R, k = 1, 15, 40, 100
τ_total = τ_A + τ_I + τ_R
R_DM_ϵ_m(w) = (τ_A/τ_total)*((k+w)/(k + w*(τ_A + (τ_I + τ_R)/ϵ_FM_m)/τ_total ))

# normalization of mRNA level
x = 1:1e4
R_FM_ϵ_m_norm = (R_FM_ϵ_m.(x) - minimum(R_FM_ϵ_m.(x)))/(maximum(R_FM_ϵ_m.(x)) - minimum(R_FM_ϵ_m.(x)))
R_DM_ϵ_m_norm = (R_DM_ϵ_m.(x) - minimum(R_DM_ϵ_m.(x)))/(maximum(R_DM_ϵ_m.(x)) - minimum(R_DM_ϵ_m.(x)))
plot(log10(x), R_FM_ϵ_m_norm,  label  = "Refractory FM model",line=(4,:black))
plot!(log10(x), R_DM_ϵ_m_norm, label  = "Refractory DM model",line=(4,:dash))
title!("Model Comparison : Dose Response")
xaxis!("Log[TF]")
yaxis!("Normalized mRNA Level")
savefig("my_model_ϵ.png")


# 3 states FM DM models with modifed ϵ  |  with Cooperativity (changing Hill coef)
# --------- My Modified models (Refractory cycling model)    |   τ_AB = ϵ_DM_m * τ_A,  τ_AU = τ_A
ϵ_FM_m = 0.5

anim = @animate for i=1:0.5:20
    τ_A, τ_I, τ_R, k = 15, 1, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_FM_ϵ_m(w) = (τ_A/τ_total)*((k^i + w^i/ϵ_FM_m)/(k^i + w^i*(τ_I + (τ_A + τ_R)/ϵ_FM_m)/τ_total ))

    τ_A, τ_I, τ_R, k = 1, 15, 40, 100
    τ_total = τ_A + τ_I + τ_R
    R_DM_ϵ_m(w) = (τ_A/τ_total)*((k^i+w^i)/(k^i + w^i*(τ_A + (τ_I + τ_R)/ϵ_FM_m)/τ_total ))

    # normalization of mRNA level
    x = 1:1e4
    R_FM_ϵ_m_norm = (R_FM_ϵ_m.(x) - minimum(R_FM_ϵ_m.(x)))/(maximum(R_FM_ϵ_m.(x)) - minimum(R_FM_ϵ_m.(x)))
    R_DM_ϵ_m_norm = (R_DM_ϵ_m.(x) - minimum(R_DM_ϵ_m.(x)))/(maximum(R_DM_ϵ_m.(x)) - minimum(R_DM_ϵ_m.(x)))
    plot(log10(x), R_FM_ϵ_m_norm,  label  = "Refractory FM model",line=(4,:black))
    plot!(log10(x), R_DM_ϵ_m_norm, label  = "Refractory DM model",line=(4,:dash))
    title!("Model Comparison : Dose Response")
    xaxis!("Log[TF]")
    yaxis!("Normalized mRNA Level")
end
gif(anim, "/tmp/ϵ_m_cooperativity.gif", fps = 30)





























































# @gif for i=1:20
#     τ_A, τ_I, τ_R, k = 2, 1, 2, 100
#     τ_total = τ_A + τ_I + τ_R
#     R_FM(w) = (τ_A/τ_total)*(w^i/(k^i*τ_I/τ_total + w^i))
#
#     τ_A, τ_I, τ_R, k = 4, 1, 2, 100
#     τ_total = τ_A + τ_I + τ_R
#     R_DM(w) = (τ_A/τ_total)*((k^i+w^i)/(k^i + w^i*τ_A/τ_total))
#
#     R_FM_norm = (R_FM.(x) - minimum(R_FM.(x)))/(maximum(R_FM.(x)) - minimum(R_FM.(x)))
#     R_DM_norm = (R_DM.(x) - minimum(R_DM.(x)))/(maximum(R_DM.(x)) - minimum(R_DM.(x)))
#     plot(log10(x), R_FM_norm,  label  = "Refractory FM model",line=(4,:black))
#     plot!(log10(x), R_DM_norm, label  = "Refractory DM model",line=(4,:dash))
# end every 1
