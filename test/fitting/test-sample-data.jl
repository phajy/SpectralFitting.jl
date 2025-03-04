using SpectralFitting
using Test

# might have to download the model data first
SpectralFitting.download_model_data(XS_Laor)
SpectralFitting.download_model_data(PhotoelectricAbsorption)

function prepare_data!(data, low, high)
    drop_bad_channels!(data)
    regroup!(data)
    normalize!(data)
    mask_energies!(data, low, high)
    data
end

# path to the data directory
data1 = SpectralFitting.XmmData(joinpath(testdir, "xmm/pn_spec_grp.fits"))
prepare_data!(data1, 0.8, 10.0)

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())
# todo: check the chi2 calculations
@test result.χ2 ≈ 484.5 atol = 0.1

@test result.u ≈ [8.9215e-05, 6.5673, 0.010446, 1.8422, 0.13760] rtol = 1e-3
xspec_u = [9.25e-5, 6.559, 1.048e-2, 1.84862, 0.140182]
@test result.u ≈ xspec_u rtol = 1e-2

data2 = NuStarData(joinpath(testdir, "nustar/nu60001047002A01_sr_grp_simple.pha"))
prepare_data!(data2, 2.0, 14.0)

prob = FittingProblem(model, data2)
result = fit(prob, LevenbergMarquadt())

@test result.χ2 ≈ 191.0 atol = 0.1
@test result.u ≈ [
    0.0002059507411466161,
    6.499205973061382,
    0.01989636801877023,
    1.9972614257180774,
    0.33187746370113763,
] rtol = 1e-2

# test joint
prob = FittingProblem(model => data1, model => data2)
result = fit(prob, LevenbergMarquadt())

@test result[1].χ2 ≈ 484.546 atol = 0.1
@test result[2].χ2 ≈ 190.999 atol = 0.1

# with binding
bind!(prob, :K_1)
result = fit(prob, LevenbergMarquadt())

@test result[1].χ2 ≈ 485.0180163473964 atol = 0.1
@test result[2].χ2 ≈ 202.1704312745641 atol = 0.1

# todo: with background subtraction
data1_nobkg = deepcopy(data1)
set_units!(data1_nobkg, u"counts")
subtract_background!(data1_nobkg)
set_units!(data1_nobkg, u"counts / (s * keV)")

prob = FittingProblem(model, data1_nobkg)
result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 496.0927976691264 atol = 0.1
@test result.u ≈ [9.2113e-05, 6.5597, 0.010478, 1.8483, 0.13960] rtol = 1e-3

# try different domain entirely
set_domain!(data1, collect(range(0.01, 20.0, 1000)))

model = PhotoelectricAbsorption() * XS_PowerLaw() + XS_Laor()

# construct the model and data problem
prob = FittingProblem(model, data1)

result = fit(prob, LevenbergMarquadt())
# these have been checked and are the same as XSPEC
@test result.χ2 ≈ 466.6477135605372 atol = 0.1
@test result.u ≈ [
    9.098183992394181e-5,
    6.574126894811102,
    0.01052894049667675,
    1.8472152777352222,
    0.13435297425331832,
] rtol = 1e-3
