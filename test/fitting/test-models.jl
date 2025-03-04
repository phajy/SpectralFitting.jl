using SpectralFitting
using Test


model = BremsStrahlung(kT = FitParam(3.0)) + GaussianLine()

domain = collect(range(0.1, 10.0, 200))
data = simulate(domain, model, var = 1e-5, seed = 5)

model = BremsStrahlung() + GaussianLine()

# smoke test
result = fit(FittingProblem(model => data), LevenbergMarquadt())
