using LsqFit

  # a two-parameter exponential model
  # x: array of independent variables
  # p: array of model parameters
  model(x, p) = p[1]*exp.(-x.*p[2])
  
  # some example data
  # xdata: independent variables
  # ydata: dependent variable
  xdata = range(0, stop=10, length=20)
  ydata = model(xdata, [100.0 2.0]) + 0.01*randn(length(xdata))
  p0 = [0.5, 0.5]
  
  fit = curve_fit(model, xdata, ydata, p0)
 a = fit.resid
 dof(fit)
 coef(fit)

 sum(a.^2)