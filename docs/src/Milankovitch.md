# Milankovitch Cycles

## Variations in orbital parameters
```@example
using Insolation #hide
using Plots #hide

dt = collect(-300e3:100:300e3);
ϖ_spline, γ_spline, e_spline = orbital_params_spline();
ϖ, γ, e = ϖ_spline.(dt), γ_spline.(dt), e_spline.(dt);

p1 = plot(dt ./ (1e3), sin.(ϖ), legend=false);
ylabel!("sin(ϖ)");
p2 = plot(dt ./ (1e3), γ, legend=false);
ylabel!("γ");
p3 = plot(dt ./ (1e3), e, legend=false);
ylabel!("e");
xlabel!("time (kY)")
plot(p1, p2, p3, layout = grid(3,1), size=(600,400), dpi=250);
savefig("orbital_params.png")
```
![](orbital_params.png)
