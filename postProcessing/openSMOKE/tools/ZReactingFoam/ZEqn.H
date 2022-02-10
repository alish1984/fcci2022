// mixture fraction equation assuming Le = 1, Sc = 1
fvScalarMatrix ZEqn
(
    fvm::ddt(rho, Z)
  + fvm::div(phi, Z)
  - fvm::laplacian(turbulence->muEff(), Z)
 ==
    fvOptions(rho, Z)
);
ZEqn.solve();
Info << "mixture fraction min/max  = "
<< min(Z).value() << ", "
<< max(Z).value() << endl;
Z.max(0.0);
Z.min(1.0);
Info << "mixture fraction min/max  = "
<< min(Z).value() << ", "
<< max(Z).value() << endl;
