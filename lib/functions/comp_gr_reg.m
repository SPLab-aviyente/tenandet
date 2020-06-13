function val = comp_gr_reg(La, Phi)

La = Lunfold(La);
val = trace(La'*La*Phi);
end