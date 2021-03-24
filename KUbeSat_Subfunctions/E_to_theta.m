function th = E_to_theta(e,E)
  th = 2*atan(sqrt((1+e)/(1-e))*tan(0.5*E));
  th = mod(th,2*pi);
end