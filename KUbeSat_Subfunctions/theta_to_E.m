function E = theta_to_E(e,theta)
  E = 2*atan(sqrt((1-e)/(1+e))*tan(0.5*theta));
  E = mod(E, 2*pi);
end