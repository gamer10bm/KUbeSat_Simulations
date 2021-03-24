% SOLUTION TO KEPLER'S EQUATION FOR ELLIPSE (Laguerre's Method)
% By: Brian Kaplinger, Ph.D.
% Changed: February 21, 2014
function E = kepler_ellipse_L(e, M)
  % Input Variables:
  %   e: orbit eccentricity (e < 1)
  %   M: mean anomaly (radians)
  % Output Variables:
  %   E: eccentric anomaly (radians)

  tol = 1e-12;  % set tolerance for solution value
  maxiter = 1000;  % set maximum number of iterations
  
  % set initial guess (in radians)
  % from (4.29), Chobotov, 2002
  B = cos(e) - (0.5*pi - e)*sin(e);
  E = M + e*sin(M)/(B + M*sin(e));

  % begin iteration
  for i = 1:maxiter
      % compute function value
      % Kepler's equation M = E - e*sin(E)
      % f(E) = E - e*sinE(E) - M
      f = E - e*sin(E) - M;
      
      % compute derivative value at current guess
      df = 1 - e*cos(E);
      
      % compute second derivative value
      ddf = e*sin(E);
      
      % set new iteration value
      E_new = E - 5*f/(df+2*sign(df)*sqrt(abs(4*df^2-5*f*ddf)));
      
      % check convergence
      if (abs(f) < tol)
          E = E_new;
          break;
      else
          E = E_new;
      end
  end
  
  % display answer
  %disp(['Eccentricity is ', num2str(e), '.']);
  %disp(['Mean Anomaly is ', num2str(M), ' radians.']);
  %disp(['Eccentric Anomaly is ', num2str(E), ' radians.']);
  %disp(['Solution took ', num2str(i), ' iterations.']);
end