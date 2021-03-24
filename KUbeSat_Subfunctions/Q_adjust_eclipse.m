function QQ = Q_adjust_eclipse(Qmax)

n = length(Qmax);
% check for eclipses
eclipse = zeros(n,1);
Ymax = Qmax(1,:);
% check for eclipse at the start time
if (abs(Ymax(2)-Ymax(1)) < 0.000001)
    N_eclipse = 1; 
else
    N_eclipse = 0;
end
for i = 2:n
    if (abs(Ymax(i)-Ymax(i-1)) < 0.000001)
        eclipse(i) = 1;
        % check if this is a new eclipse
        if (i > 2 && abs(Ymax(i-1)-Ymax(i-2)) > 0.000001)
            N_eclipse = N_eclipse + 1;
        end
    end
end
% adjust array in case of eclipse(s)
if (N_eclipse > 0)
    E_start = zeros(N_eclipse,1); E_end = zeros(N_eclipse,1);
    if (eclipse(1) > 0)
        E_start(1) = 1;
    else
        for i = 2:n
            if (eclipse(i) > 0 && eclipse(i-1) < 1)
                E_start(1) = i-1;
                break;
            end
        end
    end
    for i = E_start(1)+1:n
        if (eclipse(i) < 1 && eclipse(i-1) > 0)
            E_end(1) = i;
            break;
        end
    end
    if (N_eclipse > 1)
        for j = 2:N_eclipse
            for i = E_end(j-1)+1:n
                if (eclipse(i) > 0 && eclipse(i-1) < 1)
                    E_start(j) = i-1;
                    break;
                end
            end
            for i = E_start(j)+1:n
                if (eclipse(i) < 1 && eclipse(i-1) > 0)
                    E_end(j) = i;
                    break;
                end
            end
        end
    end
    if (eclipse(n) > 0)
        E_end(N_eclipse) = n;
    end
    % adjust quaternions angles
    for j = 1:N_eclipse
        QL = Qmax(:,E_start(j));
        QR = Qmax(:,E_end(j));
        nn = E_end(j) - E_start(j);
        for i = (E_start(j)+1):(E_end(j)-1)
            Qmax(:,i) = QL + (QR-QL).*(i-E_start(j))./nn;
        end
    end     
   
end

QQ = Qmax;

end