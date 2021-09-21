function rho = HighAltDrag(h)
    hKm = h/1000;
    
    if(hKm > 600)
        hKm = 600;
    end
    
    h0 = -1;
    rho0 = -1;
    H = -1;
    
    if (hKm < 25) 
        h0 = 0;
        rho0 = 1.225;
        H = 7.249;
    elseif (hKm < 30)
        h0 = 25;
        rho0 = 3.899e-2;
        H = 6.349;
    elseif (hKm < 40)
        h0 = 30;
        rho0 = 1.774e-2;
        H = 6.682;
    elseif (hKm < 50)
        h0 = 40;
        rho0 = 3.972e-3;
        H = 7.554;
    elseif (hKm < 60)
        h0 = 50;
        rho0 = 1.057e-3;
        H = 8.382;
    elseif (hKm < 70)
        h0 = 60;
        rho0 = 3.206e-4;
        H = 7.714;
    elseif (hKm < 80)
        h0 = 70;
        rho0 = 8.770e-5;
        H = 6.549;
    elseif (hKm < 90)
        h0 = 80;
        rho0 = 1.905e-5;
        H = 5.799;
    elseif (hKm < 100)
        h0 = 90;
        rho0 = 3.396e-6;
        H = 5.382;
    elseif (hKm < 110)
        h0 = 100;
        rho0 = 5.297e-7;
        H = 5.877;
    elseif (hKm < 120)
        h0 = 100;
        rho0 = 9.661e-7;
        H = 7.263;
    elseif (hKm < 130)
        h0 = 120;
        rho0 = 2.438e-8;
        H = 9.473;
    elseif (hKm < 140)
        h0 = 130;
        rho0 = 8.484e-9;
        H = 12.636;
    elseif (hKm < 150)
        h0 = 140;
        rho0 = 3.845e-9;
        H = 16.149;
    elseif (hKm < 180)
        h0 = 150;
        rho0 = 2.070e-9;
        H = 22.523;
    elseif (hKm < 200)
        h0 = 180;
        rho0 = 5.464e-10;
        H = 29.740;
    elseif (hKm < 250)
        h0 = 200;
        rho0 = 2.789e-10;
        H = 37.105;
    elseif (hKm < 300)
        h0 = 250;
        rho0 = 7.248e-11;
        H = 45.546;
    elseif (hKm < 350)
        h0 = 300;
        rho0 = 2.418e-11;
        H = 53.628;
    elseif (hKm < 400)
        h0 = 350;
        rho0 = 9.518e-12;
        H = 53.298;
    elseif (hKm < 450)
        h0 = 400;
        rho0 = 3.725e-12;
        H = 58.515;
    elseif (hKm < 500)
        h0 = 450;
        rho0 = 1.585e-12;
        H = 60.828;
    elseif (hKm <= 600)
        h0 = 500;
        rho0 = 6.967e-13;
        H = 63.822;
    end

    rho = rho0*exp(-(hKm - h0)/H);
    
end