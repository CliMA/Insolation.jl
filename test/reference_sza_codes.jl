using Dates

function inst_sza(date, tz, dLongitude, dLatitude)
    dEarthMeanRadius = 6371.01
    dAstronomicalUnit = 149597890

    ###################################################################
    # Calculate difference in days between the current Julian Day
    # and JD 2451545.0, which is noon 1 January 2000 Universal Time
    ###################################################################
    julian_day_abs = datetime2julian(date)
    dElapsedJulianDays = julian_day_abs - 2451545.0

    ###################################################################
    # Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    # ecliptic in radians but without limiting the angle to be less than 2*Pi
    # (i.e., the result may be greater than 2*Pi)
    ###################################################################
    dOmega = 2.1429 - 0.0010394594 * dElapsedJulianDays
    dMeanLongitude = 4.8950630 + 0.017202791698 * dElapsedJulianDays  # Radians
    dMeanAnomaly = 6.2400600 + 0.0172019699 * dElapsedJulianDays
    dEclipticLongitude = dMeanLongitude + 0.03341607 * sin(dMeanAnomaly) + 0.00034894 * sin(2. * dMeanAnomaly) - 0.0001134 - 0.0000203 * sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * dElapsedJulianDays + 0.0000396 * cos(dOmega)

    ###################################################################
    # Calculate celestial coordinates ( right ascension and declination ) in radians
    # but without limiting the angle to be less than 2*Pi (i.e., the result may be
    # greater than 2*Pi)
    ###################################################################
    dSin_EclipticLongitude = sin(dEclipticLongitude)
    dY = cos(dEclipticObliquity) * dSin_EclipticLongitude
    dX = cos(dEclipticLongitude)
    dRightAscension = atan(dY, dX)
    if dRightAscension < 0.0
        dRightAscension = dRightAscension + 2.0 * π
    end
    dDeclination = asin(sin(dEclipticObliquity) * dSin_EclipticLongitude)

    ###################################################################
    # Calculate local coordinates ( azimuth and zenith angle ) in degrees
    ###################################################################
    dHours, dMinutes, dSeconds = Dates.hour(date), Dates.minute(date), Dates.second(date)
    dDecimalHours = dHours + (dMinutes + dSeconds / 60.) / 60. - tz
    dGreenwichMeanSiderealTime = 6.6974243242 + 0.0657098283 * dElapsedJulianDays + dDecimalHours
    dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15. + dLongitude) * (π / 180.)
    dHourAngle = dLocalMeanSiderealTime - dRightAscension
    #println(dLocalMeanSiderealTime, "\t old")
    println(dHourAngle, "\t old")
    dLatitudeInRadians = dLatitude * (π / 180.)
    dCos_Latitude = cos(dLatitudeInRadians)
    dSin_Latitude = sin(dLatitudeInRadians)
    dCos_HourAngle = cos(dHourAngle)
    dZenithAngle = acos(dCos_Latitude * dCos_HourAngle * cos(dDeclination) + sin(dDeclination) * dSin_Latitude)

    # Parallax Correction
    dParallax = (dEarthMeanRadius / dAstronomicalUnit) * sin(dZenithAngle)
    dZenithAngle = (dZenithAngle + dParallax) / (π / 180.)
    
    # Set to a max of 90 deg.
    if dZenithAngle > 90.0
        dZenithAngle = 90.0
    end

    return dZenithAngle
end

function daily_sza(days_since_equinox, lat, gamma, pi, e)
    # convert inputs from degrees to radians
    lat = lat * (2*π/360.0)
    gamma = gamma * (2*π/360.0)
    pi = pi * (2*π/360.0)
    
    # constants
    Ya = 365.26 # days
    
    # step 1, calculate the mean anomaly at vernal equinox
    beta = sqrt(1-e^2)
    M_VE = -pi + (e + e^3/4.0)*(1 + beta)*sin(pi)
    
    # step 2, calculate the mean anomaly
    M = (2*π*(days_since_equinox))/(Ya) + M_VE
    
    # step 3, calculate the true anomaly
    A = M + (2*e - e^3/4.0)*sin(M)
    
    # step 4, calculate the distance to the sun
    d = (1-e^2)/(1+e*cos(A))
    
    # step 5, calculate the solar longitude
    L_s = A + pi
    
    # step 6, calculate the declination angle 
    delta = asin(sin(gamma) * sin(L_s))
    
    # step 7, calculate the sunrise/sunset angle
    T = tan(lat) * tan(delta)
    if T >= 1
        eta_d = π
    elseif T <= -1
        eta_d = 0.0
    else
        eta_d = acos(-1 * T)
    end
    
    # step 8, calculate the daily averaged cos(zenith angle)
    c1 = eta_d*sin(lat)*sin(delta)
    c2 = cos(lat)*cos(delta)*sin(eta_d)
    cosbar = (1/π)*(c1+c2)

    sza = acos(cosbar)
    sza = rad2deg(sza)

    return sza
end