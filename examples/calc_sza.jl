using Dates
using Insolation.SZA

function main()
    PST = -7.0
    California = [-105.0, 40.0]

    tz = PST
    lon, lat = California

    date = Dates.DateTime(2020,3,21,12,0,0)
    sza = instantaneous_zenith_angle(date, tz, lon, lat)
    println(rad2deg(sza))

    date = Dates.DateTime(2020,6,21,12,0,0)
    sza = instantaneous_zenith_angle(date, tz, lon, lat)
    println(rad2deg(sza))

    date = Dates.DateTime(2010,9,21,12,0,0)
    sza = instantaneous_zenith_angle(date, tz, lon, lat)
    println(rad2deg(sza))

    date = Dates.DateTime(2010,12,21,12,0,0)
    sza = instantaneous_zenith_angle(date, tz, lon, lat)
    println(rad2deg(sza))

    date = Dates.now()
    sza = instantaneous_zenith_angle(date, tz, lon, lat)
    println(rad2deg(sza))

    szabar = daily_zenith_angle(date, lon, lat)
    println(rad2deg(szabar))
end

main()