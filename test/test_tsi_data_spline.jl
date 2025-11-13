import Plots
import Dates

tsi_data_spline = TSIDataSpline(FT)
monthly_dates, tsi_data = Insolation._get_tsi_data()

# Evaluate exactly at each date
for (date, tsi_val) in zip(monthly_dates, tsi_data)
    @test evaluate(tsi_data_spline, date) ≈ tsi_val
end

@test isnan(evaluate(tsi_data_spline, Dates.DateTime(0)))
@test isnan(evaluate(tsi_data_spline, Dates.DateTime(4200)))

dates = collect(first(monthly_dates):Dates.Day(1):last(monthly_dates))
tsi = evaluate.(tsi_data_spline, dates)
tsi_data_spline_plot = Plots.plot(dates, tsi)
Plots.savefig(tsi_data_spline_plot, "tsi_data_spline.png")

# Solar cycle 24
# Approximate dates for solar min and max
solar_min_date = Dates.DateTime(2007, 12, 15, 12)
solar_max_date = Dates.DateTime(2014, 3, 15, 12)
solar_min = evaluate(tsi_data_spline, solar_min_date)
solar_max = evaluate(tsi_data_spline, solar_max_date)

@test solar_max - solar_min > 1

dates = collect(solar_min_date:Dates.Day(1):solar_max_date)
tsi = evaluate.(tsi_data_spline, dates)
tsi_data_spline_plot = Plots.plot(dates, tsi)
Plots.savefig(tsi_data_spline_plot, "solar_cycle_24.png")
