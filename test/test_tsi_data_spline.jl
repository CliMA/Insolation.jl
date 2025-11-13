import Artifacts

tsi_data_spline = TSIDataSpline()
monthly_dates, tsi_data = Insolation._get_tsi_data()

# Evaluate exactly at each date
for (date, tsi_val) in zip(monthly_dates, tsi_data)
    @test evaluate(tsi_data_spline, date) ≈ tsi_val
end

# Evaluate between two different dates
prev_date = monthly_dates[1]
next_date = monthly_dates[2]
midway = prev_date + (next_date - prev_date) / 2

@test evaluate(tsi_data_spline, midway) ≈
      (evaluate(tsi_data_spline, prev_date) + evaluate(tsi_data_spline, next_date)) / 2
