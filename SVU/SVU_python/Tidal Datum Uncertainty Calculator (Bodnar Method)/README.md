Tidal Datum Uncertainty Calculator (Bodnar Method)
Overview
This script calculates the tidal datum uncertainty (standard deviation) between pairs of control and subordinate stations in the Pacific region. It processes station metadata and tidal datums to compute both 3-month and long-record (>= 12 months) uncertainties in a single, streamlined pipeline using the empirical Bodnar equations for the Pacific region (1981 NOAA Technical Report by A. Nicholas Bodnar, titled: "Estimating Accuracies of Tidal Datums from Short Term Observations, https://tidesandcurrents.noaa.gov/publications/NOAA_Technical_Report_NOS_COOPS_077.pdf, Page 9, Table 3). Note that within the report, Table 1 is for the East Coast and Table 2 is for the Gulf Coast.
Prerequisites
Ensure you have the following Python libraries installed:
pandas
numpy
pyproj
1. Computing the Timing Proxy (T)
The first step of the function calculates a timing proxy (measured in hours) that represents the tidal phase difference between stations.
Reads Data: It loads the ControlStation and SubordinateStation sheets from your chosen Excel file (like DartBuoy_Datums.xlsx or IOC_Datums.xlsx).
Cleans Data: It normalizes station IDs and ensures variables like Mean Range (MN), Diurnal High Water Inequality (DHQ), and Diurnal Low Water Inequality (DLQ) are formatted as numeric values.
Math: It calculates the ratio r = (DHQ - DLQ)/MN, and then multiplies it by a constant Khours (6.21) to get the absolute time difference Thr.
2. Computing Ellipsoidal Distance
Next, the script figures out the physical distance between each pair of control and subordinate stations.
Reads Metadata: It pulls latitude and longitude coordinates for all stations from the separate metadata file (Obs_Pacific.csv).
Calculates Distance: It maps the coordinates to the station pairs generated in the first step and uses the pyproj library (WGS-84 ellipsoid) to calculate the exact distance in kilometers.
Conversions: It converts the kilometers into nautical miles (dist_nm) and also calculates the square root of that distance (D_sqrt_nm), which is a required input for the 3-month Bodnar formula.
3. Calculating Uncertainty (The 3-Month Bodnar Method)
The script then gathers all the necessary pieces to calculate the 3-month datum uncertainty.
Gathers Mean Range (MN): It extracts the MN values for both the control stations (Mnc) and subordinate stations (MNs) from the original Excel file and converts them from meters to feet.
Computes the Mean Range Ratio (MNR): It calculates the absolute difference in Mean Range between the two stations, divided by the control station's Mean Range.
Applies the Bodnar Formula: It plugs the three calculated variables into the empirical 3-month equation for the Pacific:
 
S (ft) = 0.0043(dt)+ 0.0036√D+ 0.0255 |MNc-MNs/MNc| + 0.029
Where dt is Time difference (Hours), D is Distance (Nautical Miles), MNc is Tidal range in control station, MNs is Tidal range in subordinate station
4. Calculating Uncertainty (Long-Record IOC/ UHSLC Stations)
Immediately after, the script calculates the long-record uncertainty for stations with 12 months or more of data using the alternate regression values.
Calculates SRSMN: It calculates a new variable, which is the square root of the sum of the control and subordinate mean ranges (in feet).
Applies the Long-Record Formula: It calculates the uncertainty using the long-record regression values:
S (ft) = 0.0045√ (MNc + MNs) + 0.0128 |MNc-MNs/MNc| + 0.025
5. Final Comprehensive Output
The script converts both the 3-month and long-record uncertainties back into meters. It packages all variables (IDs, distances, time proxies, mean ranges, and both final uncertainty values) together and exports them directly into one single, comprehensive CSV file.
