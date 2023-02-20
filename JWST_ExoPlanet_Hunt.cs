using System;
using System.Collections.Generic;
using System.Net;
using System.Net.Security;
using TLS;

namespace JWST_Exoplanet_Search
{
    class Program
    {
        static void Main(string[] args)
        {
            // Define the time range for the JWST observations
            DateTime startTime = new DateTime(2023, 2, 18, 0, 0, 0, DateTimeKind.Utc);
            DateTime endTime = new DateTime(2023, 2, 18, 1, 0, 0, DateTimeKind.Utc);

            // Query the JWST for its position during the observation time
            double[] jwstPos = GetJWSTPosition(startTime);

            // Define the TLS parameters
            int transitMinDuration = 1; // minimum transit duration in hours
            int transitMaxDuration = 5; // maximum transit duration in hours
            double transitDepthMin = 0.001; // minimum transit depth as a fraction of the star's flux
            double transitDepthMax = 0.1; // maximum transit depth as a fraction of the star's flux
            int nTransitsMin = 2; // minimum number of transits required for a detection
            double maxPeriod = 5.0; // maximum period to search for in days

            // Load the target star's light curve from a file or database
            double[] times = LoadLightCurveTimes();
            double[] fluxes = LoadLightCurveFluxes();

            // Detrend the light curve using the Biweight midvariance statistic
            double[] detrendedFluxes = TLS.detrend(times, fluxes, method: "biweight");

            // Perform the TLS search
            Dictionary<string, object> tlsResults = TLS.tls_search(times, detrendedFluxes, jwstPos[0], jwstPos[1], maxPeriod, 
                                                                   transit_min_duration_hours: transitMinDuration, 
                                                                   transit_max_duration_hours: transitMaxDuration, 
                                                                   transit_depth_min: transitDepthMin, 
                                                                   transit_depth_max: transitDepthMax, 
                                                                   n_transits_min: nTransitsMin);

            // Print the results
            Console.WriteLine("TLS Results:");
            Console.WriteLine("Number of planet candidates found: {0}", (int)tlsResults["n_planet_candidates"]);
            Console.WriteLine("Periods (days): {0}", string.Join(", ", (double[])tlsResults["periods"]));
            Console.WriteLine("Epochs (BJD): {0}", string.Join(", ", (double[])tlsResults["epochs"]));
            Console.WriteLine("Depths (fraction of flux): {0}", string.Join(", ", (double[])tlsResults["depths"]));
        }
    static double[] LoadLightCurveTimes() //MOCK DATA
    {
    // Define the start and end times of the observation
    DateTime startTime = new DateTime(2023, 2, 18, 0, 0, 0);
    DateTime endTime = new DateTime(2023, 2, 19, 0, 0, 0);

    // Define the time interval between observations
    TimeSpan interval = TimeSpan.FromSeconds(30);

    // Generate a sequence of times spaced by the interval
    List<DateTime> timeList = new List<DateTime>();
    for (DateTime time = startTime; time < endTime; time += interval)
        {
        timeList.Add(time);
        }

    // Convert the sequence of times to an array of Julian dates
    double[] times = new double[timeList.Count];
    for (int i = 0; i < timeList.Count; i++)
        {
        double jd = JulianDate(timeList[i]);
        times[i] = jd;
        }

    return times;
    }
    static double[] LoadLightCurveFluxes()
{
    // Generate mock data for the target star's light curve
    double[] fluxes = new double[1440]; // 1 day of observations, 1 observation every 30 seconds
    Random rand = new Random();
    for (int i = 0; i < fluxes.Length; i++)
    {
        fluxes[i] = 1 + 0.1 * rand.NextDouble(); // add some random noise to the fluxes
    }

    return fluxes;
}

static double JulianDate(DateTime time)
{
    // Calculate the Julian date for a given date and time
    double a = (14 - time.Month) / 12;
    double y = time.Year + 4800 - a;
    double m = time.Month + 12 * a - 3;
    double jd = time.Day + ((153 * m + 2) / 5) + 365 * y + (y / 4) - (y / 100) + (y / 400) - 32045;
    jd += ((time.TimeOfDay.TotalSeconds + (time.Millisecond / 1000.0)) / 86400.0);

    return jd;
}

        // This function retrieves the position of the JWST at a specific time using the Horizons System
        static double[] GetJWSTPosition(DateTime time)
{
            // Set up a WebClient to query the Horizons System
            WebClient client = new WebClient();

            // Construct the query URL
            string url = $"https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='499856'&MAKE_EPHEM='YES'" +
                        $"&TABLE_TYPE='OBSERVER'&START_TIME='{time.ToString("yyyy-MM-dd HH:mm")}'" +
                        $"&STOP_TIME='{time.AddMinutes(1).ToString("yyyy-MM-dd HH:mm")}'&STEP_SIZE='1%20min'&CSV_FORMAT='YES'";

            // Query the Horizons System and obtain the response as a string
            string response = client.DownloadString(url);

            // Parse the response to obtain the JWST's position
            string[] lines = response.Split('\n');
            string[] fields = lines[4].Split(',');
            double ra = double.Parse(fields[3]);
            double dec = double.Parse(fields[4]);
            double distance = double.Parse(fields[5]);

            // Return the JWST position as an array
            return new double[] { ra, dec, distance };
        }
    }
}

