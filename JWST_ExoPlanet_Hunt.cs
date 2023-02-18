using System;
using System.Collections.Generic;
using System.Net;
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

        // This function retrieves the position of the JWST at a specific time using the Horizons System
        static double[] GetJWSTPosition(DateTime time)
        {
            // Set up a WebClient to query the Horizons System
            WebClient client = new WebClient();

            // Construct the query URL
            string url = "https://ssd.jpl.nasa.gov/horizons_batch.cgi?batch=1&COMMAND='499856'&MAKE_EPHEM='YES'" +
                         "&TABLE_TYPE='OBSERVER'&START_TIME='2023-02-18%2000:00'" +
                         "&STOP_TIME='2023-02-18%2000:01'&STEP_SIZE='1%20min'&CSV_FORMAT='YES'";

            // Query the Horizons System and obtain the response as a string
            string response = client.DownloadString(url);

            // Parse the response to obtain the JWST's position
            string[] lines = response.Split('\n');
            string[] fields = lines[4].Split(',');
            double ra = double.Parse(fields[3]);
            double dec = double.Parse(fields[4]);
            double distance = double.Parse(fields[5]);

            // Print the results
            Console.WriteLine("JWST Position:");
            Console.WriteLine("RA: {0} deg", ra);
            Console.WriteLine("Dec: {0} deg", dec);
            Console.WriteLine("Distance: {0} km", distance);
        }
    }
}

