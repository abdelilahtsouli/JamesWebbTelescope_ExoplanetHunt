using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearRegression;

namespace TLS
{
    public class TLS
    {
        static double[] detrend(double[] times, double[] fluxes, string method = "biweight")
        {
            // Copy the input fluxes to avoid modifying the original array
            double[] detrendedFluxes = (double[])fluxes.Clone();

            // Compute the residuals using the specified method
            double[] residuals;
            switch (method.ToLower())
            {
                case "biweight":
                    residuals = ComputeBiweightResiduals(times, detrendedFluxes);
                    break;
                case "linear":
                    residuals = ComputeLinearResiduals(times, detrendedFluxes);
                    break;
                default:
                    throw new ArgumentException("Invalid detrending method specified.");
            }

            // Remove the residuals from the original fluxes
            for (int i = 0; i < fluxes.Length; i++)
            {
                detrendedFluxes[i] -= residuals[i];
            }

            return detrendedFluxes;
        }

        static double[] ComputeBiweightResiduals(double[] times, double[] fluxes, int maxIterations = 10, double c = 6.0)
        {
            // Initialize the residuals to zero
            double[] residuals = new double[fluxes.Length];

            // Compute the initial weights
            double[] weights = ComputeBiweightWeights(residuals, c);

            // Perform the iterative biweight calculation
            for (int iteration = 1; iteration <= maxIterations; iteration++)
            {
                // Compute the weighted location and scale
                double location = WeightedMedian(fluxes, weights);
                double scale = WeightedScale(fluxes, weights, location);

                // Compute the standardized residuals
                for (int i = 0; i < fluxes.Length; i++)
                {
                    residuals[i] = (fluxes[i] - location) / scale;
                }

                // Update the weights
                weights = ComputeBiweightWeights(residuals, c);
            }

            return residuals;
        }

        static double[] ComputeLinearResiduals(double[] times, double[] fluxes)
        {
            // Compute the linear fit using the least-squares method
            double slope; double intercept;
            var ans = SimpleRegression.Fit(times, fluxes);
            slope = ans.Item1;
            intercept = ans.Item2;
            // Compute the predicted fluxes
            double[] predictedFluxes = new double[fluxes.Length];
            for (int i = 0; i < fluxes.Length; i++)
            {
                predictedFluxes[i] = slope * times[i] + intercept;
            }

            // Compute the residuals
            double[] residuals = new double[fluxes.Length];
            for (int i = 0; i < fluxes.Length; i++)
            {
                residuals[i] = fluxes[i] - predictedFluxes[i];
            }

            return residuals;
        }

        static double[] ComputeBiweightWeights(double[] residuals, double c)
        {
            // Compute the median absolute deviation (MAD)
            double mad = MedianAbsoluteDeviation(residuals, c);

            // Initialize the weights to zero
            double[] weights = new double[residuals.Length];

            // Compute the nonzero weights
            for (int i = 0; i < residuals.Length; i++)
            {
                double u = Math.Abs(residuals[i]) / (c * mad);
                if (u <= 1.0)
                {
                    weights[i] = Math.Pow(1.0 - u * u, 2);
                }
            }

            return weights;
        }

        static double WeightedMedian(double[] values, double[] weights)
        {
            // Check that the inputs have the same length
            if (values.Length != weights.Length)
                throw new ArgumentException("The input arrays must have the same length");

            // Sort the values and weights in ascending order
            Array.Sort(values, weights);

            // Compute the sum of the weights
            double sumWeights = weights.Sum();

            // Compute the median value
            double median = 0.0;
            double sum = 0.0;
            for (int i = 0; i < values.Length; i++)
            {
                sum += weights[i];
                if (2 * sum >= sumWeights)
                {
                    median = values[i];
                    break;
                }
            }

            return median;
        }


        static double CalculateMedian(double[] values)
        {
            double[] sortedValues = values.OrderBy(v => v).ToArray();
            int middleIndex = sortedValues.Length / 2;

            if (sortedValues.Length % 2 == 0)
            {
                return (sortedValues[middleIndex - 1] + sortedValues[middleIndex]) / 2.0;
            }
            else
            {
                return sortedValues[middleIndex];
            }
        }
        static double WeightedScale(double[] fluxes, double[] weights, double location)
        {
            int N = fluxes.Length;
            double[] deviations = new double[N];
            for (int i = 0; i < N; i++)
            {
                deviations[i] = Math.Abs(fluxes[i] - location);
            }
            double mad = WeightedMedian(deviations, weights);
            return 1.4826 * mad;
        }


        static double CalculateResidualsVariance(double[] values, double median)
        {
            double[] residuals = new double[values.Length];
            for (int i = 0; i < values.Length; i++)
            {
                residuals[i] = values[i] - median;
            }

            double[] squaredResiduals = residuals.Select(r => r * r).ToArray();
            double variance = squaredResiduals.Sum() / (squaredResiduals.Length - 1);
            return variance;
        }

        static double CalculateBiweight(double[] values, double median, double residualsVariance)
        {
            double MAD = CalculateMAD(values, median);
            double biweightDenominator = 6 * Math.Sqrt(residualsVariance);
            double biweightNumerator = values.Select(v => BiweightWeight(v - median, biweightDenominator / MAD)).Zip(values, (w, v) => w * v).Sum();
            double biweight = biweightNumerator / values.Select(v => BiweightWeight(v - median, biweightDenominator / MAD)).Sum();
            return biweight;
        }
        public static double MedianAbsoluteDeviation(double[] data, double median)
        {
            double[] absoluteDeviations = new double[data.Length];
            for (int i = 0; i < data.Length; i++)
            {
                absoluteDeviations[i] = Math.Abs(data[i] - median);
            }

            double mad = CalculateMedian(absoluteDeviations);
            return mad;
        }


        static double CalculateMAD(double[] values, double median)
        {
            double[] deviations = values.Select(v => Math.Abs(v - median)).ToArray();
            double MAD = CalculateMedian(deviations);
            return MAD;
        }

        static double BiweightWeight(double x, double c)
        {
            if (Math.Abs(x) > c)
            {
                return 0;
            }
            else
            {
                double numerator = 1 - (x * x) / (c * c);
                double denominator = 1 - (x * x) / (c * c);
                return Math.Pow(numerator / denominator, 2);
            }
        }
       

        public static Dictionary<string, object> tls_search(double[] times, double[] fluxes, double ra, double dec,
                                                            double max_period, int transit_min_duration_hours = 1,
                                                            int transit_max_duration_hours = 5,
                                                            double transit_depth_min = 0.001, double transit_depth_max = 0.1,
                                                            int n_transits_min = 2)
        {
            // Detrend the light curve
            double[] detrendedFluxes = detrend(times, fluxes, "biweight");

            // Calculate the transit depths and durations to search over
            double[] transit_durations = Enumerable.Range(transit_min_duration_hours, transit_max_duration_hours - transit_min_duration_hours + 1)
                .Select(hours => hours * 3600.0)
                .ToArray();
            double[] transit_depths = Enumerable.Range(0, 21)
                .Select(i => transit_depth_min + (i / 20.0) * (transit_depth_max - transit_depth_min))
                .ToArray();

            // Define the period grid
            double[] periods = new double[] { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0 };
            int n_bins = (int)Math.Ceiling(max_period / periods[0]);
            periods = Enumerable.Range(0, n_bins).Select(i => (i + 0.5) * periods[0]).Where(p => p <= max_period).ToArray();

            // Calculate the phased light curve for each period
            // Calculate the phased light curve for each period
            Dictionary<double, double[]> phased_light_curves = new Dictionary<double, double[]>();
            foreach (double period in periods)
            {
                double[] phases = times.Select(time => ((time / period) % 1.0 + 1.0) % 1.0)
                                    .Select(phase => phase > 0.5 ? phase - 1.0 : phase).ToArray();
                double[] phased_fluxes = phases.Zip(detrendedFluxes, (phase, flux) => (phase, flux))
                                            .OrderBy(t => t.phase)
                                            .Select(t => t.flux).ToArray();
                phased_light_curves[period] = phased_fluxes;
            }

            // Initialize the best transit dictionary
            Dictionary<string, object> best_transit = new Dictionary<string, object>()
            {
                { "period", 0.0 },
                { "duration", 0.0 },
                { "depth", 0.0 },
                { "transits", new List<Dictionary<string, object>>() },
                { "sde", 0.0 }
            };

            // Perform the search
            foreach (double period in periods)
            {
                double[] phases = times.Select(time => ((time / period) % 1.0 + 1.0) % 1.0)
                                    .Select(phase => phase > 0.5 ? phase - 1.0 : phase).ToArray();

                foreach (double transit_depth in transit_depths)
                {
                    foreach (double transit_duration in transit_durations)
                    {
                        double t0 = 0.5 * transit_duration;

                        double[] in_transit_mask = phases.Select(phase => Math.Abs(phase - (t0 / period)))
                                                        .Select(dt => dt < (0.5 * transit_duration / period) ? 1.0 : 0.0)
                                                        .ToArray();
                    }
                }
            }

            // Calculate the signal-to-noise ratio (SNR) for each transit
            double[] snrs = transitDepths.Zip(transitErrors, (depth, error) => depth / error).ToArray();

            // Identify the transit events that meet the SNR and duration criteria
            List<int> goodTransitIndices = new List<int>();
            for (int i = 0; i < transitDepths.Length; i++)
            {
                if (transitDepths[i] > 0 && transitDepths[i] >= transitDepthMin && transitDepths[i] <= transitDepthMax &&
                    transitDurations[i] >= transitMinDurationHours && transitDurations[i] <= transitMaxDurationHours &&
                    snrs[i] >= snrThreshold && nTransits[i] >= nTransitsMin)
                {
                    goodTransitIndices.Add(i);
                }
            }

            // Initialize arrays to hold the results
            int nCandidates = goodTransitIndices.Count;
            double[] candidate_periods = new double[nCandidates];
            double[] candidate_epochs = new double[nCandidates];
            double[] candidate_depths = new double[nCandidates];

            // Extract the periods, epochs, and depths of the candidate transiting planets
            for (int i = 0; i < nCandidates; i++)
            {
                int transitIndex = goodTransitIndices[i];
                candidate_periods[i] = transitPeriods[transitIndex];
                candidate_epochs[i] = transitEpochs[transitIndex];
                candidate_depths[i] = transitDepths[transitIndex];
            }
        }
    }
}


