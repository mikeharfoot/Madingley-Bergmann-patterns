using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;
using System.IO;

namespace Madingley
{
    /// <summary>
    /// A formulation of the process of dispersal
    /// </summary>
    public partial class DiffusiveDispersal: IDispersalImplementation
    {
        /// <summary>
        /// The time units associated with this implementation of dispersal
        /// </summary>
        private const string _TimeUnitImplementation = "month";
        /// <summary>
        /// Get the time units associated with this implementation of dispersal
        /// </summary>
        public string TimeUnitImplementation { get { return _TimeUnitImplementation; } }

        /// <summary>
        /// Scalar relating dispersal speed to individual body mass
        /// </summary>
        private const double _DispersalSpeedBodyMassScalar = 0.0278;
        /// <summary>
        /// Get the scalar relating dispersal speed to individual body mass
        /// </summary>
        public double DispersalSpeedBodyMassScalar { get { return _DispersalSpeedBodyMassScalar; } }

        /// <summary>
        /// Body-mass exponent of the relationship between disperal speed and individual body mass
        /// </summary>
        private const double _DispersalSpeedBodyMassExponent = 0.48;
        /// <summary>
        /// Get the body-mass exponent of the relationship between disperal speed and individual body mass
        /// </summary>
        public double DispersalSpeedBodyMassExponent { get { return _DispersalSpeedBodyMassExponent; } }

        /// <summary>
        /// Write out the values of the parameters to an output file
        /// </summary>
        /// <param name="sw">A streamwriter object to write the parameter values to</param>
        public void WriteOutParameterValues(StreamWriter sw)
        {
            // Write out parameters
            sw.WriteLine("Diffusive Dispersal\tTimeUnitImplementation\t" + Convert.ToString(_TimeUnitImplementation));
            sw.WriteLine("Diffusive Dispersal\tDispersalSpeedBodyMassScalar_per_g\t" + Convert.ToString(_DispersalSpeedBodyMassScalar));
            sw.WriteLine("Diffusive Dispersal\tDispersalSpeedBodyMassExponent\t" + Convert.ToString(_DispersalSpeedBodyMassExponent));
        }
        
        /// <summary>
        /// Calculates the average diffusive dispersal speed of individuals in a cohort given their body mass
        /// </summary>
        /// <param name="bodyMass">The current body mass of individuals in the cohort</param>
        /// <returns>The average dispersal speed, in km per month</returns>
        private double CalculateDispersalSpeed(double bodyMass)
        {
                return _DispersalSpeedBodyMassScalar * Math.Pow(bodyMass,_DispersalSpeedBodyMassExponent);
        }

        /// <summary>
        /// Calculates the probability of diffusive dispersal given average individual dispersal speed
        /// </summary>
        /// <param name="madingleyGrid">The model grid</param>
        /// <param name="latIndex">The latitude index of the grid cell to check for dispersal</param>
        /// <param name="lonIndex">The longitude index of the grid cell to check for dispersal</param>
        /// <param name="dispersalSpeed">The average speed at which individuals in this cohort move around their environment, in km per month</param>
        /// <returns>A six element array. 
        /// The first element is the probability of dispersal.
        /// The second element is the probability of dispersing in the u (longitudinal) direction
        /// The third element is the probability of dispersing in the v (latitudinal) direction
        /// The fourth element is the probability of dispersing in the diagonal direction
        /// The fifth element is the u velocity modified by the random diffusion component
        /// The sixth element is the v velocity modified by the random diffusion component
        /// Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.</returns>
        private double[] CalculateDispersalProbability(ModelGrid madingleyGrid, uint latIndex, uint lonIndex, double dispersalSpeed)
        {
            // Check that the u speed and v speed are not greater than the cell length. If they are, then rescale them; this limits the max velocity
            // so that cohorts cannot be advected more than one grid cell per time step
            double LatCellLength = madingleyGrid.CellHeightsKm[latIndex];
            double LonCellLength = madingleyGrid.CellWidthsKm[latIndex];

            // Pick a direction at random
            double RandomDirection = RandomNumberGenerator.GetUniform() * 2 * Math.PI;

            // Calculate the u and v components given the dispersal speed
            double uSpeed = dispersalSpeed * Math.Cos(RandomDirection);
            double vSpeed = dispersalSpeed * Math.Sin(RandomDirection);
            
            // Calculate the area of the grid cell that is now outside in the diagonal direction
            double AreaOutsideBoth = Math.Abs(uSpeed * vSpeed);

            // Calculate the area of the grid cell that is now outside in the u direction (not including the diagonal)
            double AreaOutsideU = Math.Abs(uSpeed * LatCellLength) - AreaOutsideBoth;

            // Calculate the proportion of the grid cell that is outside in the v direction (not including the diagonal
            double AreaOutsideV = Math.Abs(vSpeed * LonCellLength) - AreaOutsideBoth;

            // Get the cell area, in kilometres squared
            double CellArea = madingleyGrid.GetCellEnvironment(latIndex, lonIndex)["Cell Area"][0];

            // Convert areas to a probability
            double DispersalProbability = (AreaOutsideU + AreaOutsideV + AreaOutsideBoth) / CellArea;

            // Check that the whole cell hasn't moved out. This could happen if dispersal speed was high enough
            if (DispersalProbability >= 1)
            {
                Debug.Fail("Dispersal probability in diffusion should always be <= 1");
            }

            double[] NewArray = { DispersalProbability, AreaOutsideU / CellArea, AreaOutsideV / CellArea, AreaOutsideBoth / CellArea, uSpeed, vSpeed };

            return NewArray;
        }

        /// <summary>
        /// Generate a random value to see if a cohort disperses
        /// </summary>
        /// <param name="dispersalProbability">The probability of dispersal</param>
        /// <returns>Returns either the random value, if it less than dispersal probability, or -1</returns>
        private double CheckForDispersal(double dispersalProbability)
        {
            // Randomly check to see if dispersal occurs
            double RandomValue = RandomNumberGenerator.GetUniform();
            if (dispersalProbability >= RandomValue)
            {
                return RandomValue;
            }
            else
            {
                return -1.0;
            }
        }

        // Determine to which cell the cohort disperses
        // Note that if the direction is not dispersable, then it doesn't happen
        private uint[] CellToDisperseTo(ModelGrid madingleyGrid, uint latIndex, uint lonIndex, double[] dispersalArray, double RandomValue, double uSpeedIncDiffusion, double vSpeedIncDiffusion)
        {
            uint[] DestinationCell;

            // Check to see in which axis the cohort disperses
            // Longitudinally
            if (RandomValue <= dispersalArray[1])
            {
                // Work out whether dispersal is to the cell to the E or the W
                if (uSpeedIncDiffusion > 0)
                {
                    DestinationCell = madingleyGrid.CheckDispersalEast(latIndex, lonIndex);
                }
                else
                {
                    DestinationCell = madingleyGrid.CheckDispersalWest(latIndex, lonIndex);
                }

            }
            else
            {
                // Latitudinally
                if (RandomValue <= (dispersalArray[1] + dispersalArray[2]))
                {
                    // Work out whether dispersal is to the cell to the N or the S
                    if (vSpeedIncDiffusion > 0)
                    {
                        DestinationCell = madingleyGrid.CheckDispersalNorth(latIndex, lonIndex);
                    }
                    else
                    {
                        DestinationCell = madingleyGrid.CheckDispersalSouth(latIndex, lonIndex);
                    }

                }
                else
                {
                    // Diagonally
                    if (RandomValue <= (dispersalArray[1] + dispersalArray[2] + dispersalArray[3]))
                    {
                        // Work out to which cell dispersal occurs
                        if (uSpeedIncDiffusion > 0)
                        {
                            if (vSpeedIncDiffusion > 0)
                            {
                                DestinationCell = madingleyGrid.CheckDispersalNorthEast(latIndex, lonIndex);
                            }
                            else
                            {
                                DestinationCell = madingleyGrid.CheckDispersalSouthEast(latIndex, lonIndex);
                            }

                        }
                        else
                        {
                            if (vSpeedIncDiffusion > 0)
                            {
                                DestinationCell = madingleyGrid.CheckDispersalNorthWest(latIndex, lonIndex);
                            }
                            else
                            {
                                DestinationCell = madingleyGrid.CheckDispersalSouthWest(latIndex, lonIndex);
                            }
                        }
                    }
                    else
                    {
                        //INSERT ERROR HANDLING CODE HERE("Error with advection when determining cell to disperse to");
                        DestinationCell = new uint[2] { 9999999, 9999999 };
                    }
                }

            }
            return DestinationCell;
        }
    }
}
