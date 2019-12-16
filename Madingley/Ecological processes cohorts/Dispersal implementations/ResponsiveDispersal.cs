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
    public partial class ResponsiveDispersal : IDispersalImplementation
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
        /// Density threshold below which adult individuals may move to look for other adults of the same cohort
        /// </summary>
        /// <remarks>The density scales with cohort weight via: Min Density = DensityThresholdScaling / MatureMass (g)</remarks>
        private const double _DensityThresholdScaling = 50000;
        /// <summary>
        /// Get the density threshold below which adult individuals may move to look for other adults of the same cohort
        /// </summary>
        public double DensityThresholdScaling { get { return _DensityThresholdScaling; } }

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
        /// The proportion of body mass loss at which the cohort will try to disperse every time during a time step
        /// </summary>
        private const double _StarvationDispersalBodyMassThreshold = 0.8;
        /// <summary>
        /// Get the proportion of body mass loss at which the cohort will try to disperse every time during a time step
        /// </summary>
        public double StarvationDispersalBodyMassThreshold { get { return _StarvationDispersalBodyMassThreshold; } }

        /// <summary>
        /// Write out the values of the parameters to an output file
        /// </summary>
        /// <param name="sw">A streamwriter object to write the parameter values to</param>
        public void WriteOutParameterValues(StreamWriter sw)
        {
            // Write out parameters
            sw.WriteLine("Responsive Dispersal\tTimeUnitImplementation\t" + Convert.ToString(_TimeUnitImplementation));
            sw.WriteLine("Responsive Dispersal\t_DispersalSpeedBodyMassScalar\t" + Convert.ToString(_DispersalSpeedBodyMassScalar));
            sw.WriteLine("Responsive Dispersal\t_DispersalSpeedBodyMassExponent\t" + Convert.ToString(_DispersalSpeedBodyMassExponent));
            sw.WriteLine("Responsive Dispersal\tTDensityThresholdScaling\t" + Convert.ToString(_DensityThresholdScaling));
            sw.WriteLine("Responsive Dispersal\tStarvationDispersalBodyMassThreshold\t" + Convert.ToString(_StarvationDispersalBodyMassThreshold));
        }

        private bool CheckStarvationDispersal(ModelGrid gridForDispersal, uint latIndex, uint lonIndex, Cohort cohortToDisperse, int functionalGroup, int cohortNumber)
        {
            // A boolean to check whether a cohort has dispersed
            bool CohortHasDispersed = false;

            // Check for starvation driven dispersal
            // What is the present body mass of the cohort?
            // Note that at present we are just tracking starvation for adults
            double IndividualBodyMass = cohortToDisperse.IndividualBodyMass;
            double AdultMass = cohortToDisperse.AdultMass;

            // Assume a linear relationship between probability of dispersal and body mass loss, up to _StarvationDispersalBodyMassThreshold
            // at which point the cohort will try to disperse every time step
            if (IndividualBodyMass < AdultMass)
            {
                double ProportionalPresentMass = IndividualBodyMass / AdultMass;

                // If the body mass loss is greater than the starvation dispersal body mass threshold, then the cohort tries to disperse
                if (ProportionalPresentMass < _StarvationDispersalBodyMassThreshold)
                {
                    // Cohort tries to disperse
                    double[] DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(AdultMass));
                    double CohortDispersed = CheckForDispersal(DispersalArray[0]);
                    if (CohortDispersed > 0)
                    {
                        uint[] DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, CohortDispersed, DispersalArray[4], DispersalArray[5]);
                        
                        // Update the delta array of cells to disperse to, if the cohort moves
                        if (DestinationCell[0] < 999999)
                        {
                            // Update the delta array of cohorts
                            gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex, lonIndex].Add((uint)functionalGroup);
                            gridForDispersal.DeltaCohortNumberDispersalArray[latIndex, lonIndex].Add((uint)cohortNumber);
                        
                            // Update the delta array of cells to disperse to
                            gridForDispersal.DeltaCellToDisperseToArray[latIndex, lonIndex].Add(DestinationCell);
                        }
                    }

                    // Note that regardless of whether or not it succeeds, if a cohort tries to disperse, it is counted as having dispersed for the purposes of not then allowing it to disperse
                    // based on its density.
                    CohortHasDispersed = true;
                }

                // Otherwise, the cohort has a chance of trying to disperse proportional to its mass lass
                else
                {
                    // Cohort tries to disperse with a particular probability
                    // Draw a random number
                    if (((1.0 - ProportionalPresentMass) / (1.0 - _StarvationDispersalBodyMassThreshold)) > RandomNumberGenerator.GetUniform())
                    {
                        // Cohort tries to disperse
                        double[] DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(AdultMass));
                        double CohortDispersed = CheckForDispersal(DispersalArray[0]);
                        if (CohortDispersed > 0)
                        {
                            uint[] DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, CohortDispersed, DispersalArray[4], DispersalArray[5]);

                            // Update the delta array of cells to disperse to, if the cohort moves
                            if (DestinationCell[0] < 999999)
                            {
                                // Update the delta array of cohorts
                                gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex, lonIndex].Add((uint)functionalGroup);
                                gridForDispersal.DeltaCohortNumberDispersalArray[latIndex, lonIndex].Add((uint)cohortNumber);

                                // Update the delta array of cells to disperse to
                                gridForDispersal.DeltaCellToDisperseToArray[latIndex, lonIndex].Add(DestinationCell);
                            }
                        }

                        CohortHasDispersed = true;
                    }
                }

            }
            return CohortHasDispersed;
        }

        private void CheckDensityDrivenDispersal(ModelGrid gridForDispersal, uint latIndex, uint lonIndex, Cohort cohortToDisperse, int functionalGroup, int cohortNumber)
        {
            // Check the population density
            double NumberOfIndividuals = cohortToDisperse.CohortAbundance;

            // Get the cell area, in kilometres squared
            double CellArea = gridForDispersal.GetCellEnvironment(latIndex, lonIndex)["Cell Area"][0];

            // If below the density threshold
            if ((NumberOfIndividuals / CellArea) < _DensityThresholdScaling / cohortToDisperse.AdultMass)
            {
                // Check to see if it disperses (based on the same movement scaling as used in diffusive movement)
                // Calculate dispersal speed for that cohort
                double DispersalSpeed = CalculateDispersalSpeed(cohortToDisperse.IndividualBodyMass);

                // Cohort tries to disperse
                double[] DispersalArray = CalculateDispersalProbability(gridForDispersal, latIndex, lonIndex, CalculateDispersalSpeed(cohortToDisperse.AdultMass));
                
                double CohortDispersed = CheckForDispersal(DispersalArray[0]);
                
                if (CohortDispersed > 0)
                {
                    uint[] DestinationCell = CellToDisperseTo(gridForDispersal, latIndex, lonIndex, DispersalArray, CohortDispersed, DispersalArray[4], DispersalArray[5]);

                    // Update the delta array of cells to disperse to, if the cohort moves
                    if (DestinationCell[0] < 999999)
                    {
                        // Update the delta array of cohorts
                        gridForDispersal.DeltaFunctionalGroupDispersalArray[latIndex, lonIndex].Add((uint)functionalGroup);
                        gridForDispersal.DeltaCohortNumberDispersalArray[latIndex, lonIndex].Add((uint)cohortNumber);

                        // Update the delta array of cells to disperse to
                        gridForDispersal.DeltaCellToDisperseToArray[latIndex, lonIndex].Add(DestinationCell);
                    }
                }
            }
        }

        /// <summary>
        /// Calculate the average diffusive dispersal speed of individuals in a cohort given their body mass
        /// </summary>
        /// <param name="bodyMass">The current body mass of an individual in the cohort</param>
        /// <returns>The (average) dispersal speed in kilometres per month</returns>
        private double CalculateDispersalSpeed(double bodyMass)
        {
            return _DispersalSpeedBodyMassScalar * Math.Pow(bodyMass, _DispersalSpeedBodyMassExponent);
        }

        /// <summary>
        /// Calculates the probability of responsive dispersal given average individual dispersal speed and grid cell
        /// </summary>
        /// <param name="madingleyGrid">The model grid</param>
        /// <param name="latIndex">The latitude index of the grid cell to check for dispersal</param>
        /// <param name="lonIndex">The longitude index of the grid cell to check for dispersal</param>
        /// <param name="dispersalSpeed">The average dispersal speed of individuals in the acting cohort</param>
        /// <returns>A six element array. 
        /// The first element is the probability of dispersal.
        /// The second element is the probability of dispersing in the u (longitudinal) direction
        /// The third element is the probability of dispersing in the v (latitudinal) direction
        /// The fourth element is the probability of dispersing in the diagonal direction
        /// The fifth element is the u velocity
        /// The sixth element is the v velocity
        /// Note that the second, third, and fourth elements are always positive; thus, they do not indicate 'direction' in terms of dispersal.</returns>
        private double[] CalculateDispersalProbability(ModelGrid madingleyGrid, uint latIndex, uint lonIndex, double dispersalSpeed)
        {
            double LatCellLength = madingleyGrid.CellHeightsKm[latIndex];
            double LonCellLength = madingleyGrid.CellWidthsKm[latIndex];

            // Pick a direction at random
            double RandomDirection = RandomNumberGenerator.GetUniform() * 2 * Math.PI;

            // Calculate the u and v components given the dispersal speed
            double uSpeed = dispersalSpeed * Math.Cos(RandomDirection);
            double vSpeed = dispersalSpeed * Math.Sin(RandomDirection);

            // Check that the whole cell hasn't moved out (i.e. that dispersal speed is not greater than cell length). 
            // This could happen if dispersal speed was high enough; indicates a need to adjust the time step, or to slow dispersal
            if ((uSpeed > LonCellLength) || (vSpeed > LatCellLength))
            {
                Debug.Fail("Dispersal probability should always be <= 1");
            }

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

            // Check that we don't have any issues
            if (DispersalProbability >= 1)
            {
                //Debug.Fail("Dispersal probability should always be <= 1");
                DispersalProbability = 1.0;
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
        // Here we make the assumption that if the cell in the direction chosen is unsuitable, that the dispersal does not happen 
        // (i.e. that the cohorts 'bumps up' against unsuitable habitat
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
                    // Diagonally. Also allow for rounding errors here, otherwise we can get sent to the debug.fail incorrectly
                    // Only an issue here, because otherwise the random value is always larger than the sum of the elements of the dispersal array
                    if (RandomValue <= (dispersalArray[1] + dispersalArray[2] + dispersalArray[3]) + 0.0000000001)
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
                        Debug.Fail("Error selecting cell to disperse to in Responsivedispersal.cs");
                        DestinationCell = new uint[2] { 9999999, 9999999 };
                    }
                }

            }
            return DestinationCell;
        }

    }
}
