using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// Removes individuals from animal cohorts to simulate the effects of direct harvesting
    /// </summary>
    public class Harvesting
    {
        /// <summary>
        /// Constructor for harvesting class
        /// </summary>
        public Harvesting()
        {
        }

        /// <summary>
        /// Remove individuals lost from cohorts through direct harvesting of animals
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the current grid cell</param>
        /// <param name="harvestingScenario">The scenario of direct harvesting of animals to apply</param>
        /// <param name="currentTimestep">The current model time step</param>
        /// <param name="burninSteps">The number of time steps to spin the model up for before applying the harvesting scenario</param>
        /// <param name="impactSteps">The number of time steps to apply the scenario for</param>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        /// <param name="impactCell">The index of the cell, within the list of all cells to run, to apply the scenario for</param>
        public void RemoveHarvestedIndividuals(GridCellCohortHandler gridCellCohorts,
            Tuple<string, double, double> harvestingScenario, uint currentTimestep, uint burninSteps, uint impactSteps, uint totalSteps,
            SortedList<string,double[]> cellEnvironment, Boolean impactCell)
        {
            if (impactCell)
            {
                if (harvestingScenario.Item1 == "no")
                {
                    // Do not apply any harvesting
                }
                else if (harvestingScenario.Item1 == "constant")
                {
                    // If the burn-in period has been completed, then apply
                    // the harvesting scenario
                    if (currentTimestep > burninSteps)
                    {
                        ApplyHarvesting(gridCellCohorts, harvestingScenario.Item2, cellEnvironment);
                    }
                }
                else if (harvestingScenario.Item1 == "temporary")
                {
                    // If the burn-in period has been completed and the period of impact has not elapsed,
                    // then apply the harvesting scenario
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        ApplyHarvesting(gridCellCohorts, harvestingScenario.Item2, cellEnvironment);
                    }
                }
                else if (harvestingScenario.Item1 == "escalating")
                {
                    // If the spin-up period has been completed, then apply a level of harvesting
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if (currentTimestep > burninSteps)
                    {
                        // Calculate the target biomass for harvesting based on the number of time steps that have elapsed since the spin-up
                        double TargetBiomass = (Math.Min(50000, (((currentTimestep - burninSteps) / 12.0) * harvestingScenario.Item2)));

                        // Apply the harvesting scenario using the calculated target biomass
                        ApplyHarvesting(gridCellCohorts, TargetBiomass, cellEnvironment);
                    }

                }
                else if (harvestingScenario.Item1 == "temp-escalating-declining")
                {
                    // If the spin-up period has been completed, then apply a level of harvesting
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        // Calculate the target biomass for harvesting based on the number of time steps that have elapsed since the spin-up
                        double TargetBiomass = (Math.Min(50000, (((currentTimestep - burninSteps) / 12.0) * harvestingScenario.Item2)));

                        // Apply the harvesting scenario using the calculated target biomass
                        ApplyHarvesting(gridCellCohorts, TargetBiomass, cellEnvironment);
                    }
                    else
                    {
                        // Calculate the target biomass for harvesting based on the number of time steps that have elapsed since the spin-up
                        double TargetBiomass = (Math.Min(50000, ((-(totalSteps - currentTimestep) / 12.0) * harvestingScenario.Item2)));

                        // Apply the harvesting scenario using the calculated target biomass
                        ApplyHarvesting(gridCellCohorts, TargetBiomass, cellEnvironment);
                    }

                }
                else if (harvestingScenario.Item1 == "tempescalating")
                {
                    // If the spin-up period has been completed and the period of impact has not elapsed, 
                    // then remove a proportion of plant matter
                    // according to the number of time-steps that have elapsed since the spin-up ended
                    if ((currentTimestep > burninSteps) && (currentTimestep <= (burninSteps + impactSteps)))
                    {
                        // Calculate the target biomass for harvesting based on the number of time steps that have elapsed since the spin-up
                        double TargetBiomass = (Math.Min(50000, (((currentTimestep - burninSteps) / 12.0) * harvestingScenario.Item2)));

                        // Apply the harvesting scenarion using the calculated target biomass
                        ApplyHarvesting(gridCellCohorts, TargetBiomass, cellEnvironment);
                    }
                }
                else
                {
                    Debug.Fail("There is no method for the harvesting scenario specified");
                }

            }
        }

        /// <summary>
        /// Apply the results of direct harvesting of animals to the cohorts in a grid cell
        /// </summary>
        /// <param name="gridCellCohorts">The cohorts in the current grid cell</param>
        /// <param name="targetBiomass">The target biomass to be harvested</param>
        /// <param name="cellEnvironment">The environment in the current grid cell</param>
        public void ApplyHarvesting(GridCellCohortHandler gridCellCohorts, double targetBiomass,
            SortedList<string, double[]> cellEnvironment)
        {
            // Create variable to hold total available animal biomass
            double TotalAvailableBiomass = 0.0;

            // Create variable to hold estimate of preference for a given cohort
            double Preference;

            // Create variable to hold the biomass of a cohort actually harvested
            double BiomassHarvested;

            // Create jagged arrays mirroring the cohort handler to hold calculations
            double[][] AvailableBiomass = new double[gridCellCohorts.Count][];
            double[][] BiomassTimesPreference = new double[gridCellCohorts.Count][];

            // Loop over functional groups and initialise rows in the jagged arrays
            for (int i = 0; i < gridCellCohorts.Count; i++)
            {
                AvailableBiomass[i] = new double[gridCellCohorts[i].Count];
                BiomassTimesPreference[i] = new double[gridCellCohorts[i].Count];
            }

            // Convert target biomass from kg per km squared to g per cell
            targetBiomass *= 1000; //kg to g conversion
            targetBiomass *= cellEnvironment["Cell Area"][0];

            // Loop over cohorts and calculate available biomass and biomass times preference in each cohort, and total available biomass
            for (int fg = 0; fg < gridCellCohorts.Count; fg++)
            {
                for (int c = 0; c < gridCellCohorts[fg].Count; c++)
                {
                    TotalAvailableBiomass += (gridCellCohorts[fg][c].IndividualBodyMass * gridCellCohorts[fg][c].CohortAbundance);
                    AvailableBiomass[fg][c] = (gridCellCohorts[fg][c].IndividualBodyMass * gridCellCohorts[fg][c].CohortAbundance);
                    Preference = 1 / (1 + Math.Exp(-(-8 + 0.8 * Math.Log(gridCellCohorts[fg][c].IndividualBodyMass))));
                    BiomassTimesPreference[fg][c] = AvailableBiomass[fg][c] * Preference;
                }
            }

            // Loop over cohorts again, and calculate and apply the actual amount of biomass harvested
            for (int fg = 0; fg < gridCellCohorts.Count; fg++)
            {
                for (int c = 0; c < gridCellCohorts[fg].Count; c++)
                {
                    BiomassHarvested = Math.Min(AvailableBiomass[fg][c], targetBiomass * BiomassTimesPreference[fg][c] / TotalAvailableBiomass);
                    gridCellCohorts[fg][c].CohortAbundance -= (BiomassHarvested / gridCellCohorts[fg][c].IndividualBodyMass);
                }
            }

        }

    }
}
