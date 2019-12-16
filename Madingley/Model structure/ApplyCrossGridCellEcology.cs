using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace Madingley
{
    /// <summary>
    /// Class for applying changes from the cross-grid cell ecological processes. These are held in matrices of lists in the modelgrid structure.
    /// We simply loop through each cell, and check to see if there are any cohorts flagged as needing to be dispersed. If so, we point the cohort list
    /// in the new grid cell to this cohort, we delete the pointer to it and to the new grid cell in the model grid delta structures, and we delete the pointer 
    /// to it in the original cell
    /// 
    /// We can also output diagnostics here (temporarily) as the whole grid needs to be completed before dispersal is enacted.
    /// </summary>
    public class ApplyCrossGridCellEcology
    {
        /// <summary>
        /// Apply all updates from the ecological processes to the properties of the acting cohort and to the environment
        /// </summary>
        public void UpdateAllCrossGridCellEcology(ModelGrid madingleyModelGrid, ref uint dispersalCounter, CrossCellProcessTracker trackCrossCellProcesses, uint currentTimeStep)
        {
                // Create an array to hold the number of cohorts dispersing in each direction from each grid cell
                uint[, ,] InboundCohorts = new uint[madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0), madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1), 8];

                // Create an array to hold the number of cohorts dispersing in each direction to each grid cell
                uint[, ,] OutboundCohorts = new uint[madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0), madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1), 8];

                // Create an list array to hold the weights of cohorts dispersing from grid cell. Dimensions are: num grid cells lon, num grid cells lat, num cohorts dispersing
                List<double>[,] OutboundCohortWeights = new List<double>[madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0), madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1)];

                for (uint ii = 0; ii < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0); ii++)
                {
                    for (uint jj = 0; jj < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1); jj++)
                    {
                        OutboundCohortWeights[ii,jj] = new List<double>();
                    }
                }


            // Loop through the delta array that holds the grid cells of the cohorts that are flagged as needing to be moved
            for (uint ii = 0; ii < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0) ; ii++)
            {
                for (uint jj = 0; jj < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1); jj++)
                {
                    if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj] != null)
                    {
                        // No cohorts to move if there are none in the delta dispersal array
                        if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj].Count == 0)
                        {
                        }
                        // Otherwise, loop through the cohorts and change the pointers/references to them one-by-one
                        else
                        {
                            for (int kk = 0; kk < madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj].Count; kk++)
                            {
                                // Find out which grid cell it is going to
                                uint[] CellToDisperseTo = madingleyModelGrid.DeltaCellToDisperseToArray[ii, jj].ElementAt(kk);

                                // Functional group is identified by the first array
                                uint CohortToDisperseFG = madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj].ElementAt(kk);

                                // Cohort number is identified by the second array
                                uint CohortToDisperseNum = madingleyModelGrid.DeltaCohortNumberDispersalArray[ii, jj].ElementAt(kk);

                                // Simmply add it to the existing cohorts in that FG in the grid cell to disperse to
                                madingleyModelGrid.AddNewCohortToGridCell(CellToDisperseTo[0], CellToDisperseTo[1], (int)CohortToDisperseFG, madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum));

                                // Update the dispersal counter
                                dispersalCounter++;

                                // If track processes is on, add this to the list of inbound and outbound cohorts
                                if (trackCrossCellProcesses.TrackCrossCellProcesses)
                                {
                                    uint[] TemporaryCell;
                                    Boolean FoundCell = false;

                                    // Calculate in which direction the cohort is exiting and add to the outbound and inbound arrays

                                    // Leaving north border
                                    TemporaryCell = madingleyModelGrid.CheckDispersalNorth(ii, jj);
                                    if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                    {
                                        // Record the cohort as leaving the outbound cell
                                        OutboundCohorts[ii, jj, 0]++;

                                        // Record the cohort as incoming (from the S border) in the inbound cell
                                        InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 4]++;

                                        // Add the weight of the outbound cohort to the weights array
                                        OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                        FoundCell = true;
                                    }

                                    if (!FoundCell)
                                    {
                                        // Leaving north-east border
                                        TemporaryCell = madingleyModelGrid.CheckDispersalNorthEast(ii, jj);
                                        if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                        {
                                            OutboundCohorts[ii, jj, 1]++;
                                            InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 5]++;

                                            // Add the weight of the outbound cohort to the weights array
                                            OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                            FoundCell = true;
                                        }

                                        if (!FoundCell)
                                        {
                                            // Leaving east border
                                            TemporaryCell = madingleyModelGrid.CheckDispersalEast(ii, jj);
                                            if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                            {
                                                OutboundCohorts[ii, jj, 2]++;
                                                InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 6]++;

                                                // Add the weight of the outbound cohort to the weights array
                                                OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                FoundCell = true;
                                            }

                                            if (!FoundCell)
                                            {
                                                // Leaving south-east border
                                                TemporaryCell = madingleyModelGrid.CheckDispersalSouthEast(ii, jj);
                                                if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                                {
                                                    OutboundCohorts[ii, jj, 3]++;
                                                    InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 7]++;

                                                    // Add the weight of the outbound cohort to the weights array
                                                    OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                    FoundCell = true;
                                                }
                                                if (!FoundCell)
                                                {
                                                    // Leaving south border
                                                    TemporaryCell = madingleyModelGrid.CheckDispersalSouth(ii, jj);
                                                    if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                                    {
                                                        OutboundCohorts[ii, jj, 4]++;
                                                        InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 0]++;

                                                        // Add the weight of the outbound cohort to the weights array
                                                        OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                        FoundCell = true;
                                                    }
                                                    if (!FoundCell)
                                                    {
                                                        // Leaving south-west border
                                                        TemporaryCell = madingleyModelGrid.CheckDispersalSouthWest(ii, jj);
                                                        if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                                        {
                                                            OutboundCohorts[ii, jj, 5]++;
                                                            InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 1]++;

                                                            // Add the weight of the outbound cohort to the weights array
                                                            OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                            FoundCell = true;
                                                        }
                                                        if (!FoundCell)
                                                        {
                                                            // Leaving west border
                                                            TemporaryCell = madingleyModelGrid.CheckDispersalWest(ii, jj);
                                                            if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                                            {
                                                                OutboundCohorts[ii, jj, 6]++;
                                                                InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 2]++;

                                                                // Add the weight of the outbound cohort to the weights array
                                                                OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                                FoundCell = true;
                                                            }
                                                            if (!FoundCell)
                                                            {
                                                                // Leaving north-west border
                                                                TemporaryCell = madingleyModelGrid.CheckDispersalNorthWest(ii, jj);
                                                                if ((TemporaryCell[0] == CellToDisperseTo[0]) && (TemporaryCell[1] == CellToDisperseTo[1]))
                                                                {
                                                                    OutboundCohorts[ii, jj, 7]++;
                                                                    InboundCohorts[CellToDisperseTo[0], CellToDisperseTo[1], 3]++;

                                                                    // Add the weight of the outbound cohort to the weights array
                                                                    OutboundCohortWeights[ii, jj].Add(madingleyModelGrid.GetGridCellIndividualCohort(ii, jj, (int)CohortToDisperseFG, (int)CohortToDisperseNum).IndividualBodyMass);

                                                                    FoundCell = true;
                                                                }
                                                                if (!FoundCell)
                                                                {
                                                                    Console.WriteLine("Error: did not identify inbound/outbound cell when filling in dispersal process tracker");
                                                                    break;
                                                                }
                                                            }

                                                        }

                                                    }

                                                }

                                            }
                                        }
                                    }

                                }

                                // So now there is a pointer in the grid cell to which it is going. We have to delete the pointers in the original cell and in the
                                // delta array, but we need to do this without messing with the list structure; i.e. wait until all cohorts have been moved
                            }
                        }
                    }
                }
            }


            // Reset the delta array and remove the pointer to the cohort in the original list
            for (uint ii = 0; ii < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(0); ii++)
            {
                for (uint jj = 0; jj < madingleyModelGrid.DeltaFunctionalGroupDispersalArray.GetLength(1); jj++)
                {
                    if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj] != null)
                    {
                        // No cohorts to move if there are none in the delta dispersal array
                        if (madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj].Count == 0)
                        {
                        }
                        // Otherwise, loop through the cohorts and change the pointers/references to them one-by-one
                        else
                        {
                            // Delete the cohorts from the original grid cell. Note that this needs to be done carefully to ensure that the correct ones 
                            // are deleted (lists shift about when an internal element is deleted.
                            madingleyModelGrid.DeleteGridCellIndividualCohorts(ii, jj, madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj], madingleyModelGrid.DeltaCohortNumberDispersalArray[ii, jj]);

                            // Reset the lists in the delta dispersal arrays
                            madingleyModelGrid.DeltaFunctionalGroupDispersalArray[ii, jj] = new List<uint>();
                            madingleyModelGrid.DeltaCohortNumberDispersalArray[ii, jj] = new List<uint>();

                            // Reset the list in the grid cells to disperse to array
                            madingleyModelGrid.DeltaCellToDisperseToArray[ii, jj] = new List<uint[]>();
                        }
                    }
                }
            }

            if (trackCrossCellProcesses.TrackCrossCellProcesses)
            {
                // If we are tracking dispersal, then write out how many cohorts have moved to a file
                trackCrossCellProcesses.RecordDispersalForACell(InboundCohorts, OutboundCohorts, OutboundCohortWeights, currentTimeStep, madingleyModelGrid);
            }
        }

    }
}
