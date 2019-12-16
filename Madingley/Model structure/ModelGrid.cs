using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.CSV;
using Microsoft.Research.Science.Data.Imperative;
using System.Diagnostics;

using Timing;

using System.Threading;
using System.Threading.Tasks;


namespace Madingley
{
    /// <summary>
    /// A class containing the model grid (composed of individual grid cells) along with grid attributes.
    /// The model grid is referenced by [Lat index, Lon index]\
    /// <todoD>Check Set and Get state variable methods</todoD>
    /// </summary>
    public class ModelGrid
    {
        // Private variable to make sure that not more than one grid is instantiated
        private uint NumGrids = 0;

        // Model grid standardised missing value (applied to all grid cells)
        private double _GlobalMissingValue = -9999;

        /// <summary>
        /// Get the global missing value
        /// </summary>
        public double GlobalMissingValue
        {
            get { return _GlobalMissingValue; }
        }


        // Field to hold minimum latitude of the grid
        private float _MinLatitude;    
        /// <summary>
        /// Get the lower latitude of the lowest cell of the grid
        /// </summary>
        public float MinLatitude
        {
            get { return _MinLatitude; }
        }

        // Field to hold minumum longitude of the grid
        private float _MinLongitude;
        /// <summary>
        /// Get the leftmost longitude of the leftmost cell of the grid
        /// </summary>
        public float MinLongitude
        {
            get { return _MinLongitude; }
        }
        
        // Field to hold maximum latitude of the grid
        private float _MaxLatitude;
        /// <summary>
        /// Get the lowest latitude of the highest cell in the grid
        /// </summary>
        public float MaxLatitude
        {
            get { return _MaxLatitude; }
        }
        
        // Field to hold maximum longitude of the grid
        private float _MaxLongitude;
        /// <summary>
        /// Get the leftmost longitude of the rightmost cell in the grid
        /// </summary>
        public float MaxLongitude
        {
            get { return _MaxLongitude; }
        }
        
        // Field to hold latitude resolution of each grid cell
        private float _LatCellSize;
        /// <summary>
        /// Get the latitudinal length of each grid cell. Currently assumes all cells are equal sized.
        /// </summary>
        public float LatCellSize
        {
            get { return _LatCellSize; }
        }
        
        // Field to hold longitude resolution of each grid cell
        private float _LonCellSize;
        /// <summary>
        /// Get the longitudinal length of each grid cell. Currently assumes all cells are equal sized. 
        /// </summary>
        public float LonCellSize
        {
            get { return _LonCellSize; }
        }

        /// <summary>
        /// The rarefaction of grid cells to be applied to active cells in the model grid
        /// </summary>
        private int _GridCellRarefaction;
        /// <summary>
        /// Get the rarefaction of grid cells to be applied to active cells in the model grid
        /// </summary>
        public int GridCellRarefaction
        { get { return _GridCellRarefaction; } }
        
        /// <summary>
        /// The number of latitudinal cells in the model grid
        /// </summary>
        private UInt32 _NumLatCells;
        /// <summary>
        /// Get the number of latitudinal cells in the model grid
        /// </summary>
        public UInt32 NumLatCells
        {
            get { return _NumLatCells; }
        }
        
        /// <summary>
        /// The number of longitudinal cells in the model grid
        /// </summary>
        private UInt32 _NumLonCells;
        /// <summary>
        /// Get the number of longitudinal cells in the model grid
        /// </summary>
        public UInt32 NumLonCells
        {
            get { return _NumLonCells; }
        }
        
        /// <summary>
        /// The bottom (southern-most) latitude of each row of grid cells
        /// </summary>
        private float[] _Lats;
        /// <summary>
        /// Get the bottom (southern-most) latitude of each row of grid cells
        /// </summary>
        public float[] Lats
        {
            get { return _Lats; }
        }
        
        /// <summary>
        /// The left (western-most) longitude of each column of grid cells
        /// </summary>
        private float[] _Lons;
        /// <summary>
        /// Get the left (western-most) longitude of each column of grid cells
        /// </summary>
        public float[] Lons
        {
            get { return _Lons; }
        }

        /// <summary>
        /// Array of grid cells
        /// </summary>
        GridCell[,] InternalGrid;

        /// <summary>
        /// An array of lists of the functional group indices of each cohort to disperse. Array corresponds to grid cells. The lists correspond to individual cohorts to disperse.
        /// </summary>
        public List<uint>[,] DeltaFunctionalGroupDispersalArray;

        /// <summary>
        /// An array of lists of the positions within functional groups of each cohort to disperse. Array corresponds 
        /// to grid cells. The lists correspond to individual cohorts to disperse.
        /// </summary>
        public List<uint>[,] DeltaCohortNumberDispersalArray;

        /// <summary>
        /// An array of lists of paired longitude and latitude indices for the grid cells that each cohort will 
        /// to. Array corresponds to grid cells. The lists correspond to paired latitude and longitude indices that 
        /// each cohort will disperse to.
        /// </summary>
        public List<uint[]>[,] DeltaCellToDisperseToArray;
        
        /// <summary>
        /// An array of lists of cells that cohorts in a given grid cell can potentially disperse to (i.e. adjacent cells
        /// in the same realm). Array corresponds to focal grid cells. Lists correspond to cells that cohorts could
        /// disperse to from these focal cells.
        /// </summary>
        private List<uint[]>[,] CellsForDispersal;

        /// <summary>
        /// Analagous to the array of lists CellsForDispersal, but instead of containing the identities of the cells that are dispersable to,
        /// instead each array element contains a uint list which is coded to correspond to directions:
        /// 1. N, 2. NE, 3. E, 4. SE, 5. S, 6. SW, 7 W, 8, NW.
        /// Each item in the list corresponds to the analagous item in CellsForDispersal, and indicates to which direction the cell for dispersal lies. 
        /// This is used by the advective dispersal class, in order to check whether advective dispersal in a particular direction can actually occur.
        /// </summary>
        private List<uint>[,] CellsForDispersalDirection;

        /// <summary>
        /// The heights of grid cells in each latitudinal band
        /// </summary>
        public double[] CellHeightsKm;
        /// <summary>
        /// The widths of grid cells in each latitudinal band
        /// </summary>
        public double[] CellWidthsKm;

        /// <summary>
        /// An instance of the simple random number generator class
        /// </summary>
        private NonStaticSimpleRNG RandomNumberGenerator = new NonStaticSimpleRNG();

        /// <summary>
        /// Instance of the class to perform general functions
        /// </summary>
        private UtilityFunctions Utilities;

        /// <summary>
        /// Thread-local variables for tracking extinction and production of cohorts
        /// </summary>
        /// <todo>Needs a little tidying and checking of access levels</todo>
        private class ThreadLockedParallelVariablesModelGrid
        {
            /// <summary>
            /// Thread-locked variable to track the cohort ID to assign to newly produced cohorts
            /// </summary>
            public Int64 NextCohortIDThreadLocked;

        }

        /// <summary>
        /// Constructor for model grid: assigns grid properties and initialises the grid cells
        /// </summary>
        /// <param name="minLat">Minimum grid latitude (degrees)</param>
        /// <param name="minLon">Minimum grid longitude (degrees, currently -180 to 180)</param>
        /// <param name="maxLat">Maximum grid latitude (degrees)</param>
        /// <param name="maxLon">Maximum grid longitude (degrees, currently -180 to 180)</param>
        /// <param name="latCellSize">Latitudinal resolution of grid cell</param>
        /// <param name="lonCellSize">Longitudinal resolution of grid cell</param>
        /// <param name="cellRarefaction">The rarefaction to be applied to active grid cells in the model</param>
        /// <param name="enviroStack">Environmental data layers</param>
        /// <param name="cohortFunctionalGroups">The functional group definitions for cohorts in the model</param>
        /// <param name="stockFunctionalGroups">The functional group definitions for stocks in the model</param>
        /// <param name="globalDiagnostics">Global daignostic variables</param>
        /// <param name="tracking">Whether process-tracking is enabled</param>
        /// <param name="DrawRandomly">Whether the model is set to use a random draw</param>
        /// <param name="specificLocations">Whether the model is to be run for specific locations</param>
        public ModelGrid(float minLat, float minLon,float maxLat,float maxLon,float latCellSize,float lonCellSize, int cellRarefaction, 
            SortedList<string,EnviroData> enviroStack, FunctionalGroupDefinitions cohortFunctionalGroups, FunctionalGroupDefinitions
            stockFunctionalGroups, SortedList<string, double> globalDiagnostics, Boolean tracking, Boolean DrawRandomly, 
            Boolean specificLocations)
        {
            // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
            NumGrids = NumGrids + 1;
            //Debug.Assert(NumGrids < 2, "You have initialised more than one grid on which to apply models. At present, this is not supported");

            // Initialise the utility functions
            Utilities = new UtilityFunctions();

            // Seed the random number generator
            // Set the seed for the random number generator
            RandomNumberGenerator = new NonStaticSimpleRNG();
            if (DrawRandomly)
            {
                RandomNumberGenerator.SetSeedFromSystemTime();
            }
            else
            {
                RandomNumberGenerator.SetSeed(4315);
            }

            // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
            _MinLatitude = minLat;
            _MinLongitude = minLon;
            _MaxLatitude = maxLat;
            _MaxLongitude = maxLon;
            _LatCellSize = latCellSize;
            _LonCellSize = lonCellSize;
            _GridCellRarefaction = cellRarefaction;

            // Check to see if the number of grid cells is an integer
            Debug.Assert((((_MaxLatitude - _MinLatitude) % _LatCellSize) == 0), "Error: number of grid cells is non-integer: check cell size");

            
            _NumLatCells = (UInt32)((_MaxLatitude - _MinLatitude) / _LatCellSize);
            _NumLonCells = (UInt32)((_MaxLongitude - _MinLongitude) / _LonCellSize);
            _Lats = new float[_NumLatCells];
            _Lons = new float[_NumLonCells];

            // Set up latitude and longitude vectors - lower left
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                _Lats[ii] = _MinLatitude + ii * _LatCellSize;
            }
            for (int jj = 0; jj < _NumLonCells; jj++)
            {
                _Lons[jj] = _MinLongitude + jj * _LonCellSize;
            }
            

            // Instantiate a grid of grid cells
            InternalGrid = new GridCell[_NumLatCells, _NumLonCells];

            // Instantiate the arrays of lists of cohorts to disperse
            DeltaFunctionalGroupDispersalArray = new List<uint>[_NumLatCells, _NumLonCells];
            DeltaCohortNumberDispersalArray = new List<uint>[_NumLatCells, _NumLonCells];

            // Instantiate the array of lists of grid cells to disperse those cohorts to
            DeltaCellToDisperseToArray = new List<uint[]>[_NumLatCells, _NumLonCells];

            // An array of lists of cells to which organisms in each cell can disperse to; includes all cells which contribute to the 
            // perimeter list, plus diagonal cells if they are in the same realm
            CellsForDispersal = new List<uint[]>[_NumLatCells, _NumLonCells];

            // An array of lists of directions corresponding to cells which organisms can disperse to
            CellsForDispersalDirection = new List<uint>[_NumLatCells, _NumLonCells];

            Console.WriteLine("Initialising grid cell environment:");


            // Loop through to set up model grid
            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
            {
                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
                {
                    InternalGrid[ii, jj] = new GridCell(_Lats[ii],(uint)ii, _Lons[jj],(uint)jj, LatCellSize, LonCellSize, enviroStack,
                        GlobalMissingValue, cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
                    CellsForDispersal[ii,jj] = new List<uint[]>();
                    CellsForDispersalDirection[ii, jj] = new List<uint>();
                    Console.Write("\rRow {0} of {1}", ii+1, NumLatCells/GridCellRarefaction);
                }
            }
            Console.WriteLine("");
            Console.WriteLine("");


            InterpolateMissingValues();


            // Fill in the array of dispersable perimeter lengths for each grid cell
            CalculatePerimeterLengthsAndCellsDispersableTo();

            CellHeightsKm = new double[_Lats.Length];
            CellWidthsKm = new double[_Lats.Length];

            // Calculate the lengths of widths of grid cells in each latitudinal strip
            // Assume that we are at the midpoint of each cell when calculating lengths
            for (int ii = 0; ii < _Lats.Length; ii++)
            {
                 CellHeightsKm[ii] = Utilities.CalculateLengthOfDegreeLatitude(_Lats[ii] + _LatCellSize / 2) * _LatCellSize;
                 CellWidthsKm[ii] = Utilities.CalculateLengthOfDegreeLongitude(_Lats[ii] + _LatCellSize / 2) * _LonCellSize;
            }
        }

        /// <summary>
        /// Overloaded constructor for model grid to construct the grid for specific locations
        /// </summary>
        /// <param name="minLat">Minimum grid latitude (degrees)</param>
        /// <param name="minLon">Minimum grid longitude (degrees, currently -180 to 180)</param>
        /// <param name="maxLat">Maximum grid latitude (degrees)</param>
        /// <param name="maxLon">Maximum grid longitude (degrees, currently -180 to 180)</param>
        /// <param name="latCellSize">Latitudinal size of grid cells</param>
        /// <param name="lonCellSize">Longitudinal size of grid cells</param>
        /// <param name="cellList">List of indices of active cells in the model grid</param>
        /// <param name="enviroStack">List of environmental data layers</param>
        /// <param name="cohortFunctionalGroups">The functional group definitions for cohorts in the model</param>
        /// <param name="stockFunctionalGroups">The functional group definitions for stocks in the model</param>
        /// <param name="globalDiagnostics">Global diagnostic variables</param>
        /// <param name="tracking">Whether process tracking is enabled</param>
        /// <param name="specificLocations">Whether the model is to be run for specific locations</param>
        /// <param name="runInParallel">Whether model grid cells will be run in parallel</param>
        public ModelGrid(float minLat, float minLon, float maxLat, float maxLon, float latCellSize, float lonCellSize, List<uint[]> cellList, 
            SortedList<string, EnviroData> enviroStack, FunctionalGroupDefinitions cohortFunctionalGroups,
            FunctionalGroupDefinitions stockFunctionalGroups, SortedList<string, double> globalDiagnostics, Boolean tracking, 
            Boolean specificLocations, Boolean runInParallel)
        { 
            // Add one to the counter of the number of grids. If there is more than one model grid, exit the program with a debug crash.
            NumGrids = NumGrids + 1;
            //Debug.Assert(NumGrids < 2, "You have initialised more than one grid on which to apply models. At present, this is not supported");

            // Initialise the utility functions
            Utilities = new UtilityFunctions();

            // CURRENTLY DEFINING MODEL CELLS BY BOTTOM LEFT CORNER
            _MinLatitude = minLat;
            _MinLongitude = minLon;
            _MaxLatitude = maxLat;
            _MaxLongitude = maxLon;
            _LatCellSize = latCellSize;
            _LonCellSize = lonCellSize;
            _GridCellRarefaction = 1;

            // Check to see if the number of grid cells is an integer
            Debug.Assert((((_MaxLatitude - _MinLatitude) % _LatCellSize) == 0), "Error: number of grid cells is non-integer: check cell size");


            _NumLatCells = (UInt32)((_MaxLatitude - _MinLatitude) / _LatCellSize);
            _NumLonCells = (UInt32)((_MaxLongitude - _MinLongitude) / _LonCellSize);
            _Lats = new float[_NumLatCells];
            _Lons = new float[_NumLonCells];

            // Set up latitude and longitude vectors - lower left
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                _Lats[ii] = _MinLatitude + ii * _LatCellSize;
            }
            for (int jj = 0; jj < _NumLonCells; jj++)
            {
                _Lons[jj] = _MinLongitude + jj * _LonCellSize;
            }


            // Set up a grid of grid cells
            InternalGrid = new GridCell[_NumLatCells, _NumLonCells];

            // Instantiate the arrays of lists of cohorts to disperse
            DeltaFunctionalGroupDispersalArray = new List<uint>[_NumLatCells, _NumLonCells];
            DeltaCohortNumberDispersalArray = new List<uint>[_NumLatCells, _NumLonCells];

            // Instantiate the array of lists of grid cells to disperse those cohorts to
            DeltaCellToDisperseToArray = new List<uint[]>[_NumLatCells, _NumLonCells];


            // An array of lists of cells to which organisms in each cell can disperse to; includes all cells which contribute to the 
            // perimeter list, plus diagonal cells if they are in the same realm
            CellsForDispersal = new List<uint[]>[_NumLatCells, _NumLonCells];

            // An array of lists of directions corresponding to cells which organisms can disperse to
            CellsForDispersalDirection = new List<uint>[_NumLatCells, _NumLonCells];

            Console.WriteLine("Initialising grid cell environment:");


            int Count = 0;

            int NCells = cellList.Count;

            if (!runInParallel)
            {
                // Loop over cells to set up the model grid
                for (int ii = 0; ii < cellList.Count; ii++)
                {
                    // Create the grid cell at the specified position
                    InternalGrid[cellList[ii][0], cellList[ii][1]] = new GridCell(_Lats[cellList[ii][0]], cellList[ii][0],
                        _Lons[cellList[ii][1]], cellList[ii][1], latCellSize, lonCellSize, enviroStack, _GlobalMissingValue,
                        cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
                    if (!specificLocations)
                    {
                        CellsForDispersal[cellList[ii][0], cellList[ii][1]] = new List<uint[]>();
                        CellsForDispersalDirection[cellList[ii][0], cellList[ii][1]] = new List<uint>();
                        Count++;
                        Console.Write("\rInitialised {0} of {1}", Count, NCells);
                    }
                    else
                    {
                        Console.Write("\rRow {0} of {1}", ii + 1, NumLatCells / GridCellRarefaction);
                        Console.WriteLine("");
                        Console.WriteLine("");
                    }
                }
            }
            else
            {

                // Run a parallel loop over rows
                Parallel.For(0, NCells, ii =>
                {
                    // Create the grid cell at the specified position
                    InternalGrid[cellList[ii][0], cellList[ii][1]] = new GridCell(_Lats[cellList[ii][0]], cellList[ii][0],
                        _Lons[cellList[ii][1]], cellList[ii][1], latCellSize, lonCellSize, enviroStack, _GlobalMissingValue,
                        cohortFunctionalGroups, stockFunctionalGroups, globalDiagnostics, tracking, specificLocations);
                    if (!specificLocations)
                    {
                        CellsForDispersal[cellList[ii][0], cellList[ii][1]] = new List<uint[]>();
                        CellsForDispersalDirection[cellList[ii][0], cellList[ii][1]] = new List<uint>();
                    }

                    Count++;
                    Console.Write("\rInitialised {0} of {1}", Count, NCells);
                }
                 );

            }


            if (!specificLocations)
            {
                InterpolateMissingValues();


                // Ex 1
                //ConvertToConstantValues(2, "Temperature");
                //ConvertToConstantValues(2, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(2, "DiurnalTemperatureRange");
                //ConvertToConstantValues(2, "TemperatureNPP");
                //ConvertToConstantValues(2, "NPP");
                //ConvertToConstantValues(2, "Precipitation");
                //ConvertToConstantValues(2, "Fraction Year Frost");
                //ConvertToConstantValues(2, "AWC");
                //ConvertToConstantValues(2, "AET");
                //ConvertToConstantValues(2, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();


                //Ex 2
                //ConvertToConstantValues(2, "Temperature");
                //ConvertToConstantValues(2, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(2, "DiurnalTemperatureRange");
                //ConvertToMeridionalAverage(true, "TemperatureNPP");
                //ConvertToMeridionalAverage(true, "NPP");
                //ConvertToMeridionalAverage(true, "Precipitation");
                //ConvertToMeridionalAverage(true, "Fraction Year Frost");
                //ConvertToMeridionalAverage(true, "AWC");
                //ConvertToMeridionalAverage(true, "AET");
                //ConvertToMeridionalAverage(true, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();

                //Ex 3
                //ConvertToConstantValues(2, "Temperature");
                //ConvertToConstantValues(2, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(2, "DiurnalTemperatureRange");
                //ConvertToConstantValues(2, "TemperatureNPP");
                //ConvertToConstantValues(2, "NPP");
                //ConvertToConstantValues(2, "Precipitation");
                //ConvertToConstantValues(2, "Fraction Year Frost");
                //ConvertToConstantValues(2, "AWC");
                //ConvertToConstantValues(2, "AET");
                //ConvertToConstantValues(2, "Fraction Year Fire");
                ///* DON'T calculate seasonality with these average conditions to retain 'latitudinal' seasonality patterns
                //ConvertToSeasonalConstantValues(true, "TemperatureNPP");
                //ConvertToSeasonalConstantValues(true, "NPP");
                //ConvertToSeasonalConstantValues(true, "Precipitation");
                //ConvertToSeasonalConstantValues(true, "Fraction Year Frost");
                //ConvertToSeasonalConstantValues(true, "AWC");
                //ConvertToSeasonalConstantValues(true, "AET");
                //ConvertToSeasonalConstantValues(true, "Fraction Year Fire");
                //*/
                //CalculateBreedingSeason();


                //Ex 4
                //ConvertToMeridionalAverage(true, "Temperature");
                //ConvertToMeridionalAverage(true, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToMeridionalAverage(true, "DiurnalTemperatureRange");
                //ConvertToConstantValues(2, "TemperatureNPP");
                //ConvertToConstantValues(2, "NPP");
                //ConvertToConstantValues(2, "Precipitation");
                //ConvertToConstantValues(2, "Fraction Year Frost");
                //ConvertToConstantValues(2, "AWC");
                //ConvertToConstantValues(2, "AET");
                //ConvertToConstantValues(2, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();

                ////Ex 5
                //ConvertToAveragedSeasonalValues(2, "Temperature");
                //ConvertToAveragedSeasonalValues(2, "DiurnalTemperatureRange");
                ///* DON'T do anything with SD Temp or Annual Temp. These inform the seasonality of
                //// * metabolism in each cell
                ////ConvertToAveragedSeasonalValues(2, "SDTemperature");
                //// */
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(2, "TemperatureNPP");
                //ConvertToConstantValues(2, "NPP");
                //ConvertToConstantValues(2, "Precipitation");
                //ConvertToConstantValues(2, "Fraction Year Frost");
                //ConvertToConstantValues(2, "AWC");
                //ConvertToConstantValues(2, "AET");
                //ConvertToConstantValues(2, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();


                //Ex 6
                //ConvertToConstantValues(3, "Temperature");
                //ConvertToConstantValues(3, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(3, "DiurnalTemperatureRange");
                //ConvertToConstantValues(3, "TemperatureNPP");
                //ConvertToConstantValues(3, "NPP");
                //ConvertToConstantValues(3, "Precipitation");
                //ConvertToConstantValues(3, "Fraction Year Frost");
                //ConvertToConstantValues(3, "AWC");
                //ConvertToConstantValues(3, "AET");
                //ConvertToConstantValues(3, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();

                //Ex 7
                //ConvertToConstantValues(3, "Temperature");
                //ConvertToConstantValues(3, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(3, "DiurnalTemperatureRange");
                //ConvertToMeridionalAverage(true, "TemperatureNPP");
                //ConvertToMeridionalAverage(true, "NPP");
                //ConvertToMeridionalAverage(true, "Precipitation");
                //ConvertToMeridionalAverage(true, "Fraction Year Frost");
                //ConvertToMeridionalAverage(true, "AWC");
                //ConvertToMeridionalAverage(true, "AET");
                //ConvertToMeridionalAverage(true, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();

                //Ex 8
                //ConvertToConstantValues(3, "Temperature");
                //ConvertToConstantValues(3, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToConstantValues(3, "DiurnalTemperatureRange");
                //ConvertToSeasonalConstantValues(true, "TemperatureNPP");
                //ConvertToSeasonalConstantValues(true, "NPP");
                //ConvertToSeasonalConstantValues(true, "Precipitation");
                //ConvertToSeasonalConstantValues(true, "Fraction Year Frost");
                //ConvertToSeasonalConstantValues(true, "AWC");
                //ConvertToSeasonalConstantValues(true, "AET");
                //ConvertToSeasonalConstantValues(true, "Fraction Year Fire");
                //CalculateBreedingSeason();

                //Exp 9
                //ConvertToMeridionalAverage(true, "Temperature");
                //ConvertToMeridionalAverage(true, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToMeridionalAverage(true, "DiurnalTemperatureRange");
                //ConvertToConstantValues(3, "TemperatureNPP");
                //ConvertToConstantValues(3, "NPP");
                //ConvertToConstantValues(3, "Precipitation");
                //ConvertToConstantValues(3, "Fraction Year Frost");
                //ConvertToConstantValues(3, "AWC");
                //ConvertToConstantValues(3, "AET");
                //ConvertToConstantValues(3, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();


                //Ex 10
                //ConvertToSeasonalConstantValues(true, "Temperature");
                //ConvertToSeasonalConstantValues(true, "SDTemperature");
                //CalculateAnnualMean("AnnualTemperature", "Temperature");
                //ConvertToSeasonalConstantValues(true, "DiurnalTemperatureRange");
                //ConvertToConstantValues(3, "TemperatureNPP");
                //ConvertToConstantValues(3, "NPP");
                //ConvertToConstantValues(3, "Precipitation");
                //ConvertToConstantValues(3, "Fraction Year Frost");
                //ConvertToConstantValues(3, "AWC");
                //ConvertToConstantValues(3, "AET");
                //ConvertToConstantValues(3, "Fraction Year Fire");
                //CalculateNPPSeasonality();
                //CalculateBreedingSeason();


                // Fill in the array of dispersable perimeter lengths for each grid cell
                CalculatePerimeterLengthsAndCellsDispersableTo();

                CellHeightsKm = new double[_Lats.Length];
                CellWidthsKm = new double[_Lats.Length];

                // Calculate the lengths of widths of grid cells in each latitudinal strip
                // Assume that we are at the midpoint of each cell when calculating lengths
                for (int ii = 0; ii < _Lats.Length; ii++)
                {
                    CellHeightsKm[ii] = Utilities.CalculateLengthOfDegreeLatitude(_Lats[ii] + _LatCellSize / 2) * _LatCellSize;
                    CellWidthsKm[ii] = Utilities.CalculateLengthOfDegreeLongitude(_Lats[ii] + _LatCellSize / 2) * _LonCellSize;
                }
            }

            Console.WriteLine("\n");

        }

        /// <summary>
        /// Estimates missing environmental data for grid cells by interpolation
        /// </summary>
        public void InterpolateMissingValues()
        {
            SortedList<string, double[]> WorkingCellEnvironment = new SortedList<string, double[]>();
            Boolean Changed = false;

            for (uint ii = 0; ii < _NumLatCells; ii++)
            {
                for (uint jj = 0; jj < _NumLonCells; jj++)
                {
                    WorkingCellEnvironment = GetCellEnvironment(ii, jj);

                    // If the cell environment does not contain valid NPP data then interpolate values
                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]))
                    {
                        //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
                        WorkingCellEnvironment["NPP"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "NPP", WorkingCellEnvironment["Realm"][0]);
                        
                        //Calculate NPP seasonality - for use in converting annual NPP estimates to monthly
                        WorkingCellEnvironment["Seasonality"] = InternalGrid[ii, jj].CalculateNPPSeasonality(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
                        Changed = true;
                    }
                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                    else
                    {
                        WorkingCellEnvironment["NPP"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["NPP"], WorkingCellEnvironment["Missing Value"][0]);
                    }

                    // If the cell environment does not contain valid monthly mean diurnal temperature range data then interpolate values
                    if (InternalGrid[ii, jj].ContainsMissingValue(WorkingCellEnvironment["DiurnalTemperatureRange"], WorkingCellEnvironment["Missing Value"][0]))
                    {
                        //If NPP doesn't exist the interpolate from surrounding values (of the same realm)
                        WorkingCellEnvironment["DiurnalTemperatureRange"] = FillWithInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "DiurnalTemperatureRange", WorkingCellEnvironment["Realm"][0]);

                        Changed = true;
                    }

                    // Same for u and v velocities
                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0]))
                    {
                        //If u doesn't exist the interpolate from surrounding values (of the same realm)
                        WorkingCellEnvironment["uVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "uVel", WorkingCellEnvironment["Realm"][0]);

                        Changed = true;
                    }
                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                    else
                    {
                        WorkingCellEnvironment["uVel"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["uVel"], WorkingCellEnvironment["Missing Value"][0]);
                    }

                    if (!InternalGrid[ii, jj].ContainsData(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0]))
                    {
                        //If v vel doesn't exist the interpolate from surrounding values (of the same realm)
                        WorkingCellEnvironment["vVel"] = GetInterpolatedValues(ii, jj, GetCellLatitude(ii), GetCellLongitude(jj), "vVel", WorkingCellEnvironment["Realm"][0]);

                        Changed = true;
                    }
                    // Otherwise convert the missing data values to zeroes where they exist amongst valid data eg in polar regions.
                    else
                    {
                        WorkingCellEnvironment["vVel"] = InternalGrid[ii, jj].ConvertMissingValuesToZero(WorkingCellEnvironment["vVel"], WorkingCellEnvironment["Missing Value"][0]);
                    }
                    
                    if(Changed) InternalGrid[ii, jj].CellEnvironment = WorkingCellEnvironment;
                }
            }
        }

        /// <summary>
        /// Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
        /// non missing data values
        /// </summary>
        /// <param name="latIndex">Index of the latitude cell for which the weighted average over surrounding cells is requested</param>
        /// <param name="lonIndex">Index of the longitude cell for which the weighted average over surrounding cells is requested</param>
        /// <param name="lat">Latitude of the cell for which the weighted value is requested</param>
        /// <param name="lon">Longitude of the cell for which the weighted value is requested</param>
        /// <param name="dataName">Names of the data for which weighted value is requested</param>
        /// <param name="realm">Realm of the grid cell for which data is to be averaged over</param>
        /// <returns>The weighted average value of the specified data type across surrounding grid cells of the specified realm</returns>
        private double[] GetInterpolatedValues(uint latIndex, uint lonIndex, double lat, double lon, string dataName, double realm)
        {
            SortedList<string, double[]> TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
            double[] InterpData = new double[TempCellEnvironment[dataName].Length];
            uint[] InterpCount = new uint[TempCellEnvironment[dataName].Length];

            uint LowerLatIndex = latIndex - 1;
            uint UpperLatIndex = latIndex + 1;
            uint LowerLonIndex = lonIndex - 1;
            uint UpperLonIndex = lonIndex + 1;


            if (latIndex == 0) LowerLatIndex = latIndex;
            if (lat.CompareTo(this.MaxLatitude) == 0) UpperLatIndex = latIndex;

            if (lonIndex == 0) LowerLonIndex = lonIndex;
            if (lon.CompareTo(this.MaxLongitude) == 0) UpperLonIndex = lonIndex;

            //Loop over surrounding cells in the datalayer
            for (uint ii = LowerLatIndex; ii <= UpperLatIndex; ii++)
            {
                for (uint jj = LowerLonIndex; jj < UpperLonIndex; jj++)
                {
                    if (ii < _NumLatCells && jj < _NumLonCells)
                    {
                        TempCellEnvironment = GetCellEnvironment(ii, jj);

                        for (uint hh = 0; hh < InterpData.Length; hh++)
                        {
                            //If the cell contains data then sum this and increment count
                            if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm)
                            {
                                InterpData[hh] += TempCellEnvironment[dataName][hh];
                                InterpCount[hh]++;
                            }
                        }
                    }
                }
            }

            //take the mean over surrounding valid cells for each timestep
            for (int hh = 0; hh < InterpData.Length; hh++)
            {
                if (InterpCount[hh] > 0)
                {
                    InterpData[hh] /= InterpCount[hh];
                }
                else
                {
                    InterpData[hh] = 0.0;
                }
            }
            return InterpData;
        }

        /// <summary>
        /// Calculate the weighted average of surrounding grid cell data, where those grid cells are of the specified realm and contain
        /// non missing data values
        /// </summary>
        /// <param name="latIndex">Index of the latitude cell for which the weighted average over surrounding cells is requested</param>
        /// <param name="lonIndex">Index of the longitude cell for which the weighted average over surrounding cells is requested</param>
        /// <param name="lat">Latitude of the cell for which the weighted value is requested</param>
        /// <param name="lon">Longitude of the cell for which the weighted value is requested</param>
        /// <param name="dataName">Names of the data for which weighted value is requested</param>
        /// <param name="realm">Realm of the grid cell for which data is to be averaged over</param>
        /// <returns>The weighted average value of the specified data type across surrounding grid cells of the specified realm</returns>
        private double[] FillWithInterpolatedValues(uint latIndex, uint lonIndex, double lat, double lon, string dataName, double realm)
        {
            SortedList<string, double[]> TempCellEnvironment = GetCellEnvironment(latIndex, lonIndex);
            double[] InterpData = new double[TempCellEnvironment[dataName].Length];
            uint[] InterpCount = new uint[TempCellEnvironment[dataName].Length];
            uint LowerLatIndex = latIndex - 1;
            uint UpperLatIndex = latIndex + 1;
            uint LowerLonIndex = lonIndex - 1;
            uint UpperLonIndex = lonIndex + 1;


            if (latIndex == 0) LowerLatIndex = latIndex;
            if (lat.CompareTo(this.MaxLatitude) == 0) UpperLatIndex = latIndex;

            if (lonIndex == 0) LowerLonIndex = lonIndex;
            if (lon.CompareTo(this.MaxLongitude) == 0) UpperLonIndex = lonIndex;

            for (uint hh = 0; hh < InterpData.Length; hh++)
            {
                if (TempCellEnvironment[dataName][hh] == TempCellEnvironment["Missing Value"][0])
                {
                    //Loop over surrounding cells in the datalayer
                    for (uint ii = LowerLatIndex; ii <= UpperLatIndex; ii++)
                    {
                        for (uint jj = LowerLonIndex; jj <= UpperLonIndex; jj++)
                        {
                            if (ii < _NumLatCells && jj < _NumLonCells)
                            {
                                TempCellEnvironment = GetCellEnvironment(ii, jj);

                                //If the cell contains data then sum this and increment count
                                if (TempCellEnvironment[dataName][hh] != TempCellEnvironment["Missing Value"][0] && TempCellEnvironment["Realm"][0] == realm)
                                {
                                    InterpData[hh] += TempCellEnvironment[dataName][hh];
                                    InterpCount[hh]++;
                                }

                            }
                        }
                    }
                    //take the mean over surrounding valid cells for each timestep
                    if (InterpCount[hh] > 0)
                    {
                        InterpData[hh] /= InterpCount[hh];
                    }
                    else
                    {
                        InterpData[hh] = 0.0;
                    }
                }
                else
                {
                    InterpData[hh] = TempCellEnvironment[dataName][hh];
                }
            }


            return InterpData;
        }


        public void CalculateAnnualMean(string outLayerName, string dataLayerName)
        {
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    InternalGrid[ii, jj].CellEnvironment[outLayerName] =
                        InternalGrid[ii, jj].CalculateAnnualMean(
                            InternalGrid[ii, jj].CellEnvironment[dataLayerName], _GlobalMissingValue);
                }
            }
        }

        public void CalculateBreedingSeason()
        {
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    InternalGrid[ii, jj].CellEnvironment["Breeding Season"] =
                        InternalGrid[ii, jj].CalculateBreedingSeason();
                }
            }
        }

        public void CalculateNPPSeasonality()
        {

            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                    InternalGrid[ii, jj].CellEnvironment["Seasonality"] = 
                        InternalGrid[ii,jj].CalculateNPPSeasonality(
                            InternalGrid[ii, jj].CellEnvironment["NPP"],_GlobalMissingValue);
            }
        }

        /// <summary>
        /// Convert an environmental layer to either equatorial or annual mean values in each timestep
        /// </summary>
        /// <param name="equatorial">Use equatorial as the constant value</param>
        /// <param name="layerName">The name of the layer to be converted</param>
        public void ConvertToAveragedSeasonalValues(int spatialDomain, string layerName)
        {

            int NumTimsteps = GetCellEnvironment(1, 1)[layerName].Length;
            double[,] GlobalRealm = new double[_NumLatCells, _NumLonCells];
            double[,] Average = new double[_NumLatCells, _NumLonCells];
            double[] tempVector;
            

            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    GlobalRealm[ii, jj] = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
                    //for each cell calculate the SD of temperature
                    tempVector = InternalGrid[ii, jj].CellEnvironment[layerName];
                    var subset = from i in tempVector where i != _GlobalMissingValue select i;
                    Average[ii,jj] = subset.Sum() / subset.Count();

                    if (double.IsNaN(Average[ii,jj]) || double.IsInfinity(Average[ii,jj]))
                    {
                        Average[ii,jj] = _GlobalMissingValue;
                    }

                }
            }

            double[] AverageVals = CalculateAnnualAverageEnvironment(spatialDomain, layerName, GlobalRealm);

            
            //Now transform the environmental timeseries in each cell to have the average value for the spatial domain
            //requested but maintain the temporal pattern and SD.
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    if (Average[ii, jj] != _GlobalMissingValue)
                    {
                        tempVector = InternalGrid[ii, jj].CellEnvironment[layerName];
                        for (int n = 0; n < NumTimsteps; n++)
                        {
                            if (tempVector[n] != _GlobalMissingValue)
                            {
                                InternalGrid[ii, jj].CellEnvironment[layerName][n] = AverageVals[(int)GlobalRealm[ii, jj]-1] + tempVector[n] - Average[ii, jj];
                            }
                            else
                            {
                                InternalGrid[ii, jj].CellEnvironment[layerName][n] = _GlobalMissingValue;
                            }
                        }
                    }
                    
                }
            }


            
            /*double[, ,] AnnualMnArray = new double[2, _NumLatCells, _NumLonCells];
            double[,] GlobalRealm = new double[_NumLatCells, _NumLonCells];
            double[] TempEnv;

            double[] AverageVal1 = new double[NumTimsteps];
            double[] AverageVal2 = new double[NumTimsteps];

            for (int n = 0; n < NumTimsteps; n++)
            {
                AverageVal1[n] = _GlobalMissingValue;
                AverageVal2[n] = _GlobalMissingValue;
            }

            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                    GlobalRealm[ii, jj] = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
            }

            //Calculate the average conditions to set environments to
            //Equatorial average so set the environment layer to have a mean value equal to that at the equator but with seasonality
            //As per the current environment
            if (spatialDomain == 2)
            {
                uint ii = GetLatIndex(0.0);
                double[] Temp1 = new double[_NumLonCells];
                double[] Temp2 = new double[_NumLonCells];

                for (int n = 0; n < NumTimsteps; n++)
                {


                    for (int jj = 0; jj < _NumLonCells; jj++)
                    {
                        Temp1[jj] = _GlobalMissingValue;
                        Temp2[jj] = _GlobalMissingValue;
                    }

                    //for this latitude band calculate the annual average conditions across realms
                    for (int jj = 0; jj < _NumLonCells; jj++)
                    {
                        TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];

                        if (GlobalRealm[ii, jj] == 1.0)
                        {
                            Temp1[jj] = TempEnv[n];
                        }
                        else
                        {
                            Temp2[jj] = TempEnv[n];
                        }

                    }

                    var subset1 = from i in Temp1.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                    AverageVal1[n] = subset1.Sum() / subset1.Count();
                    if (double.IsNaN(AverageVal1[n]) || double.IsInfinity(AverageVal1[n]))
                        AverageVal1[n] = GlobalMissingValue;

                    subset1 = from i in Temp2.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                    AverageVal2[n] = subset1.Sum() / subset1.Count();
                    if (double.IsNaN(AverageVal2[n]) || double.IsInfinity(AverageVal2[n]))
                        AverageVal2[n] = GlobalMissingValue;
                }

            }
            //Equatorial average
            if (spatialDomain == 1)
            {

            }
            else
            {
                double[,] Temp1 = new double[_NumLatCells, _NumLonCells];
                double[,] Temp2 = new double[_NumLatCells, _NumLonCells];

                for (int n = 0; n < NumTimsteps; n++)
                {


                    for (int ii = 0; ii < _NumLatCells; ii++)
                    {
                        for (int jj = 0; jj < _NumLonCells; jj++)
                        {
                            Temp1[ii, jj] = _GlobalMissingValue;
                            Temp2[ii, jj] = _GlobalMissingValue;
                        }
                    }

                    for (int ii = 0; ii < _NumLatCells; ii++)
                    {
                        for (int jj = 0; jj < _NumLonCells; jj++)
                        {
                            TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];

                            if (GlobalRealm[ii, jj] == 1.0)
                            {
                                Temp1[ii, jj] = TempEnv[n];

                            }
                            else
                            {
                                Temp2[ii, jj] = TempEnv[n];

                            }

                        }

                    }
                    var subset1 = from i in Temp1.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                    AverageVal1[n] = subset1.Sum() / subset1.Count();
                    if (AverageVal1[n] == 0 || double.IsNaN(AverageVal1[n]) || double.IsInfinity(AverageVal1[n]))
                        AverageVal1[n] = GlobalMissingValue;

                    subset1 = from i in Temp2.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                    AverageVal2[n] = subset1.Sum() / subset1.Count();
                    if (AverageVal2[n] == 0 || double.IsNaN(AverageVal2[n]) || double.IsInfinity(AverageVal2[n]))
                        AverageVal2[n] = GlobalMissingValue;
                }
            }

            //Now set the cell environments to these conditions
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    for (int n = 0; n < NumTimsteps; n++)
                    {

                        if (GlobalRealm[ii, jj] == 1.0)
                        {
                            InternalGrid[ii, jj].CellEnvironment[layerName][n] =
                                AverageVal1[n];
                        }
                        else
                        {
                            InternalGrid[ii, jj].CellEnvironment[layerName][n] =
                                AverageVal2[n];
                        }

                    }
                }
            }*/

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialDomain">1 = Global, 2 = Equatorial, 3 = 45N</param>
        /// <param name="layerName">The layer to be converted</param>
        public void ConvertToConstantValues(int spatialDomain, string layerName)
        {
            int NumTimsteps = GetCellEnvironment(1, 1)[layerName].Length;
            double[,] GlobalRealm = new double[_NumLatCells, _NumLonCells];

            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                    GlobalRealm[ii, jj] = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
            }
            
            double[] AverageVals = CalculateAnnualAverageEnvironment(spatialDomain, layerName, GlobalRealm);

            //Now set the cell environments to these conditions
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    for (int n = 0; n < NumTimsteps; n++)
                    {
                     
                        if (GlobalRealm[ii, jj] == 1.0)
                        {
                            InternalGrid[ii, jj].CellEnvironment[layerName][n] = 
                                AverageVals[0];
                        }
                        else
                        {
                            InternalGrid[ii, jj].CellEnvironment[layerName][n] = 
                                AverageVals[1];
                        }

                    }
                }
            }

        }


        public double[] CalculateAnnualAverageEnvironment(int spatialDomain, string layerName, double[,] globalRealm)
        {
            double[, ,] AnnualMnArray = new double[2, _NumLatCells, _NumLonCells];
            double[] TempEnv;

            double AverageVal1 = _GlobalMissingValue;
            double AverageVal2 = _GlobalMissingValue;

            
            //Calculate the average conditions to set environments to
            //Equatorial average
            if (spatialDomain == 2)
            {
                uint ii = GetLatIndex(0.0);
                double[] Temp1 = new double[_NumLonCells];
                double[] Temp2 = new double[_NumLonCells];

                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    Temp1[jj] = _GlobalMissingValue;
                    Temp2[jj] = _GlobalMissingValue;
                }

                //for this latitude band calculate the annual average conditions across realms
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];

                    if (globalRealm[ii, jj] == 1.0)
                    {
                        var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                        Temp1[jj] = subset.Sum() / subset.Count();

                        if (double.IsNaN(Temp1[jj]) ||
                            double.IsInfinity(Temp1[jj]))
                            Temp1[jj] = _GlobalMissingValue;
                    }
                    else
                    {
                        var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                        Temp2[jj] = subset.Sum() / subset.Count();

                        if (double.IsNaN(Temp2[jj]) ||
                            double.IsInfinity(Temp2[jj]))
                            Temp2[jj] = _GlobalMissingValue;
                    }

                }

                var subset1 = from i in Temp1.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal1 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal1) || double.IsInfinity(AverageVal1))
                    AverageVal1 = GlobalMissingValue;

                subset1 = from i in Temp2.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal2 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal2) || double.IsInfinity(AverageVal2))
                    AverageVal2 = GlobalMissingValue;


            }
            // Take a high latitude average
            else if (spatialDomain == 3)
            {
                uint ii = GetLatIndex(45.0);
                double[] Temp1 = new double[_NumLonCells];
                double[] Temp2 = new double[_NumLonCells];

                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    Temp1[jj] = _GlobalMissingValue;
                    Temp2[jj] = _GlobalMissingValue;
                }

                //for this latitude band calculate the annual average conditions across realms
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];

                    if (globalRealm[ii, jj] == 1.0)
                    {
                        var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                        Temp1[jj] = subset.Sum() / subset.Count();

                        if (double.IsNaN(Temp1[jj]) ||
                            double.IsInfinity(Temp1[jj]))
                            Temp1[jj] = _GlobalMissingValue;
                    }
                    else
                    {
                        var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                        Temp2[jj] = subset.Sum() / subset.Count();

                        if (double.IsNaN(Temp2[jj]) ||
                            double.IsInfinity(Temp2[jj]))
                            Temp2[jj] = _GlobalMissingValue;
                    }

                }

                var subset1 = from i in Temp1.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal1 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal1) || double.IsInfinity(AverageVal1))
                    AverageVal1 = GlobalMissingValue;

                subset1 = from i in Temp2.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal2 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal2) || double.IsInfinity(AverageVal2))
                    AverageVal2 = GlobalMissingValue;


            }
            // Global average conditions
            else if (spatialDomain == 1)
            {
                double[,] Temp1 = new double[_NumLatCells, _NumLonCells];
                double[,] Temp2 = new double[_NumLatCells, _NumLonCells];

                for (int ii = 0; ii < _NumLatCells; ii++)
                {
                    for (int jj = 0; jj < _NumLonCells; jj++)
                    {
                        Temp1[ii, jj] = _GlobalMissingValue;
                        Temp2[ii, jj] = _GlobalMissingValue;
                    }
                }

                for (int ii = 0; ii < _NumLatCells; ii++)
                {
                    for (int jj = 0; jj < _NumLonCells; jj++)
                    {
                        TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];

                        if (globalRealm[ii, jj] == 1.0)
                        {
                            var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                            Temp1[ii, jj] = subset.Sum() / subset.Count();

                            if (double.IsNaN(Temp1[ii, jj]) ||
                                double.IsInfinity(Temp1[ii, jj]))
                                Temp1[ii, jj] = _GlobalMissingValue;
                        }
                        else
                        {
                            var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                            Temp2[ii, jj] = subset.Sum() / subset.Count();

                            if (double.IsNaN(Temp2[ii, jj]) ||
                                double.IsInfinity(Temp2[ii, jj]))
                                Temp2[ii, jj] = _GlobalMissingValue;
                        }

                    }

                }
                var subset1 = from i in Temp1.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal1 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal1) || double.IsInfinity(AverageVal1))
                    AverageVal1 = GlobalMissingValue;

                subset1 = from i in Temp2.Cast<double>().ToArray() where (i != _GlobalMissingValue) select i;
                AverageVal2 = subset1.Sum() / subset1.Count();
                if (double.IsNaN(AverageVal2) || double.IsInfinity(AverageVal2))
                    AverageVal2 = GlobalMissingValue;
            }

            return new double[] {AverageVal1,AverageVal2};
        }

        /// <summary>
        /// Convert a layer into meridional average conditions
        /// </summary>
        /// <param name="annual">Calculate the annual, meridional average</param>
        /// <param name="layerName">the name of the layer to convert</param>
        public void ConvertToMeridionalAverage(bool annual,string layerName)
        {
            int NumTimsteps = GetCellEnvironment(1, 1)[layerName].Length;
            double[,] AnnualMnArray = new double[_NumLatCells, _NumLonCells];
            double[] TempEnv;
            double[] Temp1 = new double[_NumLonCells];
            double[] Temp2 = new double[_NumLonCells];
            double[] MeridionalRealm = new double[_NumLonCells];
            


            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                double[] MeridionalAv1 = new double[NumTimsteps];
                double[] MeridionalAv2 = new double[NumTimsteps];

                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    MeridionalRealm[jj] = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
                    Temp1[jj] = _GlobalMissingValue;
                    Temp2[jj] = _GlobalMissingValue;
                }

                //Check if annual average is required
                if (annual)
                {
                    for (int jj = 0; jj < _NumLonCells; jj++)
                    {
                        TempEnv = InternalGrid[ii, jj].CellEnvironment[layerName];
                        if(MeridionalRealm[jj] == 1)
                        {
                            var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                            Temp1[jj] = subset.Sum() / subset.Count();

                            if(double.IsNaN(Temp1[jj]) || double.IsInfinity(Temp1[jj])) 
                            Temp1[jj] = _GlobalMissingValue;
                        }
                        else
                        {
                            var subset = from i in TempEnv where i != _GlobalMissingValue select i;
                            Temp2[jj] = subset.Sum() / subset.Count();

                            if(double.IsNaN(Temp2[jj]) || double.IsInfinity(Temp2[jj])) 
                            Temp2[jj] = _GlobalMissingValue;
                        }

                    }
                    // Calculate the meridional average using the annual mean values
                    for (int n = 0; n < NumTimsteps; n++)
			        {
                        var subset = from i in Temp1 where (i != _GlobalMissingValue) select i;
			            MeridionalAv1[n] = subset.Sum() / subset.Count();
                        if (double.IsNaN(MeridionalAv1[n]) || double.IsInfinity(MeridionalAv1[n]))
                            MeridionalAv1[n] = GlobalMissingValue;

                        subset = from i in Temp2 where (i != _GlobalMissingValue) select i;
			            MeridionalAv2[n] = subset.Sum() / subset.Count();
                        if (double.IsNaN(MeridionalAv2[n]) || double.IsInfinity(MeridionalAv2[n]))
                            MeridionalAv2[n] = GlobalMissingValue;
			        }        
                    
                }
                else
                {
                    // Calculate the meridional average for each timestep
                    for (int n = 0; n < NumTimsteps; n++)
                    {
                        for (int jj = 0; jj < _NumLonCells; jj++)
                        {
                            if(MeridionalRealm[jj] == 1.0)
                            {
                                Temp1[jj] = InternalGrid[ii, jj].CellEnvironment[layerName][n] ;
                            }
                            else
                            {
                                Temp2[jj] = InternalGrid[ii, jj].CellEnvironment[layerName][n] ;
                            }
                        }

                        var subset = from i in Temp1 where i != _GlobalMissingValue select i;
                        MeridionalAv1[n] = subset.Sum() / subset.Count();
                        if (double.IsNaN(MeridionalAv1[n]) || double.IsInfinity(MeridionalAv1[n]))
                            MeridionalAv1[n] = GlobalMissingValue;

                        subset = from i in Temp2 where i != _GlobalMissingValue select i;
                        MeridionalAv2[n] = subset.Sum() / subset.Count();
                        if (double.IsNaN(MeridionalAv2[n]) || double.IsInfinity(MeridionalAv2[n]))
                            MeridionalAv2[n] = GlobalMissingValue;
                    }
                }
                //Loop over all longitudes and set the cell environment equal to the meridional average conditions
                for (int jj = 0; jj < NumLonCells; jj++)
                {
                    if(MeridionalRealm[jj] == 1.0)
                    {
                        InternalGrid[ii, jj].CellEnvironment[layerName] = MeridionalAv1;
                    }
                    else
                    {
                        InternalGrid[ii, jj].CellEnvironment[layerName] = MeridionalAv2;
                    }
                    
                }
                

            }

            
        }

        /// <summary>
        /// Seed the stocks and cohorts from output from a previous simulation
        /// </summary>
        /// <param name="cellIndices">A list of the active cells in the model grid</param>
        /// <param name="cohortFunctionalGroupDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="stockFunctionalGroupDefinitions">The functional group definitions for stocks in the model</param>
        /// <param name="globalDiagnostics">A list of global diagnostic variables</param>
        /// <param name="nextCohortID">The ID number to be assigned to the next produced cohort</param>
        /// <param name="tracking">Whether process-tracking is enabled</param>
        /// <param name="DrawRandomly">Whether the model is set to use a random draw</param>
        /// <param name="dispersalOnly">Whether to run dispersal only (i.e. to turn off all other ecological processes</param>
        /// <param name="processTrackers">An instance of the ecological process tracker</param>
        public void SeedGridCellStocksAndCohorts(List<uint[]> cellIndices,
            InputModelState inputModelState,
            FunctionalGroupDefinitions cohortFunctionalGroupDefinitions,
            FunctionalGroupDefinitions stockFunctionalGroupDefinitions)
        //            FunctionalGroupDefinitions stockFunctionalGroupDefinitions, SortedList<string, double> globalDiagnostics, ref Int64 nextCohortID,
        //            Boolean tracking, Boolean DrawRandomly, Boolean dispersalOnly, string dispersalOnlyType)
        {
            int ii = 1;
            Console.WriteLine("Seeding grid cell stocks and cohorts:");


            //Check to see if the correct number of functional groups exist in the definitions file and in the input state

            if (cohortFunctionalGroupDefinitions.GetNumberOfFunctionalGroups() != inputModelState.GridCellCohorts[0, 0].Count)
            {
                Console.WriteLine("Mismatch in the number of functional groups defined in CohortFunctionalGroupDefinitions.csv set-up file and the Model State being read in");
                Environment.Exit(0);
            }

            int[] TerrestrialStockFunctionalIndices = stockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Terrestrial", false);
            int[] MarineStockFunctionalIndices = stockFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Marine", false);

            int[] TerrestrialCohortFunctionalIndices = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Terrestrial", false);
            int[] MarineCohortFunctionalIndices = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Marine", false);


            foreach (uint[] cellIndexPair in cellIndices)
            {

                for (int i = 0; i < inputModelState.GridCellCohorts[cellIndexPair[0], cellIndexPair[1]].Count; i++)
                {
                    InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellCohorts[i] = new List<Cohort>();
                }
                //Check which cohorts should be initialised for each cell
                if (InternalGrid[cellIndexPair[0], cellIndexPair[1]].CellEnvironment["Realm"][0] == 1)
                {
                    //This is a terrestrial cell so only add the terrestrial stocks
                    foreach (int fg in TerrestrialCohortFunctionalIndices)
                    {
                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellCohorts[fg] = inputModelState.GridCellCohorts[cellIndexPair[0], cellIndexPair[1]][fg];
                    }
                }
                else
                {
                    // this is a marine cell so only add the marine stocks
                    foreach (int fg in MarineCohortFunctionalIndices)
                    {
                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellCohorts[fg] = inputModelState.GridCellCohorts[cellIndexPair[0], cellIndexPair[1]][fg];
                    }
                }



                for (int i = 0; i < inputModelState.GridCellStocks[cellIndexPair[0], cellIndexPair[1]].Count; i++)
                {
                    InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellStocks[i] = new List<Stock>();
                }


                //Check which stocks should be initialised for each cell
                if (InternalGrid[cellIndexPair[0], cellIndexPair[1]].CellEnvironment["Realm"][0] == 1)
                {
                    //This is a terrestrial cell so only add the terrestrial stocks
                    foreach (int fg in TerrestrialStockFunctionalIndices)
                    {
                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellStocks[fg] = inputModelState.GridCellStocks[cellIndexPair[0], cellIndexPair[1]][fg];
                    }
                }
                else
                {
                    // this is a marine cell so only add the marine stocks
                    foreach (int fg in MarineStockFunctionalIndices)
                    {
                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].GridCellStocks[fg] = inputModelState.GridCellStocks[cellIndexPair[0], cellIndexPair[1]][fg];
                    }
                }

                Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);

            }

            Console.WriteLine("");
            Console.WriteLine("");
        }



        /// <summary>
        /// Seed the stocks and cohorts for all active cells in the model grid
        /// </summary>
        /// <param name="cellIndices">A list of the active cells in the model grid</param>
        /// <param name="cohortFunctionalGroupDefinitions">The functional group definitions for cohorts in the model</param>
        /// <param name="stockFunctionalGroupDefinitions">The functional group definitions for stocks in the model</param>
        /// <param name="globalDiagnostics">A list of global diagnostic variables</param>
        /// <param name="nextCohortID">The ID number to be assigned to the next produced cohort</param>
        /// <param name="tracking">Whether process-tracking is enabled</param>
        /// <param name="DrawRandomly">Whether the model is set to use a random draw</param>
        /// <param name="dispersalOnly">Whether to run dispersal only (i.e. to turn off all other ecological processes</param>
        /// <param name="dispersalOnlyType">For dispersal only runs, the type of dispersal to apply</param>
        public void SeedGridCellStocksAndCohorts(List<uint[]> cellIndices, FunctionalGroupDefinitions cohortFunctionalGroupDefinitions, 
            FunctionalGroupDefinitions stockFunctionalGroupDefinitions, SortedList<string,double> globalDiagnostics, ref Int64 nextCohortID,
            Boolean tracking, Boolean DrawRandomly, Boolean dispersalOnly, string dispersalOnlyType)
        {
            int ii = 1;
            Console.WriteLine("Seeding grid cell stocks and cohorts:");

            //Work out how many cohorts are to be seeded in each grid cell - split by realm as different set of cohorts initialised by realm
            double TotalTerrestrialCellCohorts = 0;
            double TotalMarineCellCohorts = 0;

            int[] TerrestrialFunctionalGroups = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Terrestrial", false);
            if (TerrestrialFunctionalGroups == null)
            {
                TotalTerrestrialCellCohorts = 0;
            }
            else
            {
                foreach (int F in TerrestrialFunctionalGroups)
                {
                    TotalTerrestrialCellCohorts += cohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup("Initial number of GridCellCohorts", F);
                }
            }


            int[] MarineFunctionalGroups = cohortFunctionalGroupDefinitions.GetFunctionalGroupIndex("Realm", "Marine", false);
            if (MarineFunctionalGroups == null)
            {
                TotalMarineCellCohorts = 0;
            }
            else
            {
                foreach (int F in MarineFunctionalGroups)
                {
                    TotalMarineCellCohorts += cohortFunctionalGroupDefinitions.GetBiologicalPropertyOneFunctionalGroup("Initial number of GridCellCohorts", F);
                }
            }

            foreach (uint[] cellIndexPair in cellIndices)
            {
                if (dispersalOnly)
                {
                    if (dispersalOnlyType == "diffusion")
                    {
                        // Diffusive dispersal

                        if ((cellIndexPair[0] == 90) && (cellIndexPair[1] == 180))
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, false);
                        }
                        else if ((cellIndexPair[0] == 95) && (cellIndexPair[1] == 110))
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, false);
                        }
                        else
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, true);
                        }
                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
                    }
                    else if (dispersalOnlyType == "advection")
                    {
                        // Advective dispersal
                        /*
                        if ((cellIndexPair[0] == 58) && (cellIndexPair[1] == 225))
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, false);
                        }
                        else if ((cellIndexPair[0] == 95) && (cellIndexPair[1] == 110))
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, false);
                        }
                        else
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                            stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                            DrawRandomly, true);
                        }
                        */
                        if (InternalGrid[cellIndexPair[0], cellIndexPair[1]].CellEnvironment["Realm"][0] == 1.0)
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(
                                cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, globalDiagnostics, 
                                ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts, 
                                DrawRandomly, true);
                        }
                        else
                        {
                            InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(
                                cohortFunctionalGroupDefinitions,stockFunctionalGroupDefinitions, globalDiagnostics, 
                                ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts, 
                                DrawRandomly, false);
                        }
                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
                    }
                    else if (dispersalOnlyType == "responsive")
                    {
                        // Responsive dispersal

                        InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(cohortFunctionalGroupDefinitions,
                        stockFunctionalGroupDefinitions, globalDiagnostics, ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                        DrawRandomly, true);

                        Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
                    }
                    else
                    {
                        Debug.Fail("Dispersal only type not recognized from initialisation file");
                    }
                }

                else
                {
                    InternalGrid[cellIndexPair[0], cellIndexPair[1]].SeedGridCellCohortsAndStocks(
                        cohortFunctionalGroupDefinitions, stockFunctionalGroupDefinitions, globalDiagnostics, 
                        ref nextCohortID, tracking, TotalTerrestrialCellCohorts, TotalMarineCellCohorts,
                        DrawRandomly, false);
                    Console.Write("\rGrid Cell: {0} of {1}", ii++, cellIndices.Count);
                }
            }
            Console.WriteLine("");
            Console.WriteLine("");
        }

        /// <summary>
        /// Returns the stocks within the specified grid cell
        /// </summary>
        /// <param name="latIndex">Latitude index</param>
        /// <param name="lonIndex">Longitude index</param>
        /// <returns>The stock handler for the specified grid cell</returns>
        public GridCellStockHandler GetGridCellStocks(uint latIndex, uint lonIndex)
        {
            return InternalGrid[latIndex, lonIndex].GridCellStocks;
        }

        /// <summary>
        /// Sets the stocks in the specified grid cell to the passed stocks
        /// </summary>
        /// <param name="newGridCellStocks">New stocks for the grid cell</param>
        /// <param name="latIndex">Latitude index</param>
        /// <param name="lonIndex">Longitude index</param>
        public void SetGridCellStocks(GridCellStockHandler newGridCellStocks, uint latIndex, uint lonIndex)
        {
            InternalGrid[latIndex, lonIndex].GridCellStocks = newGridCellStocks;
        }

        /// <summary>
        /// Returns the array (indexed by functional group) of lists of gridCellCohorts for the specified grid cell
        /// </summary>
        /// <param name="latIndex">Latitude index of grid cell</param>
        /// <param name="lonIndex">Longitude index of grid cell</param>
        /// <returns>Arry (indexed by functional group) of lists of gridCellCohorts</returns>
        public GridCellCohortHandler GetGridCellCohorts(uint latIndex, uint lonIndex)
        {
            return InternalGrid[latIndex, lonIndex].GridCellCohorts;
        }

        /// <summary>
        /// Extracts an individual cohort from a particular grid cell
        /// </summary>
        /// <param name="latIndex">Latitude index of grid cell</param>
        /// <param name="lonIndex">Longitude index of grid cell</param>
        /// <param name="functionalGroup">Functional group of cohort</param>
        /// <param name="positionInList">Index of cohort position in the list</param>
        /// <returns></returns>
        public Cohort GetGridCellIndividualCohort(uint latIndex, uint lonIndex, int functionalGroup, int positionInList)
        {
            return InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].ElementAt(positionInList);
        }

        // NOTE TO SELF: These need more error checking, and also the access levels more tightly controlled
        /// <summary>
        /// Remove an individual cohort from a functionall group; necessary due to dispersal moving cohorts from one cell to another
        /// </summary>
        /// <param name="latIndex">Grid cell latitude index</param>
        /// <param name="lonIndex">Grid cell longitude index</param>
        /// <param name="functionalGroup">Cohort functional group</param>
        /// <param name="positionInList">Position of cohort in the list of that functional group</param>
        public void DeleteGridCellIndividualCohort(uint latIndex, uint lonIndex, int functionalGroup, int positionInList)
        {
            InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].RemoveAt(positionInList);
        }

        /// <summary>
        /// Delete a specified list of cohorts from a grid cell
        /// </summary>
        /// <param name="latIndex">The latitudinal index of the grid cell to delete cohorts from</param>
        /// <param name="lonIndex">The longitudinal index of the grid cell to delete cohorts from</param>
        /// <param name="cohortFGsToDelete">A list of the functional groups that each cohort to delete belongs to</param>
        /// <param name="cohortNumbersToDelete">A list of the positions with each functional group that each cohort to delete occupies</param>
        /// <remarks>This is inefficient and needs double-checking for errors</remarks>
        public void DeleteGridCellIndividualCohorts(uint latIndex, uint lonIndex, List<uint> cohortFGsToDelete, List<uint> cohortNumbersToDelete)
        {
            
            // Get the unique functional groups that have cohorts to be removed
            uint[] TempList = cohortFGsToDelete.Distinct().ToArray();

            // Loop over these unique functional  groups
            for (int ii = 0; ii < TempList.Length; ii++)
			{
                // Get the functional group index of the current functional group
                int FG = (int)TempList[ii];

                // Create a local list to hold the positions of the cohorts to delete from this functional group
			    List<uint> CohortIndexList = new List<uint>();
                // Loop over all cohorts to be deleted
                for (int jj = 0; jj < cohortFGsToDelete.Count; jj++)
			    {
                    // Check whether the functional group correpsonds with the functional group currently being processed 
                    if (cohortFGsToDelete.ElementAt((int)jj) == FG)
                    {
                        // Add the cohort to the list of cohorts to delete
                        CohortIndexList.Add(cohortNumbersToDelete[jj]);
                    }
			    }

                // Sort the list of positions of the cohorts to delete in this functional group
                CohortIndexList.Sort();
                // Reverse the list so that the highest positions come first
                CohortIndexList.Reverse();

                // Loop over cohorts and delete in turn, starting with cohorts in the highest positions
                for (int kk = 0; kk < CohortIndexList.Count; kk++)
			    {
			        InternalGrid[latIndex, lonIndex].GridCellCohorts[FG].RemoveAt((int)CohortIndexList[kk]);                  
			    }

	        }
            
        }

        /// <summary>
        /// Replace the gridCellCohorts in a grid cell with a new list of gridCellCohorts
        /// </summary>
        /// <param name="newGridCellCohorts">The new list of gridCellCohorts</param>
        /// <param name="latIndex">Grid cell latitude index</param>
        /// <param name="lonIndex">Grid cell longitude index</param>
        public void SetGridCellCohorts(GridCellCohortHandler newGridCellCohorts, uint latIndex, uint lonIndex)
        {
            InternalGrid[latIndex, lonIndex].GridCellCohorts = newGridCellCohorts;
        }

        /// <summary>
        /// Add a new cohort to an existing list of cohorts in the grid cell - or create a new list if there is not one present
        /// </summary>
        /// <param name="latIndex">Latitude index of the grid cell</param>
        /// <param name="lonIndex">Longitude index of the grid cell</param>
        /// <param name="functionalGroup">Functional group of the cohort (i.e. array index)</param>
        /// <param name="cohortToAdd">The cohort object to add</param>
        public void AddNewCohortToGridCell(uint latIndex, uint lonIndex, int functionalGroup, Cohort cohortToAdd)
        {
            InternalGrid[latIndex, lonIndex].GridCellCohorts[functionalGroup].Add(cohortToAdd);
        }

        /// <summary>
        /// Return the value of a specified environmental layer from an individual grid cell
        /// </summary>
        /// <param name="variableName">The name of the environmental lyaer</param>
        /// <param name="timeInterval">The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <param name="latCellIndex">The latitudinal cell index</param>
        /// <param name="lonCellIndex">The longitudinal cell index</param>
        /// <param name="variableExists">Returns false if the environmental layer does not exist, true if it does</param>
        /// <returns>The value of the environmental layer, or a missing value if the environmental layer does not exist</returns>
        public double GetEnviroLayer(string variableName, uint timeInterval, uint latCellIndex, uint lonCellIndex, out bool variableExists)
        {
            return InternalGrid[latCellIndex, lonCellIndex].GetEnviroLayer(variableName,timeInterval, out variableExists);
        }

        /// <summary>
        /// Set the value of a specified environmental layer in an individual grid cell
        /// </summary>
        /// <param name="variableName">The name of the environmental layer</param>
        /// <param name="timeInterval">The time interval within the environmental variable to set (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <param name="setValue">The value to set</param>
        /// <param name="latCellIndex">The latitudinal cell index</param>
        /// <param name="lonCellIndex">The longitudinal cell index</param>
        /// <returns>True if the value is set successfully, false otherwise</returns>
        public bool SetEnviroLayer(string variableName, uint timeInterval, double setValue, uint latCellIndex, uint lonCellIndex)
        {
            return InternalGrid[latCellIndex, lonCellIndex].SetEnviroLayer(variableName,timeInterval, setValue);
        }

        
        /// <summary>
        /// Set the value of a given delta type for the specified ecological process within the specified grid cell
        /// </summary>
        /// <param name="deltaType">The type of delta value to set (e.g. 'biomass', 'abundance' etc.)</param>
        /// <param name="ecologicalProcess">The name of the ecological process to set the value of delta for</param>
        /// <param name="setValue">The value to set</param>
        /// <param name="latCellIndex">The latitudinal index of the cell</param>
        /// <param name="lonCellIndex">The longitudinal index of the cell</param>
        /// <returns>True if the value is set successfully, false otherwise</returns>
        public bool SetDeltas(string deltaType, string ecologicalProcess, double setValue, uint latCellIndex, uint lonCellIndex)
        {
            return InternalGrid[latCellIndex, lonCellIndex].SetDelta(deltaType, ecologicalProcess, setValue);
        }
        
        /// <summary>
        /// Get the total of a state variable for specific cells
        /// </summary>
        /// <param name="variableName">The name of the variable</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model intialisation</param>
        /// <returns>Summed value of variable over whole grid</returns>
        /// <todo>Overload to work with vector and array state variables</todo>
        public double StateVariableGridTotal(string variableName, string traitValue, int[] functionalGroups, List<uint[]> cellIndices, 
            string stateVariableType, MadingleyModelInitialisation initialisation)
        {

            double tempVal = 0;

            double[,] TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

            // Loop through and sum values across a grid, excluding missing values
            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                tempVal += TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]];
            }

            return tempVal;
        }

        /// <summary>
        /// Gets a state variable for specified functional groups of specified entity types in a specified grid cell
        /// </summary>
        /// <param name="variableName">The name of the variable to get: 'biomass' or 'abundance'</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">The functional group indices to get the state variable for</param>
        /// <param name="latCellIndex">The latitudinal index of the cell</param>
        /// <param name="lonCellIndex">The longitudinal index of the cell</param>
        /// <param name="stateVariableType">The type of entity to return the state variable for: 'stock' or 'cohort'</param>
        /// <param name="modelInitialisation">The Madingley Model initialisation</param>
        /// <returns>The state variable for specified functional groups of specified entity types in a specified grid cell</returns>
        public double GetStateVariable(string variableName, string traitValue, int[] functionalGroups, uint latCellIndex, uint lonCellIndex, 
            string stateVariableType, MadingleyModelInitialisation modelInitialisation)
        {

            double returnValue = 0.0;

            switch (stateVariableType.ToLower())
            {
                case "cohort":
                    
                    GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex, lonCellIndex].GridCellCohorts;

                    switch (variableName.ToLower())
                    {
                        case "biomass":
                            if (traitValue != "Zooplankton")
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                    }
                                }
                            }
                            else
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                    }
                                }
                            }
                            break;

                        case "abundance":
                            if (traitValue != "Zooplankton")
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        returnValue += item.CohortAbundance;
                                    }
                                }
                            }
                            else
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                        returnValue += item.CohortAbundance;
                                    }
                                }
                            }
                            break;

                        default:
                            Debug.Fail("For cohorts, state variable name must be either 'biomass' or 'abundance'");
                            break;
                    }
                    break;

                case "stock":
                    GridCellStockHandler TempStocks = InternalGrid[latCellIndex, lonCellIndex].GridCellStocks;

                    switch (variableName.ToLower())
                    {
                        case "biomass":
                            foreach (int f in functionalGroups)
                            {
                                foreach (var item in TempStocks[f])
                                {
                                    returnValue += item.TotalBiomass;
                                }
                            }
                            break;
                        default:
                            Debug.Fail("For stocks, state variable name must be 'biomass'");
                            break;
                    }
                    break;

                default:
                    Debug.Fail("State variable type must be either 'cohort' or 'stock'");
                    break;

            }

            





            return returnValue;
        }

        /// <summary>
        /// Gets a state variable density for specified functional groups of specified entity types in a specified grid cell
        /// </summary>
        /// <param name="variableName">The name of the variable to get: 'biomass' or 'abundance'</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">The functional group indices to get the state variable for</param>
        /// <param name="latCellIndex">The latitudinal index of the cell</param>
        /// <param name="lonCellIndex">The longitudinal index of the cell</param>
        /// <param name="stateVariableType">The type of entity to return the state variable for: 'stock' or 'cohort'</param>
        /// <param name="modelInitialisation">The Madingley Model initialisation</param>
        /// <returns>The state variable density for specified functional groups of specified entity types in a specified grid cell</returns>
        public double GetStateVariableDensity(string variableName, string traitValue, int[] functionalGroups, uint latCellIndex, 
            uint lonCellIndex, string stateVariableType, MadingleyModelInitialisation modelInitialisation)
        {

            double returnValue = 0.0;

            switch (stateVariableType.ToLower())
            {
                case "cohort":

                    GridCellCohortHandler TempCohorts = InternalGrid[latCellIndex, lonCellIndex].GridCellCohorts;

                    switch (variableName.ToLower())
                    {
                        case "biomass":
                            if (traitValue != "Zooplankton (all)")
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                    }
                                }
                            }
                            else
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                            returnValue += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                    }
                                }
                            }
                            break;

                        case "abundance":
                            if (traitValue != "Zooplankton (all)")
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        returnValue += item.CohortAbundance;
                                    }
                                }
                            }
                            else
                            {
                                foreach (int f in functionalGroups)
                                {
                                    foreach (var item in TempCohorts[f])
                                    {
                                        if (item.IndividualBodyMass <= modelInitialisation.PlanktonDispersalThreshold)
                                            returnValue += item.CohortAbundance;
                                    }
                                }
                            }
                            break;

                        default:
                            Debug.Fail("For cohorts, state variable name must be either 'biomass' or 'abundance'");
                            break;
                    }
                    break;

                case "stock":
                    GridCellStockHandler TempStocks = InternalGrid[latCellIndex, lonCellIndex].GridCellStocks;

                    switch (variableName.ToLower())
                    {
                        case "biomass":
                            foreach (int f in functionalGroups)
                            {
                                foreach (var item in TempStocks[f])
                                {
                                    returnValue += item.TotalBiomass;
                                }
                            }
                            break;
                        default:
                            Debug.Fail("For stocks, state variable name must be 'biomass'");
                            break;
                    }
                    break;

                default:
                    Debug.Fail("State variable type must be either 'cohort' or 'stock'");
                    break;

            }


            return returnValue / (InternalGrid[latCellIndex, lonCellIndex].CellEnvironment["Cell Area"][0]);
        }

        /// <summary>
        /// Get the mean density of a state variable for specific cells
        /// </summary>
        /// <param name="variableName">The name of the variable</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model initialisation</param>
        /// <returns>Mean density of variable over whole grid</returns>
        public double StateVariableGridMeanDensity(string variableName, string traitValue, int[] functionalGroups, List<uint[]> cellIndices, 
            string stateVariableType, MadingleyModelInitialisation initialisation)
        {

            double tempVal = 0;

            double[,] TempStateVariable = this.GetStateVariableGridDensityPerSqKm(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

            // Loop through and sum values across a grid, excluding missing values
            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                tempVal += TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]];
            }

            return tempVal / cellIndices.Count;
        }


        /// <summary>
        /// Return an array of values for a single state variable over specific cells
        /// </summary>
        /// <param name="variableName">Variable name</param>
        /// <param name="traitValue">The trait values of functional groups to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model initialisation</param>
        /// <returns>Array of state variable values for each grid cell</returns>
        public double[,] GetStateVariableGrid(string variableName, string traitValue, int[] functionalGroups, List<uint[]> cellIndices, 
            string stateVariableType, MadingleyModelInitialisation initialisation)
        {
            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];

            switch (variableName.ToLower())
            {
                case "biomass":
                    for (int ii = 0; ii < cellIndices.Count; ii++)
                    {
                        // Check whether the state variable concerns cohorts or stocks
                        if (stateVariableType.ToLower() == "cohort")
                        {
                            if (traitValue != "Zooplankton")
                            {
                                // Check to make sure that the cell has at least one cohort
                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                                {
                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
                                    {
                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                        {
                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
                                            {
                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                // Check to make sure that the cell has at least one cohort
                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                                {
                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
                                    {
                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                        {
                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
                                            {
                                                if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += ((item.IndividualBodyMass + item.IndividualReproductivePotentialMass) * item.CohortAbundance);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else if (stateVariableType.ToLower() == "stock")
                        {
                            // Check to make sure that the cell has at least one stock
                            if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks != null)
                            {
                                for (int nn = 0; nn < functionalGroups.Length; nn++)
                                {
                                    if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]] != null)
                                    {
                                        foreach (Stock item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellStocks[functionalGroups[nn]].ToArray())
                                        {
                                            TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += (item.TotalBiomass);

                                        }
                                    }

                                }
                            }
                        }
                        else
                        {
                            Debug.Fail("Variable 'state variable type' must be either 'stock' 'or 'cohort'");
                        }
                        
                    }
                    break;
                case "abundance":
                    for (int ii = 0; ii < cellIndices.Count; ii++)
                    {
                        // Check whether the state variable concerns cohorts or stocks
                        if (stateVariableType.ToLower() == "cohort")
                        {
                            if (traitValue != "Zooplankton")
                            {
                                // Check to make sure that the cell has at least one cohort
                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                                {
                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
                                    {
                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                        {
                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
                                            {
                                                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += item.CohortAbundance;
                                            }
                                        }
                                    }
                                }
                            }
                            else
                            {
                                // Check to make sure that the cell has at least one cohort
                                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts != null)
                                {
                                    for (int nn = 0; nn < functionalGroups.Length; nn++)
                                    {
                                        if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]] != null)
                                        {
                                            foreach (Cohort item in InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].GridCellCohorts[functionalGroups[nn]].ToArray())
                                            {
                                                if (item.IndividualBodyMass <= initialisation.PlanktonDispersalThreshold)
                                                    TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] += item.CohortAbundance;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        else
                        {
                            Debug.Fail("Currently abundance cannot be calculated for grid cell stocks");
                        }
                    }
                    break;
                default:
                    Debug.Fail("Invalid search string passed for cohort property");
                    break;
            }

            return TempStateVariable;

        }
        
        /// <summary>
        /// Return an array of values for a single state variable over specific cells, given in densities per km^2
        /// </summary>
        /// <param name="variableName">Variable name</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model initialisation</param>
        /// <returns>Array of state variable values for each grid cell</returns>
        public double[,] GetStateVariableGridDensityPerSqKm(string variableName, string traitValue, int[] functionalGroups, 
            List<uint[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
        {
            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
            double CellArea;

            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] /= CellArea;
            }

            return TempStateVariable;
        }


        /// <summary>
        /// Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells
        /// </summary>
        /// <param name="variableName">The name of the variable</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model initialisation</param>
        /// <returns>Array of log(state variable values +1 ) for each grid cell</returns>
        public double[,] GetStateVariableGridLog(string variableName, string traitValue, int[] functionalGroups, List<uint[]> cellIndices, 
            string stateVariableType, MadingleyModelInitialisation initialisation)
        {

            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];

            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);
            
            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] = Math.Log(TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]]+1);
            }

            return TempStateVariable;
        }


        /// <summary>
        /// Return an array of log(values + 1) for a state variable for particular functional groups over specific cells. State variable (currently only biomass or abundance) must be >= 0 in all grid cells
        /// </summary>
        /// <param name="variableName">The name of the variable</param>
        /// <param name="traitValue">The functional group trait value to get data for</param>
        /// <param name="functionalGroups">A vector of functional group indices to consider</param>
        /// <param name="cellIndices">List of indices of active cells in the model grid</param>
        /// <param name="stateVariableType">A string indicating the type of state variable; 'cohort' or 'stock'</param>
        /// <param name="initialisation">The Madingley Model intialisation</param>
        /// <returns>Array of log(state variable values +1 ) for each grid cell</returns>
        public double[,] GetStateVariableGridLogDensityPerSqKm(string variableName, string traitValue, int[] functionalGroups, 
            List<uint[]> cellIndices, string stateVariableType, MadingleyModelInitialisation initialisation)
        {

            double[,] TempStateVariable = new double[this.NumLatCells, this.NumLonCells];
            double CellArea;

            TempStateVariable = this.GetStateVariableGrid(variableName, traitValue, functionalGroups, cellIndices, stateVariableType, initialisation);

            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                CellArea = GetCellEnvironment(cellIndices[ii][0], cellIndices[ii][1])["Cell Area"][0];
                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] /= CellArea;
                TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]] = Math.Log(TempStateVariable[cellIndices[ii][0], cellIndices[ii][1]]+1);
            }

            return TempStateVariable;

        }

        /// <summary>
        /// Returns, for a given longitude, the appropriate longitude index in the grid
        /// ASSUMES THAT LONGITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!
        /// </summary>
        /// <param name="myLon">Longitude, in degrees</param>
        /// <returns>longitude index in the model grid</returns>
        public uint GetLonIndex(double myLon)
        {
            Debug.Assert((myLon >= _MinLongitude && myLon < _MaxLongitude), "Error: latitude out of range");

            return (uint)Math.Floor((myLon - _MinLongitude) / _LonCellSize);

        }

        /// <summary>
        /// Return the longitude of a cell at a particular lon. index
        /// </summary>
        /// <param name="cellLonIndex">The longitudinal index (i.e. row) of the cell</param>
        /// <returns>Returns the longitude of the bottom of the cell, in degrees</returns>
        public double GetCellLongitude(uint cellLonIndex)
        {
            Debug.Assert((cellLonIndex <= (_NumLonCells - 1)), "Error: Cell index out of range when trying to find the longitude for a particular cell");

            double TempLongitude = double.MaxValue;

            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                if (InternalGrid[ii, cellLonIndex] != null)
                    TempLongitude = InternalGrid[ii, cellLonIndex].Longitude;
            }

            Debug.Assert(TempLongitude != double.MaxValue, "Error trying to find cell longitude - no grid cells have been initialised for this latitude index: " + cellLonIndex.ToString());

            return TempLongitude;
        }

        /// <summary>
        /// Return the latitude of a cell at a particular lat. index
        /// </summary>
        /// <param name="cellLatIndex">The latitudinal index (i.e. row) of the cell</param>
        /// <returns>Returns the latitude of the bottom of the cell, in degrees</returns>
        public double GetCellLatitude(uint cellLatIndex)
        {
            Debug.Assert((cellLatIndex <= (_NumLatCells - 1)), "Error: Cell index out of range when trying to find the latitude for a particular cell");

            double TempLatitude = double.MaxValue;

            for (int jj = 0; jj < _NumLonCells ; jj++)
            {
                if (InternalGrid[cellLatIndex, jj] != null)
                {
                    TempLatitude = InternalGrid[cellLatIndex, jj].Latitude;
                    break;
                }
            }

            Debug.Assert(TempLatitude != double.MaxValue, "Error trying to find cell latitude - no grid cells have been initialised for this latitude index: " + cellLatIndex.ToString());

            return TempLatitude;

        }


        /// <summary>
        /// Returns, for a given latitude, the appropriate latitude index in the grid
        /// ASSUMES THAT LATITUDES IN THE MODEL GRID OBJECT REFER TO LOWER LEFT CORNERS!!!
        /// </summary>
        /// <param name="myLat">Latitude, in degrees</param>
        /// <returns>latitude index in the model grid</returns>
        public uint GetLatIndex(double myLat)
        {

            Debug.Assert((myLat >= _MinLatitude && myLat < _MaxLatitude), "Error: latitude out of range");

            return (uint)Math.Floor((myLat - _MinLatitude) / _LatCellSize);

        }

        /// <summary>
        /// A method to return the values for all environmental data layers for a particular grid cell
        /// </summary>
        /// <param name="cellLatIndex">Latitude index of grid cell</param>
        /// <param name="cellLonIndex">Longitude index of grid cell</param>
        /// <returns>A sorted list containing environmental data layer names and values</returns>
        public SortedList<string, double[]> GetCellEnvironment(uint cellLatIndex, uint cellLonIndex)
        {
            return InternalGrid[cellLatIndex, cellLonIndex].CellEnvironment;
        }

        /// <summary>
        /// A method to return delta values for the specified delta type in a particular grid cell
        /// </summary>
        /// <param name="deltaType">The delta type to return</param>
        /// <param name="cellLatIndex">Latitude index of grid cell</param>
        /// <param name="cellLonIndex">Longitude index of grid cell</param>
        /// <returns>A sorted list containing deltas</returns>
        public Dictionary<string, double> GetCellDeltas(string deltaType, uint cellLatIndex, uint cellLonIndex)
        {
            return InternalGrid[cellLatIndex, cellLonIndex].Deltas[deltaType];
        }

        /// <summary>
        /// A method to return all delta values in a particular grid cell
        /// </summary>
        /// <param name="cellLatIndex">Latitude index of grid cell</param>
        /// <param name="cellLonIndex">Longitude index of grid cell</param>
        /// <returns>A sorted list of sorted lists containing deltas</returns>
        public Dictionary<string, Dictionary<string, double>> GetCellDeltas(uint cellLatIndex, uint cellLonIndex)
        {
            return InternalGrid[cellLatIndex, cellLonIndex].Deltas;
        }

        
        /// <summary>
        /// Get a grid of values for an environmental data layer
        /// </summary>
        /// <param name="enviroVariable"> The name of the environmental data layer</param>
        /// <param name="timeInterval">The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <returns>The values in each grid cell</returns>
        public double[,] GetEnviroGrid(string enviroVariable,uint timeInterval)
        {
            // Check to see if environmental variable exists
            for (int ii = 0; ii < _NumLatCells; ii++)
            {
                for (int jj = 0; jj < _NumLonCells; jj++)
                {
                    if(InternalGrid[ii,jj] != null)
                        Debug.Assert(InternalGrid[ii, jj].CellEnvironment.ContainsKey(enviroVariable), "Environmental variable not found when running GetEnviroGrid");
                }
            }

            double[,] outputData = new double[_NumLatCells, _NumLonCells];

            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
            {
                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
                {
                    outputData[ii, jj] = InternalGrid[ii, jj].CellEnvironment[enviroVariable][timeInterval];
                }
            }

            return outputData;
        }

        /// <summary>
        /// Get a grid of values for an environmental data layer in specific cells
        /// </summary>
        /// <param name="enviroVariable">The name of the environmental data layer to return</param>
        /// <param name="timeInterval">The desired time interval for which to get data (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <param name="cellIndices">List of active cells in the model grid</param>
        /// <returns>The values in each grid cell</returns>
        public double[,] GetEnviroGrid(string enviroVariable, uint timeInterval, List<uint[]> cellIndices)
        {
            // Check to see if environmental variable exists
            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                if (InternalGrid[cellIndices[ii][0], cellIndices[ii][1]] != null)
                        Debug.Assert(InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].CellEnvironment.ContainsKey(enviroVariable), 
                            "Environmental variable not found when running GetEnviroGrid");
                
            }

            // Create grid to hold the data to return
            double[,] outputData = new double[_NumLatCells, _NumLonCells];

            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                outputData[cellIndices[ii][0], cellIndices[ii][1]] = InternalGrid[cellIndices[ii][0], cellIndices[ii][1]].CellEnvironment
                    [enviroVariable][timeInterval];
            }

            return outputData;
        }

        /// <summary>
        /// Return the total over the whole grid for an environmental variable
        /// </summary>
        /// <param name="enviroVariable">The environmental variable</param>
        /// <param name="timeInterval">The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <returns>The total of the variable over the whole grid</returns>
        public double GetEnviroGridTotal(string enviroVariable, uint timeInterval)
        {
            double[,] enviroGrid = GetEnviroGrid(enviroVariable,timeInterval);
            double enviroTotal = 0.0;

            for (int ii = 0; ii < _NumLatCells; ii+=GridCellRarefaction)
            {
                for (int jj = 0; jj < _NumLonCells; jj+=GridCellRarefaction)
                {
                    enviroTotal += enviroGrid[ii, jj];
                }
            }

            return enviroTotal;
        }

        /// <summary>
        /// Return the sum of an environmental variable over specific cells
        /// </summary>
        /// <param name="enviroVariable">The environmental variable</param>
        /// <param name="timeInterval">The desired time interval within the environmental variable (i.e. 0 if it is a yearly variable
        /// or the month index - 0=Jan, 1=Feb etc. - for monthly variables)</param>
        /// <param name="cellIndices">List of active cells in the model grid</param>
        /// <returns>The total of the variable over the whole grid</returns>
        public double GetEnviroGridTotal(string enviroVariable, uint timeInterval, List<uint[]> cellIndices)
        {
            double[,] enviroGrid = GetEnviroGrid(enviroVariable,timeInterval, cellIndices);
            double enviroTotal = 0.0;

            for (int ii = 0; ii < cellIndices.Count; ii++)
            {
                enviroTotal += enviroGrid[cellIndices[ii][0], cellIndices[ii][1]];
            }

            return enviroTotal;
        }

        /// <summary>
        /// Check to see if the top perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckTopPerimeterTraversable(uint latCell, uint lonCell, double gridCellRealm)
        {
            // Check to see if top perimeter is traversable
            if (InternalGrid[latCell + 1, lonCell].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell + 1), (lonCell) });

                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(1);
            }
        }

        /// <summary>
        /// Check to see if the top right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckTopRightPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {

            // Check to see if right perimeter is traversable
            if (InternalGrid[latCell + 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell + 1), (lonCellToGoTo) });
                
                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(2);
            }
        }

        /// <summary>
        /// Check to see if the right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckRightPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {
            // Check to see if right perimeter is traversable
            if (InternalGrid[latCell, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell), (lonCellToGoTo) });

                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(3);

            }
        }


        /// <summary>
        /// Check to see if the bottom right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckBottomRightPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {
            // Check to see if bottom right perimeter is traversable
            if (InternalGrid[latCell - 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell - 1), (lonCellToGoTo) });
                
                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(4);
            }
        }


        /// <summary>
        /// Check to see if the right perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckBottomPerimeterTraversable(uint latCell, uint lonCell, double gridCellRealm)
        {
            // Check to see if top perimeter is traversable
            if (InternalGrid[latCell - 1, lonCell].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell - 1), (lonCell) });

                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(5);
            
            }
        }


        /// <summary>
        /// Check to see if the bottom left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckBottomLeftPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {
            // Check to see if bottom right perimeter is traversable
            if (InternalGrid[latCell - 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell - 1), (lonCellToGoTo) });
                
                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(6);

            }
        }


        /// <summary>
        /// Check to see if the left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckLeftPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {
            // Check to see if left perimeter is traversable
            if (InternalGrid[latCell, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell), (lonCellToGoTo) });

                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(7);

            }
        }


        /// <summary>
        /// Check to see if the top left perimeter of the cell is traversable for dispersal (i.e. is from the same realm)
        /// </summary>
        /// <param name="latCell">The latitudinal cell index</param>
        /// <param name="lonCell">The longitudinal cell index</param>
        /// <param name="lonCellToGoTo">The index of the cell to go to (needs to take into account grid wrapping)</param>
        /// <param name="gridCellRealm">The grid cell realm</param>
        private void CheckTopLeftPerimeterTraversable(uint latCell, uint lonCell, uint lonCellToGoTo, double gridCellRealm)
        {
            // Check to see if bottom right perimeter is traversable
            if (InternalGrid[latCell + 1, lonCellToGoTo].CellEnvironment["Realm"][0] == gridCellRealm)
            {
                // Add the cell above to the list of cells that are dispersable to
                CellsForDispersal[latCell, lonCell].Add(new uint[2] { (latCell + 1), (lonCellToGoTo) });

                // Also add it to the directional list
                CellsForDispersalDirection[latCell, lonCell].Add(8);

            }
        }
 
        // Currently assumes that the grid does not run from -90 to 90 (in which case there would be transfer at top and bottom latitude)
        // Also needs checking to see if it works with a sub-grid
        /// <summary>
        /// Calculate the dispersable perimeter lengths of each of the grid cells
        /// </summary>
        private void CalculatePerimeterLengthsAndCellsDispersableTo()
        {
            int counter = 1;
            // Loop through grid cells
            for (uint ii = 0; ii < _NumLatCells; ii++)
			{
                // Bottom of the grid
			    if (ii == 0)
                {
                    // Loop through the longitude indices of each cell
                    for (uint jj = 0; jj < _NumLonCells; jj++)
			        {
                        // Get the realm of the cell (i.e. whether it is land or sea)
                        double GridCellRealm = InternalGrid[ii,jj].CellEnvironment["Realm"][0];
                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
                        {
                            Console.Write("\r{0} cells classified as neither land nor sea",counter);
                            counter++;
                            break;
                        }

                        // Check to see if we are at the left-most edge
			            if (jj == 0)
                        {
                            // Are we on a grid that spans the globe?
                            if ((_MaxLongitude - _MinLongitude) > 359.9)
                            {
                                // Check to see if the top perimeter is dispersable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                
                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);

                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);                  
                            }

                            // Otherwise, we are simply on a non-wrappable boundary. 
                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
                            else
                            {
                                // Check to see if the top perimeter is traversable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
                                
                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                
                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                            }
                        }
                        // Check to see if we are at the right-most edge
                        else if (jj == (_NumLonCells - 1))
                        {
                            // Are we on a grid that spans the globe?
                            if ((_MaxLongitude - _MinLongitude) > 359.9)
                            {
                                // Check to see if the top perimeter is traversable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
                                
                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }
                            // Otherwise, we are simply on a non-wrappable boundary. 
                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
                            else
                            {
                                // Check to see if the top perimeter is traversable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }
                        }

                        // Otherwise internal in the grid longitudinally
                        else
                        {
                            // Check to see if the top perimeter is traversable
                            CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                            // Check to see if the top right perimeter is dispersable
                            CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                
                            // Check to see if the right perimeter is dispersable
                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the left perimeter is dispersable
                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                            // Check to see if the top left perimeter is dispersable
                            CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                        }
			        }
                }

                // Top of the grid
                else if (ii == (_NumLatCells -1))
                {
                    // Loop through the longitude indices of each cell
                    for (uint jj = 0; jj < _NumLonCells; jj++)
			        {
                        // Get the realm of the cell (i.e. whether it is land or sea)
                        double GridCellRealm = InternalGrid[ii,jj].CellEnvironment["Realm"][0];
                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
                        {
                            Console.Write("\r{0} cells classified as neither land nor sea", counter);
                            counter++;
                            break;
                        }
                    
                        // Check to see if we are at the left-most edge
                        if (jj == 0)
                        {
                              // Are we on a grid that spans the globe?
                                if ((_MaxLongitude - _MinLongitude) > 359.9)
                                {
                                    // Check to see if the right perimeter is dispersable
                                    CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                    // Check to see if the bottom right perimeter is dispersable
                                    CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                
                                    // Check to see if the bottom perimeter is dispersable
                                    CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);

                                    // Check to see if the bottom left perimeter is dispersable
                                    CheckBottomLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);

                                    // Check to see if the left perimeter is dispersable
                                    CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
                                }
                                // Otherwise, we are simply on a non-wrappable boundary. 
                                // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
                                else
                                {
                                    // Check to see if the right perimeter is dispersable
                                    CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                    
                                    // Check to see if the bottom right perimeter is dispersable
                                    CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                    // Check to see if the bottom perimeter is dispersable
                                    CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
                                }
                        }
                        // Check to see if we are at the right-most edge
                        else if (jj == (_NumLonCells - 1))
                        {
                            // Are we on a grid that spans the globe?
                            if ((_MaxLongitude - _MinLongitude) > 359.9)
                            {
                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
                                
                                // Check to see if the bottom right perimeter is dispersable
                                CheckBottomRightPerimeterTraversable(ii, jj, 0, GridCellRealm);

                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
                                
                                // Check to see if the bottom left perimeter is dispersable
                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }

                            // Otherwise, we are simply on a non-wrappable boundary. 
                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
                            else
                            {
                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the bottom left perimeter is dispersable
                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }

                        }
                        // Otherwise, internal in the grid longitudinally
                        else
                        {
                            // Check to see if the right perimeter is dispersable
                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the bottom right perimeter is dispersable
                            CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the bottom perimeter is dispersable
                            CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);

                            // Check to see if the bottom left perimeter is dispersable
                            CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                            // Check to see if the left perimeter is dispersable
                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                        }
                    }

                }
                // Otherwise internal latitudinally
                else
                {
                    // Loop through the longitude indices of each cell
                    for (uint jj = 0; jj < _NumLonCells; jj++)
                    {
                        // Get the realm of the cell (i.e. whether it is land or sea)
                        double GridCellRealm = InternalGrid[ii, jj].CellEnvironment["Realm"][0];
                        if ((GridCellRealm != 1.0) && (GridCellRealm != 2.0))
                        {
                            Console.Write("\r{0} cells classified as neither land nor sea", counter);
                            counter++;
                            break;
                        }

                        // Check to see if we are at the left-most edge
                        if (jj == 0)
                        {
                            // Are we on a grid that spans the globe?
                            if ((_MaxLongitude - _MinLongitude) > 359.9)
                            {
                                // Check to see if the top perimeter is dispersable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the bottom right perimeter is dispersable
                                CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
                                
                                // Check to see if the bottom left perimeter is dispersable
                                CheckBottomLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);

                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, _NumLonCells - 1, GridCellRealm);
                            }
                            // Otherwise, we are simply on a non-wrappable boundary. 
                            // Assumes that we have a closed system on this boundary and that organisms cannot disperse through it
                            else
                            {
                                // Check to see if the top perimeter is dispersable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);
                                
                                // Check to see if the bottom right perimeter is dispersable
                                CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
                            }
                        }
                        // Check to see if we are at the rightmost edge
                        else if (jj == (_NumLonCells - 1))
                        {
                            // Are we on a grid that spans the globe?
                            if ((_MaxLongitude - _MinLongitude) > 359.9)
                            {
                                // Check to see if the top perimeter is dispersable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);
                                
                                // Check to see if the top right perimeter is dispersable
                                CheckTopRightPerimeterTraversable(ii, jj, 0, GridCellRealm);

                                // Check to see if the right perimeter is dispersable
                                CheckRightPerimeterTraversable(ii, jj, 0, GridCellRealm);
                                
                                // Check to see if the bottom right perimeter is dispersable
                                CheckBottomRightPerimeterTraversable(ii, jj, 0, GridCellRealm);

                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);
                                
                                // Check to see if the bottom left perimeter is dispersable
                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                                
                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }
                            else
                            {
                                // Check to see if the top perimeter is dispersable
                                CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the bottom perimeter is dispersable
                                CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);

                                // Check to see if the bottom left perimeter is dispersable
                                CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the left perimeter is dispersable
                                CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                                // Check to see if the top left perimeter is dispersable
                                CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                            }
                        }
                        // Otherwise internal in the grid both latitudinally and longitudinally - the easiest case
                        else
                        {
                            // Check to see if the top perimeter is dispersable
                            CheckTopPerimeterTraversable(ii, jj, GridCellRealm);

                            // Check to see if the top right perimeter is dispersable
                            CheckTopRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the right perimeter is dispersable
                            CheckRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the bottom right perimeter is dispersable
                            CheckBottomRightPerimeterTraversable(ii, jj, jj + 1, GridCellRealm);

                            // Check to see if the bottom perimeter is dispersable
                            CheckBottomPerimeterTraversable(ii, jj, GridCellRealm);

                            // Check to see if the bottom left perimeter is dispersable
                            CheckBottomLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                            // Check to see if the left perimeter is dispersable
                            CheckLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);

                            // Check to see if the top left perimeter is dispersable
                            CheckTopLeftPerimeterTraversable(ii, jj, jj - 1, GridCellRealm);
                        }
                    }
                 }
			}
            Console.WriteLine("\n");

        }

        /// <summary>
        /// Given a grid cell from where a cohort is dispersing, select at random a grid cell for it to disperse to from those that exist within the 
        /// same realm
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the cell from which the cohort is dispersing</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the cell from which the cohort is dispersing</param>
        /// <returns></returns>
        public uint[] GetRandomGridCellToDisperseTo(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            // Select a cell at random
            int CellPickedAtRandom = (int)Math.Floor(RandomNumberGenerator.GetUniform() * CellsForDispersal[fromCellLatIndex, fromCellLonIndex].Count);

            // Return the coordinates of that cell
            return CellsForDispersal[fromCellLatIndex, fromCellLonIndex][CellPickedAtRandom];
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the north of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the north of the focal grid cell</returns>
        public uint[] CheckDispersalNorth(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] NorthCell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 1)
                {
                    NorthCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                    break;
                }
            }

            return NorthCell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the east of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indices of the cell that lies to the east of the focal grid cell</returns>
        public uint[] CheckDispersalEast(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] EastCell = new uint[2] {9999999,9999999};

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 3)
                {
                    EastCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                    break;
                }
            }

            return EastCell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the south of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the south of the focal grid cell</returns>
        public uint[] CheckDispersalSouth(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] SouthCell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 5)
                {
                    SouthCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return SouthCell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the west of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the west of the focal grid cell</returns>
        public uint[] CheckDispersalWest(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] WestCell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 7)
                {
                    WestCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return WestCell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the northeast of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the northeast of the focal grid cell</returns>
        public uint[] CheckDispersalNorthEast(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] NECell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 2)
                {
                    NECell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return NECell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the southeast of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the southeast of the focal grid cell</returns>
        public uint[] CheckDispersalSouthEast(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] SECell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 4)
                {
                    SECell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return SECell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the southwest of the focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the southwest of the focal grid cell</returns>
        public uint[] CheckDispersalSouthWest(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] SWCell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 6)
                {
                    SWCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return SWCell;
        }

        /// <summary>
        /// Get the longitudinal and latitudinal indices of the cell that lies to the northwest of focal grid cell, if a viable cell to disperse to
        /// </summary>
        /// <param name="fromCellLatIndex">The latitudinal index of the focal grid cell</param>
        /// <param name="fromCellLonIndex">The longitudinal index of the focal grid cell</param>
        /// <returns>The longitudinal and latitudinal cell indcies of the cell that lies to the northwest of the focal grid cell</returns>
        public uint[] CheckDispersalNorthWest(uint fromCellLatIndex, uint fromCellLonIndex)
        {
            uint[] NWCell = new uint[2] { 9999999, 9999999 };

            for (int ii = 0; ii < CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex].Count; ii++)
            {
                if (CellsForDispersalDirection[fromCellLatIndex, fromCellLonIndex][ii] == 8)
                {
                    NWCell = CellsForDispersal[fromCellLatIndex, fromCellLonIndex][ii];
                }
            }

            return NWCell;
        }

    }
}
