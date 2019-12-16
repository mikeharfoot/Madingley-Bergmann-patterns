using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;


using Microsoft.Research.Science.Data;
using Microsoft.Research.Science.Data.Imperative;

namespace Madingley
{
    public class InputModelState
    {

        /// <summary>
        /// The handler for the cohorts in this grid cell
        /// </summary>
        private GridCellCohortHandler[,] _GridCellCohorts;
        /// <summary>
        /// Get or set the cohorts in this grid cell
        /// </summary>
        public GridCellCohortHandler[,] GridCellCohorts
        {
            get { return _GridCellCohorts; }
            set { _GridCellCohorts = value; }
        }

        /// <summary>
        /// The handler for the stocks in this grid cell
        /// </summary>
        private GridCellStockHandler[,] _GridCellStocks;
        /// <summary>
        /// Get or set the stocks in this grid cell
        /// </summary>
        public GridCellStockHandler[,] GridCellStocks
        {
            get { return _GridCellStocks; }
            set { _GridCellStocks = value; }
        }

        private Boolean _InputState = false;

        public Boolean InputState
        {
            get { return _InputState; }
            set { _InputState = value; }
        }
        

        public InputModelState(string outputPath, string filename)
        {

            //Set the input state flag to be true
            _InputState = true;
                                    
            // Construct the string required to access the file using Scientific Dataset
            string _ReadFileString = "msds:nc?file=" + outputPath + filename +".nc&openMode=readOnly";

            // Open the data file using Scientific Dataset
            DataSet StateDataSet = DataSet.Open(_ReadFileString);


            float[] Latitude = StateDataSet.GetData<float[]>("Latitude");
            float[] Longitude = StateDataSet.GetData<float[]>("Longitude");
            float[] CohortFunctionalGroup = StateDataSet.GetData<float[]>("Cohort Functional Group");
            float[] Cohort = StateDataSet.GetData<float[]>("Cohort");

            float[] StockFunctionalGroup = StateDataSet.GetData<float[]>("Stock Functional Group");
            float[] Stock = StateDataSet.GetData<float[]>("Stock");

            int NumDims;

            String TempVariableName;
            SortedList<string, double[, , ,]> ModelCohortStateList = new SortedList<string, double[, , ,]>();
            SortedList<string, double[, , ,]> ModelStockStateList = new SortedList<string, double[, , ,]>();
            



            for (int i = 0; i < StateDataSet.Variables.Count(); i++)
			{
                NumDims = StateDataSet.Variables.ElementAt(i).Dimensions.Count;

                if (NumDims == 4)
                {
                    TempVariableName = StateDataSet.Variables.ElementAt(i).Name;
                    if (TempVariableName.Contains("Cohort"))
                    {
                        ModelCohortStateList.Add(TempVariableName,StateDataSet.GetData<double[, , ,]>(TempVariableName));
                    }
                    else if(TempVariableName.Contains("Stock"))
                    {
                        ModelStockStateList.Add(TempVariableName, StateDataSet.GetData<double[, , ,]>(TempVariableName));
                    }

                }
			}

            _GridCellCohorts = new GridCellCohortHandler[Latitude.Length, Longitude.Length];
            _GridCellStocks = new GridCellStockHandler[Latitude.Length, Longitude.Length];


            for (int Lat = 0; Lat < Latitude.Length; Lat++ )
            {
                for (int Lon = 0; Lon < Longitude.Length; Lon++ )
                {
                    _GridCellCohorts[Lat, Lon] = new GridCellCohortHandler(CohortFunctionalGroup.Length);
                    _GridCellStocks[Lat, Lon] = new GridCellStockHandler(StockFunctionalGroup.Length);

                    /*
                     * For each lat and lon need to loop over all cohorts and stocks adding each cohort/stock into the appropriate handler
                     */

                    for (int fg = 0; fg < CohortFunctionalGroup.Length; fg++)
                    {
                        _GridCellCohorts[Lat, Lon][fg] = new List<Cohort>();

                        for (int c = 0; c < Cohort.Length; c++)
                        {

                            long temp = 0;

                            Cohort TempCohort = new Cohort(
                                (byte)fg,
                                ModelCohortStateList["CohortJuvenileMass"][Lat, Lon, fg, c],
                                ModelCohortStateList["CohortAdultMass"][Lat, Lon, fg, c],
                                ModelCohortStateList["CohortIndividualBodyMass"][Lat, Lon, fg, c],
                                ModelCohortStateList["CohortCohortAbundance"][Lat, Lon, fg, c],
                                Math.Exp(ModelCohortStateList["CohortLogOptimalPreyBodySizeRatio"][Lat, Lon, fg, c]),
                                Convert.ToUInt16(ModelCohortStateList["CohortBirthTimeStep"][Lat, Lon, fg, c]),
                                ModelCohortStateList["CohortProportionTimeActive"][Lat, Lon, fg, c], ref temp,
                                ModelCohortStateList["CohortTrophicIndex"][Lat, Lon, fg, c],
                                false);


                            _GridCellCohorts[Lat, Lon][fg].Add(TempCohort);

                        }

                    }

                    for (int fg = 0; fg < StockFunctionalGroup.Length; fg++)
                    {
                        _GridCellStocks[Lat, Lon][fg] = new List<Stock>();

                        for (int c = 0; c < Stock.Length; c++)
                        {

                            Stock TempStock = new Stock(
                                (byte)fg,
                                ModelStockStateList["StockIndividualBodyMass"][Lat, Lon, fg, c],
                                ModelStockStateList["StockTotalBiomass"][Lat, Lon, fg, c]);

                            _GridCellStocks[Lat, Lon][fg].Add(TempStock);

                        }
                    }

                }
            }

        }


    }
}
