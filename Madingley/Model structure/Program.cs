using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Diagnostics;

namespace Madingley
{

    /// <summary>
    /// The entry point for the model
    /// <todoM>Write model output to an output file</todoM>
    /// </summary>
    class Program
    {
        /// <summary>
        /// Starts a model run or set of model runs
        /// </summary>
        static void Main()
        {

            // Write out model details to the console
            Console.ForegroundColor = ConsoleColor.Green;
            Console.WriteLine("Madingley model v. 0.3333333\n");
            Console.ForegroundColor = ConsoleColor.White;

            // Declare an instance of RunSimulations
            RunSimulations MakeSimulations = new RunSimulations();

            // Specify the working directory
            //string OutputDir = "C:/Users/derekt/Dropbox/Madingley stuff/Model outputs/MadingleyOutputs" +
            //string OutputDir = "C:/Users/derekt/desktop/MadingleyOutputs" +
            //string OutputDir = "Ensemble" +
            string OutputDir = "Bergmann_ko" +
            //string OutputDir = "C:/Users/mikeha/Work/Research/Visual Studio 2010/Madingley/madingley outputs" +
                System.DateTime.Now.Year + "-"
                + System.DateTime.Now.Month + "-"
                + System.DateTime.Now.Day + "_"
                + System.DateTime.Now.Hour + "."
                + System.DateTime.Now.Minute + "."
                + System.DateTime.Now.Second + "/";

            // Create the working directory if this does not already exist
            System.IO.Directory.CreateDirectory(OutputDir);

            // Declare an instance of ScenarioParameterInitialisation to read in the parameters for this model run or set of runs
            ScenarioParameterInitialisation Scenarios = new ScenarioParameterInitialisation("Scenarios.csv", OutputDir);


            // Run the desired simulation or batch of simulations
            MakeSimulations.RunAllSimulations("EcosystemModelInitialisation.csv", Scenarios, OutputDir);

            Console.ReadKey();
        }

        
    }
}