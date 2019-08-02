using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;

using Newtonsoft.Json;

using static System.IO.File;
using static System.IO.Path;

namespace Prospector
{
    public class Program
    { 
        public static void Serialize(string path, object obj)
        {
            WriteAllText(path, JsonConvert.SerializeObject(obj));
        }

        public static T Deserialize<T>(string path)
        {
            return JsonConvert.DeserializeObject<T>(ReadAllText (path));
        }

        public static void Main()
        {
            Stopwatch stopwatch = new Stopwatch();
            stopwatch.Start();

            string fasta = @"P:\CRISPR\data\pyogenes.fasta";
            Sequence genome = new Sequence(fasta);

            Crisprs crisprs = Crisprs.DiscoverCrisprs(genome, 35, 39);

            Console.WriteLine("Before merging crisprs:");
            crisprs.SortByStartPos();
            Console.WriteLine(crisprs);

            Console.WriteLine("After merging crisprs:");
            crisprs.MergeCrisprs();
            crisprs.SortByStartPos();
            Console.WriteLine(crisprs);

            stopwatch.Stop();
            Console.WriteLine($"Completed in {stopwatch.ElapsedMilliseconds / 1000.0} seconds");

            Console.WriteLine("Press any key to quit.");
            Console.ReadKey();
        }
    }
}