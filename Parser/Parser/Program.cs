using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Newtonsoft.Json;

using static System.IO.File;


using static System.IO.Path;

namespace Parser
{
    public class Program
    {
        static Bio.IO.FastA.FastAParser fasta = new Bio.IO.FastA.FastAParser();

        public static void Serialize(string path, object obj)
        {
            WriteAllText(path, JsonConvert.SerializeObject(obj));
        }

        public static T Deserialize<T>(string path)
        {
            return JsonConvert.DeserializeObject<T>(ReadAllText (path));
        }

        public static void Execution(Sequence bacterium)
        {
            //Crisprs crisprs = Crisprs.DiscoverCrisprs(bacterium, Crispr.REPEAT_MIN, Crispr.REPEAT_MAX);
            //Serialize(Path.Combine(DIR, "aureus_crisprs.json"), crisprs);
            //Crisprs crisprs = Deserialize<Crisprs>(Path.Combine(DIR, "aureus_crisprs.json"));

            //Console.WriteLine("Before merging crisprs:");
            //crisprs.SortByStartPos();
            //Console.WriteLine(crisprs);

            //Console.WriteLine("After merging crisprs:");
            //crisprs.MergeCrisprs();
            //crisprs.SortByStartPos();
            //Console.WriteLine(crisprs); ;
        }


        public static IEnumerable<Bio.ISequence> ReadFasta(string path)
        {
            Stream reader = new FileStream(path, FileMode.Open);
            return fasta.Parse(reader);
        }


        //Sequence bacterium = new Sequence(bacterium_path);
        //Execution(bacterium);
        //var reader = new StreamReader(file);
        //while (reader.ReadLine().StartsWith(">")) ;
        //Seq = reader.ReadToEnd().Replace("\n", "");
        //Start = 0;

        public static void Main()
        {
            const string BACTERIA = @"P:\Honours\bacteria";
            const string PHAGES = @"P:\Honours\phage";
            string bacterium_fasta = Combine(BACTERIA, "pyogenes.fasta");
            string phages_fasta = Combine(PHAGES, "phages.fasta");
            string spacers_fasta = Combine(PHAGES, "pyogenes_spacers.fa");

            //Bio.ISequence bacterium = ReadFasta(bacterium_fasta).ElementAt(0);
            var spacers = ReadFasta(spacers_fasta);
            var phages = ReadFasta(phages_fasta);

            var aligner = new Bio.Algorithms.Alignment.SmithWatermanAligner();

        
            Quit();
        }

        public static void Quit()
        {
            Console.WriteLine("Press any key to quit.");
            Console.ReadKey();
        }


    }
}