using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;

using Newtonsoft.Json;

namespace Parser
{
    public class Program
    {
        public const string DIR = @"P:\Honours\";

        public static void Main()
        { 
            
            //File.WriteAllText(Path.Combine(DIR, "crisprs.json"), crispr.Json());



            //Sequence cas9 = new Sequence("cas9.fasta");
            //Sequence streptococcus = new Sequence("streptococcus.fasta");
            //Sequence aureus = new Sequence("aureus.fasta");s


            //List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = OrderedDyadFrequencies(aureus, REPEAT_MIN, REPEAT_MAX);
            //List<Sequence> consensuses = new List<Sequence>();
            //foreach (KeyValuePair<Sequence, int> consensus in ordered_dyad_frequencies)
            //{
            //    Console.WriteLine("Adding consensus {0} with global frequency {1}", consensus.Key, consensus.Value);
            //    consensuses.Add(consensus.Key);
            //}

            //Crisprs crisprs = DiscoverCrisprs(aureus, consensuses);
            //Console.WriteLine(crisprs);


            Console.WriteLine("Press any key to quit.");
            Console.ReadKey();
        }

        private static readonly HttpClient client = new HttpClient();

        public static void BLAST()
        {

        }

    }
}