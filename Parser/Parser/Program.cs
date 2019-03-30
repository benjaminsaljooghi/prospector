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


        //public static void PrintClusters(Sequence sequence)
        //{
        //    foreach (List<Sequence> kmers in KmerWindow(sequence, 23, 40))
        //    {
        //        List<Sequence> dyads = Dyads(kmers);
        //        Dictionary<Sequence, List<int>> arrays = new Dictionary<Sequence, List<int>>();
        //        foreach (Sequence dyad in dyads)
        //        {
        //            if (!arrays.ContainsKey(dyad))
        //            {
        //                arrays.Add(dyad, new List<int>());
        //            }
        //            arrays[dyad].Add(dyad.Pos);
        //        }

        //        var myList = arrays.ToList();

        //        myList.Sort(
        //            delegate (KeyValuePair<Sequence, List<int>> pair1,
        //            KeyValuePair<Sequence, List<int>> pair2)
        //            {
        //                return pair1.Value.Count().CompareTo(pair2.Value.Count());
        //            }
        //        );

        //        foreach (KeyValuePair<Sequence, List<int>> dyad_positions in myList)
        //        {
        //            if (dyad_positions.Value.Count() > 1)
        //            {
        //                string positions = String.Join(", ", dyad_positions.Value);
        //                Console.WriteLine("{0}: {1}", dyad_positions.Key, positions);
        //            }
        //        }
        //    }
        //}


        private static readonly HttpClient client = new HttpClient();

        public static void BLAST()
        {

        }

    }
}