using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.Net.Http;

using Microsoft.FSharp.Core;
using Microsoft.FSharp.Collections;

namespace Parser
{
    using System.Collections;
    using static Sequence;

   
    public class Algorithms
    {

        public static void Main()
        {
            Sequence cas9 = new Sequence("cas9.fasta");
            Sequence streptococcus = new Sequence("streptococcus.fasta");
            Sequence aureus = new Sequence("aureus.fasta");
            PrintClusters(aureus);
            //int array_begin = 860748;
            //Sequence upstream = streptococcus.Substring(array_begin - 10000, 10000);
            //Console.WriteLine("upstream length: " + upstream.Length);
            //Console.WriteLine("cas9 length: " + cas9.Length);
        }











        public static void LocalAlignment(Sequence big, Sequence small)
        {

        }

        public static void PrintClusters(Sequence sequence)
        {
            foreach (List<Sequence> kmers in KmerWindow(sequence, 30, 40))
            {
                List<Sequence> dyads = Dyads(kmers);
                Dictionary<Sequence, List<int>> arrays = new Dictionary<Sequence, List<int>>();
                foreach (Sequence dyad in dyads)
                {
                    if (!arrays.ContainsKey(dyad))
                    {
                        arrays.Add(dyad, new List<int>());
                    }
                    arrays[dyad].Add(dyad.Pos);
                }

                var myList = arrays.ToList();

                myList.Sort(
                    delegate (KeyValuePair<Sequence, List<int>> pair1,
                    KeyValuePair<Sequence, List<int>> pair2)
                    {
                        return pair1.Value.Count().CompareTo(pair2.Value.Count());
                    }
                );

                foreach (KeyValuePair<Sequence, List<int>> dyad_positions in myList)
                {
                    if (dyad_positions.Value.Count() > 1)
                    {
                        string positions = String.Join(", ", dyad_positions.Value);
                        Console.WriteLine("{0}: {1}", dyad_positions.Key, positions);
                    }
                }
            }
        }


        private static readonly HttpClient client = new HttpClient();

        public static void BLAST()
        {
            
        }

    }
}