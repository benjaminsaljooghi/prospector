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
        public class Crispr
        {
            public Crispr(Sequence consensus)
            {
                Consensus = consensus;
                Repeats = new List<int>();
            }

            public Sequence Consensus { get; }

            public List<int> Repeats { get; }

            //public void AddRepeats(List<int> repeat_indices)
            //{
            //    Repeats.AddRange(repeat_indices);
            //}

            public void AddRepeat(int repeat_index)
            {
                Repeats.Add(repeat_index);
            }
        }

        public class Crisprs
        {
            List<Crispr> crisprs = new List<Crispr>();
            
            public void RegisterCrispr(Crispr new_crispr)
            {
                foreach (Crispr crispr in crisprs)
                {
                    if (new_crispr.Consensus.Equals(crispr.Consensus))
                    {
                        throw new Exception("Consensus already registered.");
                    }
                }

                crisprs.Add(new_crispr);

                // Mutant check
                //foreach (Crispr crispr in crisprs)
                //{
                //    if (Mutant(new_crispr.Consensus, crispr.Consensus))
                //    {
                //        crispr.AddRepeats(new_crispr.Repeats);
                //    }
                //}
            }
        }



        public static void Main()
        {
            Sequence cas9 = new Sequence("cas9.fasta");
            Sequence streptococcus = new Sequence("streptococcus.fasta");
            Sequence aureus = new Sequence("aureus.fasta");

            Crisprs crisprs = new Crisprs();

            List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = OrderedDyadFrequencies(streptococcus, 28, 40);
            List<Sequence> consensuses = new List<Sequence>();
            foreach (KeyValuePair<Sequence, int> consensus in ordered_dyad_frequencies)
            {
                consensuses.Add(consensus.Key);
            }

            DiscoverCrisprs(streptococcus, consensuses);


            Console.WriteLine("Press enter to quit.");
            Console.ReadLine();
        }

        public static Dictionary<Sequence, int> DyadFrequencies(Sequence genome, int k)
        {
            Dictionary<Sequence, int> frequencies = new Dictionary<Sequence, int>();
            int num_kmers = genome.Length - k + 1;
            for (int i = 0; i < num_kmers; i++)
            {
                Sequence kmer = genome.Substring(i, k);
                if (!Dyad(kmer))
                {
                    continue;
                }
                if (!frequencies.ContainsKey(kmer))
                {
                    frequencies.Add(kmer, 0);
                }
                frequencies[kmer] += 1;
            }
            return frequencies;
        }

        public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Sequence genome, int k_min, int k_max)
        {
            List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = new List<KeyValuePair<Sequence, int>>();
            for (int k = k_min; k < k_max; k++)
            {
                ordered_dyad_frequencies.AddRange(OrderedDyadFrequencies(genome, k));
            }
            return ordered_dyad_frequencies;
        }

        public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Sequence genome, int k)
        {
            return OrderedDyadFrequencies(DyadFrequencies(genome, k));
        }

        public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Dictionary<Sequence, int> dyad_frequencies)
        {
            List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = dyad_frequencies.ToList();
            ordered_dyad_frequencies.Sort((a, b) => b.Value - a.Value);
            return ordered_dyad_frequencies;
        }

        public static Crisprs DiscoverCrisprs(Sequence genome, List<Sequence> consensuses)
        {
            Crisprs crisprs = new Crisprs();
            foreach (Sequence consensus in consensuses)
            {
                Crispr crispr = DiscoverCrispr(genome, consensus);
                crisprs.RegisterCrispr(crispr);
            }
            return crisprs;
        }

        public static Crispr DiscoverCrispr(Sequence genome, Sequence consensus)
        {
            Crispr crispr = new Crispr(consensus);

            int k = consensus.Length;
            int spacer_skip = 10;

            // Upstream scan
            int index = consensus.Pos + k + spacer_skip;
            const int reset = 100;
            int countdown = reset;
            while (countdown-- > 0)
            {
                Sequence kmer = genome.Substring(index++, k);
                if (Mutant(consensus, kmer))
                {
                    crispr.AddRepeat(kmer.Pos);
                    index = kmer.Pos + k + spacer_skip;
                    countdown = reset;
                }
            }

            // Downstream scan
            index = consensus.Pos - k - spacer_skip;
            countdown = reset;
            while (countdown-- > 0)
            {
                Sequence kmer = genome.Substring(index--, k);
                if (Mutant(consensus, kmer))
                {
                    crispr.AddRepeat(kmer.Pos);
                    index = kmer.Pos - k - spacer_skip;
                    countdown = reset;
                }
            }

            return crispr;
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