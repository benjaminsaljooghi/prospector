using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Net.Http;

using Newtonsoft.Json;

namespace Parser
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

        public void AddRepeat(int repeat_index)
        {
            Repeats.Add(repeat_index);
        }

        public override string ToString()
        {
            string repeats = "";
            foreach (int repeat in Repeats)
            {
                repeats += repeat + " ";
            }
            return string.Format("{0} : {1}", Consensus, repeats);
        }

        public string Json()
        {
            return JsonConvert.SerializeObject(this);
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
        }

        public override string ToString()
        {
            string result = "";
            foreach (Crispr crispr in crisprs)
            {
                result += string.Format("{0}\n", crispr);
            }
            return result;
        }
    }



    public class Algorithm
    {
        public const string dir = @"P:\Honours\";

        public const int REPEATS_MIN = 3;

        public const int REPEAT_MIN = 20;
        public const int REPEAT_MAX = 60;

        public const int SPACER_MIN = 21;
        public const int SPACER_MAX = 72;

        public static void Main()
        {
            Sequence seq = new Sequence("ACGT", 0);
            //string a = JsonConvert.SerializeObject(seq);
            //Console.WriteLine(a);

            Crispr crispr = new Crispr(seq);
            crispr.AddRepeat(3);
            crispr.AddRepeat(73);

            File.WriteAllText(Path.Combine(Sequence.dir, "crisprs.json"), crispr.Json());

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



        public static Crisprs DiscoverCrisprs(Sequence genome, List<Sequence> consensuses)
        {
            Crisprs crisprs = new Crisprs();
            foreach (Sequence consensus in consensuses)
            {
                Crispr crispr = DiscoverCrispr(genome, consensus);
                if (crispr != null)
                {
                    crisprs.RegisterCrispr(crispr);
                }
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
                try
                {
                    Sequence kmer = genome.Substring(index++, k);
                    if (Sequence.Mutant(consensus, kmer))
                    {
                        crispr.AddRepeat(kmer.Pos);
                        index = kmer.Pos + k + spacer_skip;
                        countdown = reset;
                    }
                }
                catch (ArgumentOutOfRangeException)
                {
                    Console.WriteLine("Index was out of bounds. Continuing...");
                    break;
                }

            }

            // Downstream scan
            index = consensus.Pos - k - spacer_skip;
            countdown = reset;
            while (countdown-- > 0)
            {
                try
                {
                    Sequence kmer = genome.Substring(index--, k);
                    if (Sequence.Mutant(consensus, kmer))
                    {
                        crispr.AddRepeat(kmer.Pos);
                        index = kmer.Pos - k - spacer_skip;
                        countdown = reset;
                    }
                }
                catch (ArgumentOutOfRangeException)
                {
                    Console.WriteLine("Index was out of bounds. Continuing...");
                    break;
                }
            }

            return crispr.Repeats.Count >= REPEATS_MIN ? crispr : null;
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