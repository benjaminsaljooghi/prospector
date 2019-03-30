using Newtonsoft.Json;
using System;
using System.Collections.Generic;

namespace Parser
{
    public class Crispr
    {
        public const int REPEATS_MIN = 3;

        public const int REPEAT_MIN = 20;
        public const int REPEAT_MAX = 60;

        public const int SPACER_MIN = 21;
        public const int SPACER_MAX = 72;

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



    }

    public class Crisprs
    {
        public List<Crispr> Clusters { get; } = new List<Crispr>();

        public void RegisterCrispr(Crispr new_crispr)
        {
            foreach (Crispr crispr in Clusters)
            {
                if (new_crispr.Consensus.Equals(crispr.Consensus))
                {
                    throw new Exception("Consensus already registered.");
                }
            }
            Clusters.Add(new_crispr);
        }

        
        public override string ToString()
        {
            string result = "";
            foreach (Crispr crispr in Clusters)
            {
                result += string.Format("{0}\n", crispr);
            }
            return result;
        }

        public static Crisprs DiscoverCrisprs(Sequence genome, List<Sequence> consensuses)
        {
            Crisprs crisprs = new Crisprs();
            foreach (Sequence consensus in consensuses)
            {
                Crispr crispr = Crispr.DiscoverCrispr(genome, consensus);
                if (crispr != null)
                {
                    crisprs.RegisterCrispr(crispr);
                }
            }
            return crisprs;
        }
    }
}
