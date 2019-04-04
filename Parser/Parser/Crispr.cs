using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Parser
{
    public class Crispr
    {
        public const int REPEATS_MIN = 3;

        public const int REPEAT_MIN = 20;
        public const int REPEAT_MAX = 60;

        public const int SPACER_MIN = 21;
        public const int SPACER_MAX = 72;

        public const int SCAN_DOMAIN = 100;

        public Crispr(Sequence consensus)
        {
            Consensus = consensus;
            Repeats = new List<int>();
        }

        public Sequence Consensus { get; private set; }

        public List<int> Repeats { get; }

        public int End
        {
            get
            {
                Repeats.Sort();
                int last_repeat = Repeats[Repeats.Count - 1];
                return last_repeat + Consensus.Length;
            }
        }

        public int Start
        {
            get
            {
                Repeats.Sort();
                int first_repeat = Repeats[0];
                return first_repeat;
            }
        }

        public void AddRepeat(int repeat_index)
        {
            Repeats.Add(repeat_index);
        }

        public void AddRepeats(List<int> repeats)
        {
            Repeats.AddRange(repeats);
        }

        public void UpdateConsensus(Sequence genome)
        {
            int consensus_index = Repeats.GroupBy(i => i).OrderByDescending(grp => grp.Count()).Select(grp => grp.Key).First();
            Consensus = genome.Substring(consensus_index, Consensus.Length);
        }

        public override string ToString()
        {
            string repeats = "";
            foreach (int repeat in Repeats)
            {
                repeats += repeat + " ";
            }
            return string.Format($"{Consensus,-REPEAT_MAX} : {repeats}");
        }

        public override bool Equals(object obj)
        {
            Crispr c = obj as Crispr;
            if (c == null)
            {
                return false;
            }
            bool consensus = Consensus == c.Consensus;
            bool repeats = Repeats.All(c.Repeats.Contains) && Repeats.Count == c.Repeats.Count;
            return consensus && repeats;
        }

        public override int GetHashCode()
        {
            int hc = Consensus.GetHashCode() * 7;
            foreach (int repeat in Repeats)
            {
                hc ^= repeat.GetHashCode();
            }
            return hc;
        }

        public static int MinDistance(Crispr a, Crispr b)
        {
            if (a.Start < b.Start)
            {
                return b.Start - a.End;
            }
            else if (a.Start > b.Start)
            {

            }
            else
            {
                throw new Exception("Two CRISPRs have the same start position. This should not be possible.");
            }
        }

        public static Crispr DiscoverCrispr(Sequence genome, Sequence consensus)
        {
            Crispr crispr = new Crispr(consensus);

            int k = consensus.Length;
            int spacer_skip = 10;

            // Upstream scan
            int index = consensus.Pos + k + spacer_skip;
            const int reset = SCAN_DOMAIN;
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
        public HashSet<Crispr> Clusters { get; } = new HashSet<Crispr>();

        public void RegisterCrispr(Crispr new_crispr)
        {
            //foreach (Crispr crispr in Clusters)
            //{
            //    if (new_crispr.Consensus.Equals(crispr.Consensus))
            //    {
            //        throw new Exception("Consensus already registered.");
            //        //Console.WriteLine("Tried to register a crispr with a consensus that is already registered. Merging the two CRISPRs...");
            //        //crispr.AddRepeats(new_crispr.Repeats);
            //    }
            //}
            Clusters.Add(new_crispr);
        }

        public void RegisterCrisprs(Crisprs crisprs)
        {
            Clusters.UnionWith(crisprs.Clusters);
        }

        public override string ToString()
        {
            string result = "";
            foreach (Crispr crispr in Clusters)
            {
                result += string.Format($"{crispr}\n");
            }
            return result;
        }

        public void PrintMutantConsensuses()
        {
            foreach (Crispr crispr in Clusters)
            {
                foreach (Crispr sub_crispr in Clusters)
                {
                    if (crispr.Equals(sub_crispr))
                    {
                        continue;
                    }
                    if (Sequence.Mutant(crispr.Consensus, sub_crispr.Consensus, allow_discrepant_lengths:true))
                    {
                        Console.WriteLine("Found a mutant.");
                    }
                }
            }
        }

        public static Crisprs DiscoverCrisprs(Sequence genome, int k)
        {
            Crisprs crisprs = new Crisprs();
            int num_kmers = genome.Length - k + 1;
            for (int i = 0; i < num_kmers; i++)
            {
                Console.WriteLine($"{i} at entry");

                Sequence kmer = genome.Substring(i, k);

                if (!Sequence.Dyad(kmer))
                {
                    continue;
                }

                Console.WriteLine($"found CRISPR at {i} of genome len {genome.Length}");

                Crispr crispr = Crispr.DiscoverCrispr(genome, kmer);
                if (crispr != null)
                {
                    crisprs.RegisterCrispr(crispr);
                    i = crispr.End + 1;
                }
            }
            return crisprs;
        }

        public static Crisprs DiscoverCrisprs(Sequence genome, int k_begin, int k_end)
        {
            Crisprs crisprs = new Crisprs();
            for (int k = k_begin; k <= k_end; k++)
            {
                crisprs.RegisterCrisprs(DiscoverCrisprs(genome, k));   
            }
            return crisprs;
        }
    }
}
