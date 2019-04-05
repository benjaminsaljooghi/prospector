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

        public const int SCAN_DOMAIN = 1000;

        public Crispr()
        {
            Repeats = new List<Sequence>();
        }

        public List<Sequence> Repeats
        {
            get;
        }

        public Sequence Consensus
        {
            get
            {
                return Repeats.GroupBy(i => i).OrderByDescending(grp => grp.Count()).Select(grp => grp.Key).First();
            }
        }

        public Sequence First
        {
            get
            {
                Repeats.Sort(Sequence.CompareStart);
                return Repeats[0];
            }
        }

        public Sequence Last
        {
            get
            {
                Repeats.Sort(Sequence.CompareStart);
                return Repeats[Repeats.Count - 1];
            }                
        }

        public void AddRepeat(Sequence repeat)
        {
            Repeats.Add(repeat);
        }

        public void AddRepeats(List<Sequence> repeats)
        {
            Repeats.AddRange(repeats);
        }

        public override string ToString()
        {
            Repeats.Sort(Sequence.CompareStart);
            string repeats = "";
            foreach (Sequence repeat in Repeats)
            {
                repeats += repeat.Start + " ";
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
            foreach (Sequence repeat in Repeats)
            {
                hc += repeat.GetHashCode();
            }
            return hc;
        }

        public static int MinDistance(Crispr a, Crispr b)
        {
            int a_start = a.First.Start;
            int a_end = a.Last.End;

            int b_start = b.First.Start;
            int b_end = b.Last.End;

            if (a_start == b_start)
            {
                return 0;
            }
            else if (a_start < b_start)
            {
                return b_start - a_end;
            }
            else
            {
                return a_start- b_end;
            }
        }

        public static bool Mergeable(Crispr a, Crispr b)
        {
            bool mutant = Sequence.Mutant(a.Consensus, b.Consensus, true);
            bool proximal = MinDistance(a, b) <= SCAN_DOMAIN;
            return mutant && proximal;
        }

        public static int CompareByStart(Crispr a, Crispr b)
        {
            return a.First.Start - b.First.Start;
        }

        public static Crispr DiscoverCrispr(Sequence genome, Sequence consensus)
        {
            Crispr crispr = new Crispr();
            crispr.AddRepeat(consensus);

            int k = consensus.Length;
            int spacer_skip = 10;

            // Upstream scan
            int index = consensus.Start + k + spacer_skip;
            const int reset = SCAN_DOMAIN;
            int countdown = reset;
            while (countdown-- > 0)
            {
                try
                {
                    Sequence kmer = genome.Substring(index++, k);
                    if (Sequence.Mutant(consensus, kmer))
                    {
                        crispr.AddRepeat(kmer);
                        index = kmer.Start + k + spacer_skip;
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
            index = consensus.Start - k - spacer_skip;
            countdown = reset;
            while (countdown-- > 0)
            {
                try
                {
                    Sequence kmer = genome.Substring(index--, k);
                    if (Sequence.Mutant(consensus, kmer))
                    {
                        crispr.AddRepeat(kmer);
                        index = kmer.Start - k - spacer_skip;
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
            Clusters.Add(new_crispr);
        }

        public void RegisterCrisprs(Crisprs crisprs)
        {
            Clusters.AddRange(crisprs.Clusters);
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

        public void SortByStartPos()
        {
            Clusters.Sort(Crispr.CompareByStart);
        }

        public void MergeCrisprs()
        {
            for (int i = 0; i < Clusters.Count; i++)
            {
                for (int j = 0; j < Clusters.Count; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    Crispr a = Clusters.ElementAt(i);
                    Crispr b = Clusters.ElementAt(j);
                    if (Crispr.Mergeable(a, b))
                    {
                        a.AddRepeats(b.Repeats);
                        Clusters.RemoveAt(j);
                        MergeCrisprs();
                    }
                }
            }
        }

        //public void PrintMutantConsensuses()
        //{
        //    foreach (Crispr crispr in Clusters)
        //    {
        //        foreach (Crispr sub_crispr in Clusters)
        //        {
        //            if (crispr.Equals(sub_crispr))
        //            {
        //                continue;
        //            }
        //            if (Sequence.Mutant(crispr.Consensus, sub_crispr.Consensus, allow_discrepant_lengths:true))
        //            {
        //                Console.WriteLine("Found a mutant.");
        //            }
        //        }
        //    }
        //}

        public static Crisprs DiscoverCrisprs(Sequence genome, int k)
        {
            Crisprs crisprs = new Crisprs();
            int num_kmers = genome.Length - k + 1;
            for (int i = 0; i < num_kmers; i++)
            {
                //Console.WriteLine($"{i} at entry");

                Sequence kmer = genome.Substring(i, k);

                if (!Sequence.Dyad(kmer))
                {
                    continue;
                }

                //Console.WriteLine($"found dyad at {i} of genome len {genome.Length}");

                Crispr crispr = Crispr.DiscoverCrispr(genome, kmer);
                if (crispr != null)
                {
                    Console.WriteLine($"registering CRISPR...");
                    crisprs.RegisterCrispr(crispr);
                    i = crispr.Last.End + 1;
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
