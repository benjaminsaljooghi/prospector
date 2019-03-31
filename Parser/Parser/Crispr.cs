﻿using Newtonsoft.Json;
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

        public Crispr(Sequence consensus)
        {
            Consensus = consensus;
            Repeats = new List<int>();
        }

        public Sequence Consensus { get; private set; }

        public List<int> Repeats { get; }

        public void AddRepeat(int repeat_index)
        {
            Repeats.Add(repeat_index);
        }

        public void UpdateConsensus()
        {
            throw new NotImplementedException();
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
