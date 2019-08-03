using Newtonsoft.Json;
using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Prospector
{
    public partial class Sequence
    {
        public string Seq { get; private set; }

        public int Start { get; private set; }

        [JsonIgnore]
        public int End { get { return Start + Seq.Length - 1; }  }

        [JsonIgnore]
        public int Length { get { return Seq.Length; } }

        [JsonConstructor]
        public Sequence(string seq, int start)
        {
            Seq = seq;
            Start = start;
        }

        public Sequence(char[] sequence, int start) : this(new string(sequence), start)
        {
        }

        public Sequence(string file)
        {
            var reader = new StreamReader(file);
            while (reader.ReadLine().StartsWith(">")) ;
            Seq = reader.ReadToEnd().Replace("\n", "");
            Start = 0;
        }

        public static implicit operator string(Sequence fasta)
        {
            return fasta.Seq;
        }

        public override bool Equals(object obj)
        {
            Sequence seq = obj as Sequence;
            return Seq == seq.Seq && Start == seq.Start;
        }

        public override int GetHashCode()
        {
            return Seq.GetHashCode() * 7 + Start.GetHashCode() * 7;
        }

        public override string ToString()
        {
            return Seq;
        }

        public IEnumerator<char> GetEnumerator()
        {
            return Seq.GetEnumerator();
        }

        public Sequence Clone()
        {
            return new Sequence(Seq, Start);
        }

        private char[] ToCharArray()
        {
            return Seq.ToCharArray();
        }

        public Sequence Substring(int start, int length)
        {
            return new Sequence(Seq.Substring(start, length), Start + start);
        }

        public static int CompareStart(Sequence a, Sequence b)
        {
            return a.Start - b.Start;
        }
    }

    public partial class Sequence
    {
        const int DYAD_MIN = 5;
        static Dictionary<char, char> complements = new Dictionary<char, char>
        {
            { 'A', 'T' },
            { 'T', 'A' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'N', 'N' }
        };

        public static HashSet<Sequence> Dyads(Sequence genome, int k)
        {
            HashSet<Sequence> dyads = new HashSet<Sequence>();
            int num_kmers = genome.Length - k + 1;
            for (int i = 0; i < num_kmers; i++)
            {
                Sequence kmer = genome.Substring(i, k);
                if (Dyad(kmer))
                {
                    dyads.Add(kmer);
                }
            }
            return dyads;
        }

        public static HashSet<Sequence> Dyads(Sequence genome, int k_begin, int k_end)
        {
            HashSet<Sequence> dyads = new HashSet<Sequence>();
            for (int k = k_begin; k <= k_end; k++)
            {
                dyads.UnionWith(Dyads(genome, k));
            }
            return dyads;
        }

        //public static Dictionary<Sequence, int> DyadFrequencies(Sequence genome, int k)
        //{
        //    Dictionary<Sequence, int> frequencies = new Dictionary<Sequence, int>();
        //    int num_kmers = genome.Length - k + 1;
        //    for (int i = 0; i < num_kmers; i++)
        //    {
        //        Sequence kmer = genome.Substring(i, k);
        //        if (!Sequence.Dyad(kmer))
        //        {
        //            continue;
        //        }
        //        if (!frequencies.ContainsKey(kmer))
        //        {
        //            frequencies.Add(kmer, 0);
        //        }
        //        frequencies[kmer] += 1;
        //    }
        //    return frequencies;
        //}

        //public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Sequence genome, int k_begin, int k_end)
        //{
        //    List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = new List<KeyValuePair<Sequence, int>>();
        //    for (int k = k_begin; k <= k_end; k++)
        //    {
        //        ordered_dyad_frequencies.AddRange(OrderedDyadFrequencies(genome, k));
        //    }
        //    return ordered_dyad_frequencies;
        //}

        //public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Sequence genome, int k)
        //{
        //    return OrderedDyadFrequencies(DyadFrequencies(genome, k));
        //}

        //public static List<KeyValuePair<Sequence, int>> OrderedDyadFrequencies(Dictionary<Sequence, int> dyad_frequencies)
        //{
        //    List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = dyad_frequencies.ToList();
        //    ordered_dyad_frequencies.Sort((a, b) => b.Value - a.Value);
        //    return ordered_dyad_frequencies;
        //}

        //public static List<Sequence> FrequencyOrderedDyads(Sequence genome, int k_begin, int k_end)
        //{
        //    List<KeyValuePair<Sequence, int>> ordered_dyad_frequencies = OrderedDyadFrequencies(genome, k_begin, k_end);
        //    List<Sequence> consensuses = new List<Sequence>();
        //    foreach (KeyValuePair<Sequence, int> consensus in ordered_dyad_frequencies)
        //    {
        //        Console.WriteLine("Adding consensus {0} with global frequency {1}", consensus.Key, consensus.Value);
        //        consensuses.Add(consensus.Key);
        //    }
        //    return consensuses;
        //}

        public static bool Mutant(Sequence a, Sequence b, bool allow_discrepant_lengths = false)
        {
            if (!allow_discrepant_lengths && a.Length != b.Length)
            {
                throw new Exception("Mutants must be the same length.");
            }

            int len = Math.Min(a.Length, b.Length);

            int allowed_point_mutations = a.Length / 10;
            int point_mutations = 0;

            for (int i = 0; i < len; i++)
            {
                if (a.Seq[i] != b.Seq[i] && ++point_mutations > allowed_point_mutations)
                {
                    return false;
                }
            }
            return true;
        }

        //public static List<Sequence> Kmers(Sequence seq, int k)
        //{
        //    List<Sequence> kmers = new List<Sequence>();
        //    int n = seq.Length - k + 1;
        //    for (int i = 0; i < n; i++)
        //    {
        //        kmers.Add(seq.Substring(i, k));
        //    }
        //    return kmers;
        //}

        //public static List<List<Sequence>> KmerWindow(Sequence seq, int k_start, int k_end)
        //{
        //    var kmerWindow = new List<List<Sequence>>();
        //    for (int k = k_start; k < k_end; k++)
        //    {
        //        kmerWindow.Add(Kmers(seq, k));
        //    }
        //    return kmerWindow;
        //}

        //public static Sequence ReverseComplement(Sequence seq)
        //{
        //    char[] rc = seq.ToCharArray();
        //    for (int i = 0; i < rc.Length; i++)
        //    {
        //        char rc_char = rc[i];
        //        rc[i] = complements[rc_char];
        //    }
        //    Array.Reverse(rc);
        //    return new Sequence(rc, seq.Start);
        //}

        //public static bool DyadCheck(Sequence seq, int dyad_min = DYAD_MIN)
        //{
        //    Sequence beginning = seq.Substring(0, dyad_min);
        //    Sequence end = seq.Substring(seq.Length - dyad_min, dyad_min);
        //    Sequence end_rc = ReverseComplement(end);
        //    return beginning.Equals(end_rc);

        //    //for (int i = dyad_min; i <= seq.Length / 2; i++)
        //    //{
        //    //    Sequence beginning = seq.Substring(0, i);
        //    //    Sequence end = seq.Substring(seq.Length - i, i);
        //    //    Sequence end_rc = ReverseComplement(end);
        //    //    if (beginning.Equals(end_rc))
        //    //    {
        //    //        return true;
        //    //    }
        //    //}
        //    //return false;
        //}

        public static bool Dyad(string seq, int dyad_min = DYAD_MIN)
        {
            int len = seq.Length;
            for (int i = 0; i < dyad_min; i++)
            {
                char beginning_upstream = seq[i];
                char end_downstream = seq[len - i - 1];
                char end_downstream_comp = complements[end_downstream];
                if (beginning_upstream != end_downstream_comp)
                {
                    return false;
                }
            }
            return true;
        }


        public static List<Sequence> Dyads(List<Sequence> kmers)
        {
            List<Sequence> dyads = new List<Sequence>();
            foreach (Sequence kmer in kmers)
            {
                if (Dyad(kmer))
                {
                    dyads.Add(kmer);
                }
            }
            return dyads;
        }

    }


}
