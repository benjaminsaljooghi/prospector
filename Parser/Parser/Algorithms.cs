using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Parser
{
    using System.Collections;
    using static Sequence;

    //public class Dyad : Sequence
    //{
    //    // ?
    //}

    public partial class Sequence : IEnumerable<char>, ICloneable
    {
        const int dyad_min = 5;
        static Dictionary<char, char> complements = new Dictionary<char, char>
        {
            { 'A', 'T' },
            { 'T', 'A' },
            { 'C', 'G' },
            { 'G', 'C' }
        };

        public Sequence Substring(int start, int length)
        {
            return new Sequence(Seq.Substring(start, length), Pos + start);
        }

        public static List<Sequence> Kmers(Sequence seq, int k)
        {
            List<Sequence> kmers = new List<Sequence>();
            int n = seq.Length - k + 1;
            for (int i = 0; i < n; i++)
            {
                kmers.Add(seq.Substring(i, k));
            }
            return kmers;
        }

        public static List<List<Sequence>> KmerWindow(Sequence seq, int k_start, int k_end)
        {
            var kmerWindow = new List<List<Sequence>>();
            for (int k = k_start; k < k_end; k++)
            {
                kmerWindow.Add(Kmers(seq, k));
            }
            return kmerWindow;
        }

        static Sequence ReverseComplement(Sequence seq)
        {
            char[] rc = seq.ToCharArray();
            for (int i = 0; i < rc.Length; i++)
            {
                rc[i] = complements[rc[i]];
            }
            Array.Reverse(rc);
            return new Sequence(rc, seq.Pos);
        }

        static bool Palindrome(Sequence seq)
        {
            return ReverseComplement(seq).Equals(seq);
        }

        public static bool Dyad(Sequence seq)
        {
            for (int i = dyad_min; i < seq.Length / 2; i++)
            {
                Sequence beginning = seq.Substring(0, i);
                Sequence end = seq.Substring(seq.Length - i, i);
                Sequence end_rc = ReverseComplement(end);
                if (beginning.Equals(end_rc))
                {
                    return true;
                }
            }
            return false;
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

    public class Algorithms
    {
        static void Main(string[] args)
        {
            Sequence sequence = new Sequence(@"P:\Honours\sequence.fasta");
            foreach (List<Sequence> kmers in KmerWindow(sequence, 34, 38))
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
    }
}
