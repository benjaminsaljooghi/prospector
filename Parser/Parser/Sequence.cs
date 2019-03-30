using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Parser
{
    public partial class Sequence : IEnumerable<char>, ICloneable
    {
        const string dir = @"P:\Honours\";

        public string Seq { get; }
        public int Length { get { return Seq.Length; } }
        public int Pos { get; }

        public Sequence(string sequence, int pos)
        {
            Seq = sequence;
            Pos = pos;
        }

        public Sequence(char[] sequence, int pos) : this(new string(sequence), pos)
        {

        }

        public Sequence(string file)
        {
            var reader = new StreamReader(Path.Combine(dir, file));
            while (reader.ReadLine().StartsWith(">")) ;
            Seq = reader.ReadToEnd().Replace("\n", "");
            Pos = 0;
        }

        public static implicit operator string(Sequence fasta)
        {
            return fasta.Seq;
        }

        public override bool Equals(object obj)
        {
            return Seq == (obj as Sequence).Seq;
        }

        public override int GetHashCode()
        {
            return Seq.GetHashCode();
        }

        public override string ToString()
        {
            return Seq;
        }

        public IEnumerator<char> GetEnumerator()
        {
            return Seq.GetEnumerator();
        }

        IEnumerator IEnumerable.GetEnumerator()
        {
            return Seq.GetEnumerator();
        }

        public Sequence Clone()
        {
            return new Sequence(Seq, Pos);
        }

        object ICloneable.Clone()
        {
            return Clone();
        }

        private char[] ToCharArray()
        {
            return Seq.ToCharArray();
        }
    }







    public partial class Sequence : IEnumerable<char>, ICloneable
    {
        const int dyad_min = 5;
        static Dictionary<char, char> complements = new Dictionary<char, char>
        {
            { 'A', 'T' },
            { 'T', 'A' },
            { 'C', 'G' },
            { 'G', 'C' },
            { 'N', 'N' }
        };


        public static bool Mutant(Sequence a, Sequence b)
        {
            if (a.Length != b.Length)
            {
                throw new Exception("Mutants must be the same length.");
            }

            int allowed_point_mutations = a.Length / 10;
            int point_mutations = 0;

            for (int i = 0; i < a.Length; i++)
            {
                if (a.Seq[i] != b.Seq[i] && ++point_mutations > allowed_point_mutations)
                {
                    return false;
                }
            }
            return true;
        }

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
                char rc_char = rc[i];
                rc[i] = complements[rc_char];
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
            for (int i = dyad_min; i <= seq.Length / 2; i++)
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


}
