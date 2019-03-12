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
            var reader = new StreamReader(file);
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
            return Seq == obj as Sequence;
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
}
