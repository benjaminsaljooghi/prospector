using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Parser
{
    class Program
    {
        static void FindCRISPR(string sequence)
        {

        }

        static string Substring(string sequence, int start, int end)
        {
            return sequence.Substring(start, end - start);
        }

        static string ReadFASTA(string file)
        {
            var reader = new StreamReader(file);
            while (true)
            {
                var line = reader.ReadLine();
                if (line.StartsWith(">"))
                {
                    break;
                }
            }
            return reader.ReadToEnd().Replace("\n", "");
        }

        static string ReverseComplement(string sequence)
        {
            char[] chars = sequence.ToCharArray();
            for (int i = 0; i < chars.Length; i++)
            {
                switch (chars[i])
                {
                    case 'A':
                        chars[i] = 'T';
                        break;
                    case 'C':
                        chars[i] = 'G';
                        break;
                    case 'G':
                        chars[i] = 'C';
                        break;
                    case 'T':
                        chars[i] = 'A';
                        break;
                    default:
                        break;
                }
            }
            return new string(chars);
        }

        static string Reverse(string sequence)
        {
            char[] chars = sequence.ToCharArray();
            Array.Reverse(chars);
            return new string(chars);
        }

        static bool Palindrome(string sequence)
        {
            return Reverse(ReverseComplement(sequence)).Equals(sequence);
        }

        static void Main(string[] args)
        {
            Console.WriteLine(Palindrome("ACCTAGGT"));
            //string sequence = ReadFASTA(@"P:\Honours\sequence.fasta");
            //Console.WriteLine("length: {0:n0}", sequence.Length);
            //string cas9 = Substring(sequence, 854751, 858857);
            //Console.WriteLine("cas9, len: " + cas9.Length + "\n" + cas9);
            //string crispr = Substring(sequence, 860819, 861250);
            //Console.WriteLine("cas9: ", cas9);
            //Console.WriteLine("\ncrispr: ", crispr);

        }
    }
}
