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
        



        public static void FindCRISPR(string sequence)
        {

        }

        public static string Substring(string sequence, int start, int end)
        {
            return sequence.Substring(start, end - start);
        }

        public static string ReadFASTA(string file)
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

        static void Main(string[] args)
        {
            string sequence = ReadFASTA(@"P:\Honours\sequence.fasta");
            Console.WriteLine("length: {0:n0}", sequence.Length);
            string cas9 = Substring(sequence, 854751, 858857);
            Console.WriteLine("cas9, len: " + cas9.Length + "\n" + cas9);
            //string crispr = Substring(sequence, 860819, 861250);
            //Console.WriteLine("cas9: ", cas9);
            //Console.WriteLine("\ncrispr: ", crispr);

        }
    }
}
