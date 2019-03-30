using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;


using Newtonsoft.Json;

namespace Parser
{
    public class Program
    {
        public const string DIR = @"P:\Honours\";

        public static void Write(string content, string file_name, string dir = DIR)
        {
            File.WriteAllText(Path.Combine(dir, file_name), content);
        }

        public static string Read(string file_name, string dir = DIR)
        {
            return File.ReadAllText(Path.Combine(dir, file_name));
        }

        public static void Serialize(object obj, string file_name, string dir = DIR)
        {
            Write(JsonConvert.SerializeObject(obj), file_name, dir);
        }

        public static T Deserialize<T>(string file_name, string dir = DIR)
        {
            return JsonConvert.DeserializeObject<T>(Read(file_name, dir));
        }

        public static void Main()
        {
            Sequence cas9 = new Sequence("cas9.fasta");
            Sequence pyogenes = new Sequence("pyogenes.fasta");
            Sequence aureus = new Sequence("aureus.fasta");

            //List<Sequence> dyads = Sequence.FrequencyOrderedDyads(aureus, Crispr.REPEAT_MIN, Crispr.REPEAT_MAX);
            List<Sequence> dyads = Deserialize<List<Sequence>>("aureus_dyads.json");

            Crisprs crisprs = Crisprs.DiscoverCrisprs(aureus, dyads);
            Console.WriteLine(crisprs);

            Serialize(dyads, "aureus_dyads.json");
            Serialize(crisprs, "aureus_crisprs.json");

            Console.WriteLine("Press any key to quit.");
            Console.ReadKey();
        }

    }
}