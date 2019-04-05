﻿using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

using Newtonsoft.Json;

using static System.IO.File;

namespace Parser
{
    public class Program
    {
        public const string DIR = @"P:\Honours\";

        public static void Serialize(string path, object obj)
        {
            WriteAllText(path, JsonConvert.SerializeObject(obj));
        }

        public static T Deserialize<T>(string path)
        {
            return JsonConvert.DeserializeObject<T>(ReadAllText (path));
        }

        public static void Execution(string bacterium_path)
        {
            Sequence bacterium = new Sequence(bacterium_path);
            //Crisprs crisprs = Crisprs.DiscoverCrisprs(bacterium, Crispr.REPEAT_MIN, Crispr.REPEAT_MAX);
            //Serialize(Path.Combine(DIR, "aureus_crisprs.json"), crisprs);
            Crisprs crisprs = Deserialize<Crisprs>(Path.Combine(DIR, "aureus_crisprs.json"));

            Console.WriteLine("Before merging crisprs:");
            //crisprs.SortByStartPos();
            Console.WriteLine(crisprs);


            Console.WriteLine("After merging crisprs:");

            crisprs.MergeCrisprs();

            crisprs.SortByStartPos();
            Console.WriteLine(crisprs); ;
        }

        public static void Main()
        {
            string bacterium_path = Path.Combine(DIR, "aureus.fasta");
            Execution(bacterium_path);
            Quit();
        }

        public static void Quit()
        {
            Console.WriteLine("Press any key to quit.");
            Console.ReadKey();
        }


    }
}