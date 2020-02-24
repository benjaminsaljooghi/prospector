/*  $Id: blast_demo.cpp 494083 2016-03-03 16:46:54Z gouriano $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Tom Madden
 *
 * File Description:
 *   Sample application for the running a blast search.
 *
 */



#include <ncbi_pch.hpp>
#include <corelib/ncbiapp.hpp>
#include <corelib/ncbienv.hpp>
// #include <corelib/ncbiargs.hpp>

#include <objmgr/object_manager.hpp>

#include <objects/seqalign/Seq_align_set.hpp>

#include <algo/blast/api/sseqloc.hpp>
#include <algo/blast/api/local_blast.hpp>
#include <algo/blast/api/uniform_search.hpp>
#include <algo/blast/api/blast_types.hpp>
#include <algo/blast/api/blast_aux.hpp>
#include <algo/blast/api/objmgr_query_data.hpp>
#include <algo/blast/api/blast_options_handle.hpp>
#include <algo/blast/api/blast_nucl_options.hpp>

#include <algo/blast/blastinput/blast_input.hpp>
#include <algo/blast/blastinput/blast_fasta_input.hpp>


#include "../prospector/stdafx.h"
#include "../prospector/util.h"
#include "../prospector/prospector.h"


USING_NCBI_SCOPE;
USING_SCOPE(blast);


/////////////////////////////////////////////////////////////////////////////
//  CBlastDemoApplication::


class CBlastDemoApplication : public CNcbiApplication
{
private:
    virtual void Init(void);
    virtual int  Run(void);
    virtual void Exit(void);
};


/////////////////////////////////////////////////////////////////////////////
//  Init test for all different types of arguments


void CBlastDemoApplication::Init(void)
{

}



/////////////////////////////////////////////////////////////////////////////
//  Run test (printout arguments obtained from command-line)




map<string, int> BLAST(vector<string> seqs)
{

    printf("Instantiating a BLAST program...\n");
    clock_t start = clock();

    EProgram program = ProgramNameToEnum("blastn"); 
    CRef<CBlastOptionsHandle> opts(CBlastOptionsFactory::Create(program));
    
    CRef<CObjectManager> objmgr = CObjectManager::GetInstance();
    if (!objmgr) {
         throw std::runtime_error("Could not initialize object manager");
    }
    const bool is_protein = !!Blast_QueryIsProtein(opts->GetOptions().GetProgramType());
    SDataLoaderConfig dlconfig(is_protein);
    CBlastInputSourceConfig iconfig(dlconfig);
    CScope scope(*objmgr);
    string target_db_path("/home/ben/Documents/crispr-data/bacteriophages.fasta");
    const CSearchDatabase target_db(target_db_path, CSearchDatabase::eBlastDbIsNucleotide);

    string fasta = Util::seqs_to_fasta(seqs);
    CBlastFastaInputSource fasta_input(fasta, iconfig);
    CBlastInput blast_input(&fasta_input);
    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);
    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));
    CLocalBlast blaster(query_factory, opts, target_db);

	printf("BLAST program instantiated in %.3f seconds.\n", Util::duration(start));
    printf("Running BLAST program...\n");
    start = clock();

    CSearchResultSet results = *blaster.Run();
    
    printf("BLAST program completed in %.3f seconds.\n", Util::duration(start));

    map<string, int> seq_max_scores;
    for (unsigned int i = 0; i < results.GetNumResults(); i++)
    {
        CConstRef<CSeq_align_set> sas = results[i].GetSeqAlign();
        int max_score = 0;
        for (auto alignment : sas->Get())
        {
            int score;
            ncbi::objects::CSeq_align::EScoreType score_type = ncbi::objects::CSeq_align::EScoreType::eScore_IdentityCount;
            alignment->GetNamedScore(score_type, score);
            if (score > max_score)
            {
                max_score = score;
            }
        }
        seq_max_scores[seqs[i]] = max_score;
    }
    return seq_max_scores;
}





int CBlastDemoApplication::Run(void)
{
    Util::Prospection prospection = Prospector::prospector_main();

    vector<Util::Locus> crisprs = prospection.crisprs;
    string genome = prospection.genome;
    
    


    vector<string> all_spacers;
    for (Util::Locus crispr : crisprs)
    {
        vector<string> spacers = Util::spacers(genome, crispr);
        all_spacers.insert(all_spacers.end(), spacers.begin(), spacers.end());
    }


    map<string, int> spacer_scores = BLAST(all_spacers);



    
    vector<Util::Locus> crisprs_filtered = vector<Util::Locus>();

    for (Util::Locus crispr : crisprs)
    {

        vector<string> spacers = Util::spacers(genome, crispr);

        float spacer_identity_percent_sum = 0;
        for (string spacer : spacers)
        {
            spacer_identity_percent_sum += spacer_scores[spacer] / spacer.length();
        }

        float mean_identity = spacer_identity_percent_sum / spacers.size(); 
    
        if (mean_identity >= 0.5)
        {
            crisprs_filtered.push_back(crispr);
        }
    }


    for (Util::Locus crispr : crisprs_filtered)
    {
        vector<string> repeats = Util::repeats(genome, crispr);
        vector<string> spacers = Util::spacers(genome, crispr);

        std::ostringstream string_stream;
    
        string_stream << crispr.genome_indices[0] << " " << crispr.k << endl;

        string_stream << "\t" << "repeats:" << endl;
        for (unsigned int i = 0; i < spacers.size(); i++)
        {
            string_stream << "\t\t" << crispr.genome_indices[i] << " " << repeats[i] << " " << crispr.genome_indices[i] + crispr.k << endl;
        }

        string_stream << endl;
    
        string_stream << "\t" << "spacers:" << endl;

        for (string spacer : spacers)
        {
            string_stream << "\t\t" << spacer_scores[spacer] << "/" << spacer.length() << " " << spacer << endl; 
        }

        printf("%s\n\n", string_stream.str().c_str());

    }


    // now check upstream of the arrays for cas genes

    int k = 20;
    int upstream_size = 10000;


    // get Cas9
    string cas9_seq = Util::parse_fasta("../crispr-data/addGeneSpCas9.fasta").begin()->second;


    vector<string> cas9_kmers = get_kmers(cas9_seq, k);
    


    string crisrpcasfinder_cas9 = genome.substr(854751, 858857 - 854751);
    vector<string> crisprcasfinder_cas9_kmers = get_kmers(crisrpcasfinder_cas9, k);

    for (Util::Locus crispr : crisprs_filtered)
    {
        // transform the upstream 10k into kmers
        string upstream = genome.substr(crispr.genome_indices[0] - upstream_size, upstream_size);
        vector<string> upstream_kmers = get_kmers(upstream, k);

        // where, or do, these kmers overlap? (comparing the set of the cas9_kmers against the set of the upstream kmers)

        // what quantity of kmers in the cas9_kmers exists in the upstream_kmers?

        int present = 0;

        // for (string cas9_kmer : cas9_kmers)
        for (string cas9_kmer : crisprcasfinder_cas9_kmers)
        {
            for (string upstream_kmer : upstream_kmers)
            {
                present += cas9_kmer.compare(upstream_kmer) == 0 ? 1 : 0;
            }
        }

        printf("for k size %d we have a cas9_kmer count of %d and a presence of %d\n", crispr.k, crisprcasfinder_cas9_kmers.size(), present);

    }
    




    return 0;
}






/////////////////////////////////////////////////////////////////////////////
//  Cleanup


void CBlastDemoApplication::Exit(void)
{
    // Do your after-Run() cleanup here
}


/////////////////////////////////////////////////////////////////////////////
//  MAIN


#ifndef SKIP_DOXYGEN_PROCESSING
int NcbiSys_main(int argc, ncbi::TXChar* argv[])
{
    CBlastDemoApplication().AppMain(argc, argv);
    return 0;
}
#endif /* SKIP_DOXYGEN_PROCESSING */
