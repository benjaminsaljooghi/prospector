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
    string target_db_path("/home/benjamin/proj/data/bacteriophages.fasta");
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


    vector<string> all_spacers;


    for (Util::Locus crispr : prospection.crisprs)
    {
        vector<string> spacers = Util::spacers(prospection.genome, crispr);
        all_spacers.insert(all_spacers.end(), spacers.begin(), spacers.end());
    }


    map<string, int> spacer_scores = BLAST(all_spacers);

    for (Util::Locus crispr : prospection.crisprs)
    {
        vector<string> repeats = Util::repeats(prospection.genome, crispr);
        vector<string> spacers = Util::spacers(prospection.genome, crispr);

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
