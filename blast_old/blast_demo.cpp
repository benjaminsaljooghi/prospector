
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



// #include "../prospector/stdafx.h"
#include "../util/util.h"
// #include "../prospector/prospector.h"


USING_NCBI_SCOPE;
USING_SCOPE(blast);




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


// int main()
// {
    
// }
