
#include "blast.h"

map <string, int> BLAST (set <string> _seqs)
{
    double start = omp_get_wtime();

    vector<string> seqs;
    size_t max = 0;
    for (string seq : _seqs)
    {
        seqs.push_back(seq);
        if (seq.length() > max)
        {
            max = seq.length();
        }
    }

    printf("blasting %zd sequences of max length %zd... ", _seqs.size(), max);


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
    string target_db_path("crispr-data/bacteriophages.fasta");
    const CSearchDatabase target_db(target_db_path, CSearchDatabase::eBlastDbIsNucleotide);

    string fasta = seqs_to_fasta(seqs);
    CBlastFastaInputSource fasta_input(fasta, iconfig);
    CBlastInput blast_input(&fasta_input);
    TSeqLocVector query_loc = blast_input.GetAllSeqLocs(scope);
    CRef<IQueryFactory> query_factory(new CObjMgr_QueryFactory(query_loc));
    CLocalBlast blaster(query_factory, opts, target_db);

    CSearchResultSet results = *blaster.Run();


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

    done(start);
    return seq_max_scores;
}
