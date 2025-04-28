/**
 * inspired from https://github.com/ksahlin/strobealign/blob/main/src/cmdline.cpp
*/

#include "cmdline.hpp"
#include "args.hxx"

CommandLineOptions parse_cla_segment_count(int argc, char **argv) {
    args::ArgumentParser parser("Segment Count" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "Segment Count";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters. [Required]", {'g', "genome_fasta"},args::Options::Required);
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file. [Required]", {'s', "suffix_array"},args::Options::Required);
    args::ValueFlag<int64_t> kmer_size(parser, "INT", "Kmer size to be used to construct the index. [default: 21]", {'k', "kmer_size"});
    args::Flag use_all_epsilons(parser, "BOOL", "Use consecutive epsilon values rather than powers of two (up to eps 1024) [default: false]", {'a', "all_eps"} );
    args::ValueFlag<std::string>count_fn(parser, "STRING", "File name where to save the segment count. [default: genome_fasta.rourk.woshare.out]", {'o', "out"});
    
    
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    if(genome_fasta) {opt.gn_fn = args::get(genome_fasta);}
    if(suffix_array) {opt.sa_fn = args::get(suffix_array);}
    if(kmer_size) {opt.kmer_size = args::get(kmer_size);}
    if(count_fn) {opt.count_fn = args::get(count_fn);}
    if(use_all_epsilons) {opt.use_all_epsilons = true;}
    if(opt.count_fn == "-1"){
        opt.count_fn = opt.gn_fn+".segments.txt";
    }

    return opt;
}

/*
CommandLineOptions parse_cla_kmer_type(int argc, char **argv) {
    args::ArgumentParser parser("kmer_type" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "kmer_type";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters.", 
            {'g', "genome_fasta"}, args::Options::Required);
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file.", 
            {'s', "suffix_array"}, args::Options::Required);
    
    args::ValueFlag<std::string>random_genome_fasta(parser, "STRING", "Random genome fasta file with one entry and only ACGT characters.", 
            {'r', "genome_fasta"}, args::Options::Required);
    args::ValueFlag<std::string>random_suffix_array(parser, "STRING", "Random genome suffix array of the genome in a binary file.", 
            {'t', "suffix_array"}, args::Options::Required);
    args::ValueFlag<std::string>fn_prefix(parser, "STRING", "Prefix to be used for saving the genome types.", {'p', "prefix"},
            args::Options::Required);
    
    args::ValueFlag<int64_t> kmer_size(parser, "INT", "Kmer size to be used. [default: 21]", {'k', "kmer_size"});
    
    try {
        parser.ParseCLI(argc, argv);
    }
    catch (const args::Completion& e) {
        std::cout << e.what();
        exit(EXIT_SUCCESS);
    }
    catch (const args::Help&) {
        std::cout << parser;
        exit(EXIT_SUCCESS);
    }
    catch (const args::Error& e) {
        std::cerr << parser;
        std::cerr << "Error: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    CommandLineOptions opt;

    if(genome_fasta) {opt.gn_fn = args::get(genome_fasta);}
    if(suffix_array) {opt.sa_fn = args::get(suffix_array);}
    if(random_genome_fasta) {opt.rn_gn_fn = args::get(random_genome_fasta);}
    if(random_suffix_array) {opt.rn_sa_fn = args::get(random_suffix_array);}
    if(fn_prefix) {opt.fn_prefix = args::get(fn_prefix);}
    
    if(kmer_size) {opt.kmer_size = args::get(kmer_size);}
    

    return opt;
}
*/