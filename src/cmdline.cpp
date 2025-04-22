/**
 * inspired from https://github.com/ksahlin/strobealign/blob/main/src/cmdline.cpp
*/

#include "cmdline.hpp"
#include "args.hxx"

CommandLineOptions parse_cla_count_segments(int argc, char **argv) {
    args::ArgumentParser parser("count_segments" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "count_segments";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters. [Required if point seq]", {'g', "genome_fasta"});
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file. [Required if point seq]", {'s', "suffix_array"});
    args::ValueFlag<int64_t> kmer_size(parser, "INT", "Kmer size to be used to construct the index. [default: 21]", {'k', "kmer_size"});
    args::Flag use_all_epsilons(parser, "BOOL", "Use consecutive epsilon values rather than powers of two (up to eps 1024) [default: false]", {'a', "all_eps"} );


    args::ValueFlag<std::string>in_fn(parser, "STRING", "Input File with Integer data in binary format.", {'i', "int"});
    args::Flag is_dict_point_seq(parser, "BOOL", "Boolean that says whether the points are to be considered as dictionary point set [default: false]", {'d', "is_dict_points"} );

    args::Flag is_sketch_data(parser, "BOOL", "Boolean that says whether the points are sketch data (given file is then fracminhash sketch) [default: false]", {'S', "is_sketch"} );
    args::Flag is_int_vec(parser, "BOOL", "Boolean that says whether the points are from integer vec [default: false]", {'I', "is_int_vec"} );

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
    if(is_dict_point_seq){opt.is_dict_point_seq = true;}

    if(is_sketch_data){
        opt.is_genome_data = false;
        opt.is_sketch_data = true;
    }
    if(is_int_vec){
        opt.is_genome_data = false;
        opt.is_int_vec = true;
    }
    if(in_fn){opt.in_fn = args::get(in_fn);}
    

    if(opt.count_fn == "-1"){
        if(is_sketch_data){
            opt.count_fn = opt.in_fn+".sketch.rourk.txt";
        }
        else if(is_dict_point_seq)
            opt.count_fn = opt.gn_fn+".dict.rourk.txt";
        else if(is_int_vec)
            opt.count_fn = opt.in_fn+".index.rourk.txt";
        else
            opt.count_fn = opt.gn_fn+".index.rourk.txt";
    }       

    return opt;
}

CommandLineOptions parse_cla_frac_sketch(int argc, char **argv) {
    args::ArgumentParser parser("create_fracminhash_sketch" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "create_fracminhash_sketch";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>genome_fasta(parser, "STRING", "Fasta file with one entry and only ACGT characters. [Required]", {'g', "genome_fasta"},args::Options::Required);
    args::ValueFlag<std::string>suffix_array(parser, "STRING", "Suffix array of the genome in a binary file. [Required]", {'s', "suffix_array"},args::Options::Required);
    args::ValueFlag<int64_t> kmer_size(parser, "INT", "Kmer size to be used to construct the index. [default: 21]", {'k', "kmer_size"});
    args::ValueFlag<double> theta(parser, "DOUBLE", "Theta value for sketching [default: 0.1]", {'t', "theta"});
    args::ValueFlag<std::string>sketch_name(parser, "STRING", "File name where to save the sketch (kmers). [default: genome_fasta.sketch.k_kmer.t_theta]", {'o', "sketch"});
    
    
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
    if(sketch_name) {opt.in_fn = args::get(sketch_name);}
    if(theta) {opt.theta = args::get(theta);}
    

    if(opt.in_fn == "-1"){
        opt.in_fn = opt.gn_fn+".sketch.k_"+std::to_string(opt.kmer_size)+
            ".t_"+std::to_string(opt.theta);
    }

    return opt;
}

CommandLineOptions parse_cla_space_comp(int argc, char **argv){
    args::ArgumentParser parser("space_comparison" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "space_comparison";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>sketch_name(parser, "STRING", "File name where sketch (kmers) is saved. [Required]", {'f', "sketch"}, args::Options::Required);
    args::ValueFlag<std::string>out_fn(parser, "STRING", "File name where space information is written (appended). [default: sketch_name.space.out]", {'o', "file"});
    
    
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

    if(sketch_name) {opt.in_fn = args::get(sketch_name);}
    if(out_fn) {opt.out_fn = args::get(out_fn);}
    

    if(opt.out_fn == "-1"){
        opt.out_fn = opt.in_fn+".space.out";
    }

    return opt;
}

CommandLineOptions parse_cla_pinch_point(int argc, char **argv) {
    args::ArgumentParser parser("find_pinch_point" );
    parser.helpParams.showTerminator = false;
    parser.helpParams.helpindent = 20;
    parser.helpParams.width = 90;
    parser.helpParams.programName = "find_pinch_point";
    parser.helpParams.shortSeparator = " ";

    args::HelpFlag help(parser, "help", "Print help and exit", {'h', "help"});

    // build_index
    args::ValueFlag<std::string>segment_file(parser, "STRING", "eps vs segment file. [Required]", {'s', "segment"}, args::Options::Required);
    args::ValueFlag<double>increment(parser, "DOUBLE", "increment value to use [default: 1E-6]", {'i', "increment"});
    
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

    if(segment_file) {opt.segment_file = args::get(segment_file);}
    if(increment) {opt.increment = args::get(increment);}
    
    

    return opt;
}

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