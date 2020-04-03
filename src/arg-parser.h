#ifndef ARGPARSER_H
#define ARGPARSER_H

#include "args.h"
#include "input-output.h"
#include <iostream>
#include <unordered_map>
#include <string>

using namespace std;

struct argparser{

    argparser(int argc, char* argv[]);
    unordered_map<string, string> values;
    string get(string key, string default_value);

    private:
        template<class T> void set_from_obj(args::ValueFlag<T> & opt, string key, T default_value);

};


argparser::argparser(int argc, char* argv[]){
    
    args::ArgumentParser parser("This is the Hyplex code.", "This code runs with MPI, possibly causing other flags to be passed here.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    
    args::ValueFlag<string> config_path(parser, "config", "Path to configuration .yaml file", {'c', "config"});
    args::ValueFlag<string> job_name(parser, "name", "Override job identifier name", {'n', "name"});
    
    try
    {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
        cout << parser;
        exit(0);
    }
    catch (args::ParseError e)
    {
        std::cerr << e.what() << std::endl;
        // std::cerr << parser;
        // exit(0);
    }
    catch (args::ValidationError e)
    {
        std::cerr << e.what() << std::endl;
        // std::cerr << parser;
        // exit(0);
    }

    string default_val = "";
    set_from_obj(config_path, "config", default_val);
    set_from_obj(job_name, "name", default_val);
}

template <class T>
void argparser::set_from_obj(args::ValueFlag<T> & opt, string key, T default_value){
    values[key] = opt ? args::get(opt) : default_value;
    if(opt) verbose_log( "Setting option " + key + ": " +  values[key]);
}

string argparser::get(string key, string default_value){
    string val = values[key];
    return val == "" ? default_value : val;
}


#endif