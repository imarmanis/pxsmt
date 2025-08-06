#include "llvm/IR/Module.h"
#include "llvm/AsmParser/Parser.h"
#include "llvm/Support/SourceMgr.h"
#include "LLVMModule.hpp"
#include "verifier.hpp"
#include "graph.hpp"
#include "debug.hpp"
#include "config.hpp"
#include "AnalysisInfo.hpp"
#include <ctime>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <typeinfo>

void boost::throw_exception(std::exception const & e){ std::cout << e.what() << "\n"; exit(1); }

std::string exec(std::string cmd) {
    std::array<char, 128> buffer;
    std::string result;
    auto pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        exit_msg("popen() failed!\n");
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
        result += buffer.data();
    auto rc = pclose(pipe);
    if (rc != EXIT_SUCCESS)
        exit(rc);
    return result;
}

#ifdef TRACE
std::ofstream LOGSTREAM;
#endif

#ifdef PROF
std::ofstream PROFSTREAM;
#endif

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    po::options_description desc("Allowed options");
    po::positional_options_description p;
    p.add("file", -1);
    desc.add_options()
        ("help,h", "produce help message")
        ("file,f", po::value<std::string>(), "input ll file")
        ("bound,b", po::value<int>(), "unroll bound")
        ("stats", "Emit statistics")
        ("trace,t", "enable trace")
        ("debug,d", "enable debug")
        ("no-prof,np", "do not collect profiling data")
        ("paths", po::value<unsigned>(), "propagation paths bound")
        ("pairs", po::value<unsigned>(), "unit-edge propagation bound")
        ("tso", "tso memory model")
        ("sc", "sc memory model")
        ("porf-idl", "use IDL for porf paths")
        ("icd", "use ICD instead of ITC")
        ("itc", "use ITC instead of ICD")
        ("delta", "use delta for the ICD")
        ("no-delta", "don't use delta for the ICD")
        ("idl", "only use IDL")
        ("full-crash", "fully encode crash")
        ("guard-crash", "encode crash with guard variables")
        ("all-porf", "use all porf paths")
        ("sanity,sn", "sanity : check that is is SAT w/o asserts")
        ("check-unroll,u", "check loops are fully unrolled")
    ;

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }

    if (vm.count("debug")) debug_enabled = true;
    if (vm.count("trace")) trace_enabled = true;
    if (vm.count("no-prof")) prof_enabled = false;
    if (vm.count("stats")) emit_stats = true;
    if (vm.count("tso")) mmk = TSO;
    if (vm.count("sc")) mmk = SC;
    if (vm.count("icd")) algk = ICD;
    if (vm.count("itc")) algk = ITC;
    if (vm.count("porf-idl")) porf_idl = true;
    if (vm.count("sanity")) sanity_check = true;
    if (vm.count("check-unroll")) check_unroll = true;
    if (vm.count("delta") || vm.count("no-delta")) {
        if (algk != ICD)
            std::cerr << "Delta option is only applicable with ICD algorithm.\n";
        if (vm.count("delta"))
            use_delta = true;
        if (vm.count("no-delta"))
            use_delta = false;
    }
    /* Set both */
    if (vm.count("idl")) idl = true;
    if (vm.count("full-crash")) crash_encoding = FULL;
    if (vm.count("guard-crash")) crash_encoding = GUARD;
    if (vm.count("all-porf")) all_porf_paths = true;

    if (porf_idl && (crash_encoding != PARTIAL))
        exit_msg("porf-idl flag only applicable when the PARTIAL encoding is used.\n");

    if (vm.count("bound"))
        bound = vm["bound"].as<int>();
    if (vm.count("paths"))
        max_paths = vm["paths"].as<unsigned>();
    if (vm.count("pairs"))
        max_pairs = vm["pairs"].as<unsigned>();

    std::stringstream buffer;
    std::string inp_name;
    if (vm.count("file")) {
        const boost::any& v = vm["file"].value();
        auto fname = boost::any_cast<std::string>(&v);   // returns nullptr
        inp_name = *fname;
        buffer << exec("./llvm11-bin/bin/clang-11 -S -emit-llvm -mclflushopt " + inp_name +
            //" -Xclang -enable-trivial-auto-var-init-zero-knowing-it-will-be-removed-from-clang" +
            " -o - -O0 -Xclang -disable-O0-optnone -Iinputs/include");
    } else {
        inp_name = "stdin";
        buffer << std::cin.rdbuf();
    }

    auto inp_str = buffer.str();
    auto mod = LLVMModule::getLLVMModule(inp_str);
    AnalysisInfo analysis;
    LLVMModule::transformLLVMModule(*mod, bound, analysis);

    TRACE_CODE(LLVMModule::printLLVMModule(*mod, "transformed.ll"););

    #if defined(TRACE) || defined(PROF)
    std::time_t result = std::time(nullptr);
    #endif

    TRACE_CODE(
        std::stringstream ss;
        ss << "logs/" << result << ".log";
        LOGSTREAM.open(ss.str(), std::ofstream::out | std::ofstream::trunc);
    );

    PROF_CODE(
        std::stringstream ss;
        ss << "logs/" << result << ".prof.log";
        PROFSTREAM.open(ss.str(), std::ofstream::out | std::ofstream::trunc);
    );

    std::stringstream alg_desc;
    alg_desc << "algk : " << algk;
    if (algk == ICD)
        alg_desc << " (use_delta : " << use_delta << ")";
    PLOG("input : " << inp_name << ", bound = " << bound <<
         ", mm : " << mmk << ", " << alg_desc.str() <<
         ", porf-idl : " << porf_idl << ", crash_encoding : " << crash_encoding <<
         ", paths = " << max_paths << ", pairs = " << max_pairs << "\n");


    Verifier ver(std::move(mod), analysis);
    return 0;
}
