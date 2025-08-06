extern int bound;
extern unsigned max_pairs;
extern unsigned max_paths;
extern bool emit_stats;
extern bool use_delta;
extern bool porf_idl;
extern bool idl;
extern bool all_porf_paths;
extern bool check_pers;
extern bool sanity_check;
extern bool check_unroll;

enum CE_K {
    PARTIAL,
    FULL,
    GUARD
};

enum MM_K {
    SC,
    TSO
};

enum ALG_K {
    ITC,
    ICD
};

extern CE_K crash_encoding;
extern MM_K mmk;
extern ALG_K algk;
