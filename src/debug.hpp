#ifndef _DEBUG_HPP
#define _DEBUG_HPP

#include <assert.h>
#include <fstream>

void exit_msg(const char msg[]);

extern bool debug_enabled;
extern bool trace_enabled;
extern bool prof_enabled;

#ifdef DEBUG
#define DEBUG_CODE(CODE) { if (debug_enabled){CODE} } ((void) 0)
#else
#define DEBUG_CODE(CODE) ((void) 0)
#endif

#ifdef TRACE
#define LOGSTREAM logfile
extern std::ofstream LOGSTREAM;
#define TRACE_CODE(CODE) { if (trace_enabled){CODE} } ((void) 0)
#else
#define TRACE_CODE(CODE) ((void) 0)
#endif

#ifdef PROF
#define PROFSTREAM proffile
extern std::ofstream PROFSTREAM;
#define PROF_CODE(CODE) { if (prof_enabled){CODE} } ((void) 0)
#else
#define PROF_CODE(CODE) ((void) 0)
#endif

#define SASSERT(COND) DEBUG_CODE(assert(COND);)

#define LOG(MSG) TRACE_CODE(LOGSTREAM << MSG; LOGSTREAM.flush();)
#define PLOG(MSG) PROF_CODE(PROFSTREAM << MSG; PROFSTREAM.flush();)

#endif
