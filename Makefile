CXX = clang++
CXXFLAGS = -Wall -Wextra -Wno-unused-parameter -Wimplicit-fallthrough -O2
LLVM_CXXFLAGS = `llvm-config-11 --cxxflags`
LLVM_LDFLAGS = `llvm-config-11 --ldflags --libs --system-libs`
RM = rm -f
# Sets Z3_PATH
include environment
INCLUDE = -I$(Z3_PATH)/include
LIBS = -L$(Z3_PATH)/lib -lz3 -lboost_program_options
ODIR = build
SDIR = src
_SRCS = main.cpp propagator.cpp graph.cpp verifier.cpp \
	LoopUnrollPass.cpp DeclareInternalsPass.cpp LLVMModule.cpp debug.cpp \
	EdgePriorityQueue.cpp config.cpp FAnalysisPass.cpp MAnalysisPass.cpp \
	AnalysisInfo.cpp InlineClonePass.cpp scope.cpp execution.cpp
SRCS = $(addprefix $(SRCS/,$(_SRCS))
OBJS = $(addprefix $(ODIR)/, $(_SRCS:.cpp=.o))

MAIN = pxsmt

.PHONY: clean clean-logs distclean guard-Z3 mk-build

all: guard-Z3 mk-build $(MAIN)

mk-build:
	@mkdir -p ./build

guard-Z3:
	@ if [ -z ${Z3_PATH} ] || [ ! -e ${Z3_PATH} ]; then \
		echo "Environment variable Z3_PATH not set properly in './environment'"; \
		exit 1; \
    fi

$(ODIR)/%.o: $(SDIR)/%.cpp $(wildcard $(SDIR)/*.h)
	$(CXX) $(CXXFLAGS) $(LLVM_CXXFLAGS) $(INCLUDE) -c $< -o $@

$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $(MAIN) $(OBJS) $(LLVM_LDFLAGS) $(LIBS)

clean:
	@$(RM) $(ODIR)/*.o

clean-logs:
	@find logs -type f -name '*.log' -delete

distclean: clean
	@$(RM) $(MAIN)

FORCE:
