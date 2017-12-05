# finding_stages
Finding execution stages by memory access trace analysis

This repository contains software to create and analyse a memory access trace to find execution stages in any program

## Components
1. Valgrind records the memory trace to a file
2. find_stages analyses that memory trace

## Modifications to Valgrind
The only major modifications are in the "cachegrind" module.
"cg_sim.c" contains additional methods to write the information to a file.
The "cachesim_XX_doref" methods take additional parameters to pass the code location along.
In "cg_main.c" these additional parameters are passed to the simulation.

The only other change is in coregrind/m_libcproc.c where the new method "read_nanosecond_timer" returns a timestamp in ns.

## Building and running
1. Valgrind

   The output location of the access trace is hardcoded in "cg_sim.c", it is set to "/home/schindle/cg_X.csv".
   Change it before compiling.

   To build valgrind, simply execute "valgrind-3.13.0/mybuild".
   This will compile the software and install it to "/tmp/valgrind/".

   To record a trace of any application execute:

   "sudo /tmp/valgrind/bin/valgrind --tool=cachegrind --branch-sim=yes --cache-sim=yes --I1=32768,8,64 --D1=32768,8,64 --LL=262144,8,64 <APPLICATION>"

   --I1, D1, LL specify the cache size (32KiB L1 and 256KiB L2 in this case), associativity and line size.

2. find_stages

   To build the software, a working rust toolchain is required.

   Simply move into the directory and execute "cargo build --release".

   The executable is located under ./target/release/find_stages

   It will look for the following parameters:
   - "<"file">": This is the input file. MUST be the first argument
     The file is in the csv format, with a semicolon as seperator between values
     Each line looks like this:
     "<"Operation">";"<"Miss">";"<"File">";"<"Function">";"<"Line">";"<"Address">";"<"Size">";"<"Time">"
   - filter="<"file">": Filter all functions from the input that occur in this trace. File format must be the same as the main input file.
   - long: Simulate all the lowest overlaps in each segment. This will take a long time (see code comments for more)
   - brute-force: Brute force the analysis and skip the heuristics (This will take a long time).
