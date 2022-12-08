open_project FLAT
set_top LINEAR
add_files flat.hpp
add_files flat.cpp
add_files -tb flat_tb.cpp

open_solution "solution1" -flow_target vivado
set_part {xcvu11p-flga2577-1-e}
create_clock -period 10 -name default
csim_design -clean
csynth_design
cosim_design
