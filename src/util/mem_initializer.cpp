#include "util/debug.h"
#include "util/rlimit.h"
#include "util/scoped_timer.h"
void mem_initialize() {
initialize_rlimit();
scoped_timer::initialize();
}
void mem_finalize() {
finalize_debug();
finalize_rlimit();
}
