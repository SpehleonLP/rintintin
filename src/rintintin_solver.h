#ifndef RINTINTIN_ERROR_FUNCTIONS_H
#define RINTINTIN_ERROR_FUNCTIONS_H
#include "../include/rintintin.h"
#include "rintintin_scratch.h"
	
rintintin_error_code rintintin_compute_child_table(struct rintintin_scratch_space * scratch, int32_t const* parents, uint32_t * counters);
rintintin_error_code rintintin_solve(rintintin_command * cmd, struct rintintin_metrics * dst, void * tmp_begin, void * tmp_end);

#endif // RINTINTIN_ERROR_FUNCTIONS_H
