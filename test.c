#include "stdio.h"

// test queue does not support --exclusive. It IS needed for proper runtime measurements, mind you. 
// #!/bin/csh loads a C shell, so we don't need source /etc/profile.d/modules.sh
// sbatch: error: Batch job submission failed: Job violates accounting/QOS policy (job submit limit, user's size, and/or time limits)

// I need to fix this error. Let's try a different queue.

void main() {

	printf("Hello HPC@Durham from Hamilton\n");

}
