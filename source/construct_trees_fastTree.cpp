#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(int argc, char* argv[]) {
	int g, num_threads;
        char program_path[5000];

        sscanf(argv[1], "%d", &g);
        strcpy(program_path, argv[2]);
        sscanf(argv[3], "%d", &num_threads);
	
	char system_call_string[5000];
	for (int i = 1; i <= g; i++) {
		sprintf(system_call_string, "%s -nt -gtr < group_%.2d_alignment_sampled.fasta > group_%.2d_FastTree.nwk", program_path, i, i);
		system(system_call_string);
	}
	return 0;
}

