#include<stdio.h>
#include<string.h>
#include<stdlib.h>


#define SZ 50000005
char huge_seq[2][SZ];

int C;

void read_huge_sequence(char idx, int i) {
        char split_filename[500];
        sprintf(split_filename, "seq_after_split_%d.txt", i);
        FILE* fpi = fopen(split_filename, "r");
        fgets(huge_seq[idx], sizeof(huge_seq[idx]), fpi);
        fclose(fpi);
}

int huge_seq_distance(int i, int j) {

        read_huge_sequence(0, i);
        read_huge_sequence(1, j);

        int sum = 0;
        for (int k = 0; k < C && huge_seq[0][k] && huge_seq[1][k]; k++)
                if (huge_seq[0][k] != huge_seq[1][k])
                        sum++;
        return sum;
}

int main(int argc, char * argv[]) {
	
	if (argc != 3)
		exit(1);
	int u, v;
	sscanf(argv[1], "%d", &u);
	sscanf(argv[2], "%d", &v);

	read_huge_sequence(0, 0);
	while(huge_seq[0][C])
		C++;

	printf("%d\n", C);

	printf("%d\n", huge_seq_distance(u, v));

	return 0;
}
