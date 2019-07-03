#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>

char buf[50000005];

int main(int argc, char * argv[]) {
	if (argc != 2) {
		fprintf(stderr, "Too few arguments\n");
		exit(1);
	}
	
	FILE* fpi = fopen(argv[1], "r");
	FILE* fpo;
	char filename[500];
	int cnt = -1;
	while(fgets(buf, sizeof(buf), fpi)) {
		if (buf[0] == '>') {
			if (cnt != -1)
				fclose(fpo);
			cnt++;
			sprintf(filename, "seq_after_split_%d.txt", cnt);
			fpo = fopen(filename, "w");
			continue;
		}

		for (int i = 0; buf[i]; i++)
			if (isalpha(buf[i]))
				fprintf(fpo, "%c", toupper(buf[i]));

	}	
	fclose(fpi);
	fclose(fpo);
	return 0;
}

