#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<algorithm>
#include<string>
#include<vector>
#include<map>
using namespace std;

char input_fasta_filename[1000];

map<string, int> isolate_names;
vector<string> isolate_idx_to_names;

void read_isolate_names() {
	char filename[500];
	sprintf(filename, "isolate_names_%s.txt", input_fasta_filename);
	FILE* fpi = fopen(filename, "r");
	char buf[500];
	int count = 0;
	isolate_names.clear();
	isolate_idx_to_names.clear();	
	while(fgets(buf, sizeof(buf), fpi)) {
		if (buf[0] != '>')
			continue;
		for(int i = 0; buf[i]; i++)
			if (buf[i] == '\n') {
				buf[i] = 0;
				break;
			}
		isolate_names[(string)(buf + 1)] = count;
		isolate_idx_to_names.push_back((string)(buf + 1));
		count++;
	}
	fclose(fpi);
}

void build_alignment(int idx, char* all_isolates) {
	int i;
	for (i = 0; all_isolates[i]; i++)
		if (all_isolates[i] == ',')
			all_isolates[i] = ' ';
	
	char now_isolate[500];
	vector<int> isolate_indices;

	int bs = 0, db;
	while (sscanf(all_isolates + bs, "%s%n", now_isolate, &db) == 1) {
		bs += db;
		isolate_indices.push_back(isolate_names[(string)now_isolate]);
	}
	
	FILE* fpi[10000];
	FILE* fpo[10000];
	for (i = 0; i < isolate_indices.size(); i++) {
		char read_seq_filename[500], write_seq_filename[500];
		sprintf(read_seq_filename, "seq_after_split_%d.txt", isolate_indices[i]);
		sprintf(write_seq_filename, "seq_to_merge_%d.txt", i);
		fpi[i] = fopen(read_seq_filename, "r");
		fpo[i] = fopen(write_seq_filename, "w");	
	}

	char c_now[10000];
	while (1) {
		c_now[0] = fgetc(fpi[0]);
		if (c_now[0] == EOF)
			break;
		int not_same = 0;
		for (i = 1; i < isolate_indices.size(); i++) {
			c_now[i] = fgetc(fpi[i]);

			if(c_now[i] != c_now[0]) {
				not_same = 1;
			} 
		}
		if (not_same) {
			for (i = 0; i < isolate_indices.size(); i++) {
				fprintf(fpo[i], "%c", c_now[i]);
			}
		}	
	}

	for (i = 0; i < isolate_indices.size(); i++) {
		fclose(fpi[i]);
	 	fclose(fpo[i]);
	}	

	for(i = 0; i < isolate_indices.size(); i++) {
		char write_seq_filename[500];
                sprintf(write_seq_filename, "seq_to_merge_%d.txt", i);
		fpi[i] = fopen(write_seq_filename, "r");
	}	

	char filename[500];
	sprintf(filename, "merged_alignment_group_%d_%s", idx, input_fasta_filename);
		
	FILE* fpo_m = fopen(filename, "w");

	for (i = 0; i < isolate_indices.size(); i++) {
		fprintf(fpo_m, ">%s\n", isolate_idx_to_names[isolate_indices[i]].c_str());
		int count = 0;
		while(1) {
			char now_c = fgetc(fpi[i]);
			if (now_c == EOF)	
				break;
			fprintf(fpo_m, "%c", now_c);
			if (count % 70 == 0 && count)
				fprintf(fpo_m, "\n");
			count++;
		}
		fclose(fpi[i]);
		fprintf(fpo_m, "\n");
	}	

	fclose(fpo_m);
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
                fprintf(stderr, "%s input_fasta_filename group_filename_for_low_resolution\n\n", argv[0]);
                exit(1);
        }	
	
	char all_isolates[10005], buf[20005];
	strcpy(input_fasta_filename, argv[1]);

	read_isolate_names();

        FILE* fpi = fopen(argv[2], "r");
        int g = 0;
        while (fgets(buf, sizeof(buf), fpi)) {

                sscanf(buf, "%s", all_isolates);
                g++;

                build_alignment(g, all_isolates);
        }
        fclose(fpi);	

	return 0;
}
