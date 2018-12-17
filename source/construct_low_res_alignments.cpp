#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<algorithm>
#include<vector>
#include<string>

using namespace std;

#define MXN 200
#define SZ 5000000	// 5*10^6

char buf[1000005];
int n;                  // number of isolates
int C;                  // number of columns, i.e. length of alignment
char seq[MXN][SZ];
vector<string> names_taxa;  // beginning with '>' and ending with '\n'
int sample_size;	// size of the reservoir

void read_fasta(int g) {
	char fasta_filename[150];
	sprintf(fasta_filename, "group_%.2d_alignment.fasta", g);
        FILE* fpi = fopen(fasta_filename, "r");
        names_taxa.clear();
        n = -1;
        int filled;
        while(fgets(buf, sizeof(buf), fpi)) {
                if (buf[0] == '>') {
                        names_taxa.push_back((string)buf);
                        if (n >= 0) {
                                seq[n][filled] = 0;
                                C = filled;
                        }
                        n++;
                        filled = 0;
                        continue;
                }
                for (int i = 0; buf[i]; i++)
                        if (isalpha(buf[i])) {
                                seq[n][filled] = buf[i];
                                filled++;
                        }
        }
        seq[n][filled] = 0;
        n++;
        fclose(fpi);
}

bool is_same(int idx) {
	for (int i = 1; i < n; i++) 
		if (seq[i][idx] != seq[i-1][idx])
			return false;
	return true;
}

int gene_idx(int v) {
        int ret = 0;
        int len = 0;
        int temp = v;
        while (temp > 0)
                temp /= 10, len++;
        while (len)
                ret = ret * 10 + rand() % 10,
                len--;
        return (ret % v);
}

vector<int> reservoir(int v) {
        int i;
        vector<int> ret(sample_size,0);
        for (i = 0; i < sample_size; i++)
                ret[i] = i;     // 1st sample_size columns are chosen initially
        for (; i < v; i++) {
                int idx = gene_idx(i+1);
                if (idx < sample_size)
                        ret[idx] = i;   // if randomly picked index is less than sample_size, it is replaced
        }
        return ret;
}

void sample_columns(int g) {
	read_fasta(g);
	
	vector<int> not_same;
	for (int i = 0; i < C; i++)
		if (is_same(i) == false)
			not_same.push_back(i);

	char out_filename[150];
	sprintf(out_filename, "group_%.2d_alignment_sampled.fasta", g);
	FILE* fpo = fopen(out_filename, "w");

	if (not_same.size() > sample_size) {
		vector<int> now = reservoir(not_same.size());
		for (int i = 0; i < n; i++) {
                        fprintf(fpo, "%s", names_taxa[i].c_str());
                        for (int j = 0; j < now.size(); j++) {
                                if (j && j % 70 == 0)
                                        fprintf(fpo, "\n");
                                fprintf(fpo, "%c", seq[i][not_same[now[j]]]);
                        }       
                        fprintf(fpo, "\n");
                }
	} else {
		for (int i = 0; i < n; i++) {
			fprintf(fpo, "%s", names_taxa[i].c_str());
			for (int j = 0; j < not_same.size(); j++) {
				if (j && j % 70 == 0)
					fprintf(fpo, "\n");
				fprintf(fpo, "%c", seq[i][not_same[j]]);
			}
			fprintf(fpo, "\n");	
		}
	}
	
	fclose(fpo);
}

int main(int argc, char * argv[]) {
	if (argc < 3) {
		fprintf(stderr, "%s input_fasta_filename group_filename_for_low_resolution\n\n", argv[0]);
		exit(1);
	}	
	char system_call_string[15236];
	char first_isolate[105], rest_isolates[12345];
	FILE* fpi = fopen(argv[2], "r");
	int g = 0;
	sample_size = 10000;
	while (fgets(buf, sizeof(buf), fpi)) {
		for(int i = 0; buf[i]; i++) 
			if (buf[i] == ',')	
			{
				buf[i] = ' ';
				break;
			}
		sscanf(buf, "%s %s", first_isolate, rest_isolates);
		g++;
		sprintf(system_call_string,"./partial_align_seq %s %s group_%.2d_alignment.fasta", argv[1], rest_isolates, g);
		int ret_system = system(system_call_string);
		
		sample_columns(g);
	}	
	fclose(fpi);
	return 0;
}
