#include<stdio.h>
#include<ctype.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<algorithm>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<ctime>

using namespace std;

#define SZ 50000005	// maximum length of alignment
#define MXN 20000

int num_bootstrap;  	// number of bootstrap iterations
int num_threads;	// number of threads
int n;			// number of isolates
int C;			// number of columns, i.e. length of alignment
int sample_size; 
//char seq[MXN][SZ];
char **seq;
char full_path_raxml[10005];	// path of raxml executable
vector<char> flag;
vector<string> names_taxa;  // beginning with '>' and ending with '\n' 
char input_fasta_filename[1005];	// input alignment
char fasta_filename_with_groups[1500];	// after selection of isolates, this fasta stores alingment containning only those seq with groupnames
char buff[SZ];
vector< pair<int, vector<int> > > groups;
int represent[MXN];
int user_provided_sample_size;

bool is_huge = false;
char huge_seq[2][SZ];

int seq_distance(int i, int j) {
	int sum = 0;
	for (int k = 0; k < C; k++)
		if( seq[i][k] != seq[j][k])
			sum++;
	return sum;
}

void read_huge_sequence(char idx, int i) {
	// idx is 0 or 1, to denote where it will store the huge sequence
	// i is the isolate number
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

void initialize() {
	char system_call_string[5005];
	
	sprintf(system_call_string, "grep \">\" %s | wc > size.txt", input_fasta_filename);
	int ret_system = system(system_call_string);

	sprintf(system_call_string, "wc %s >> size.txt", input_fasta_filename);
	ret_system = system(system_call_string);

	int garb[5];
	long long tot;
	FILE* fpi = fopen("size.txt", "r");
	fscanf(fpi, "%d%d%d%d%d%lld", &n, &garb[0], &garb[1], &garb[2], &garb[3], &tot);
	fclose(fpi);
	
	if (tot > 2000000000) // if greater than 2 GB, too large for memory
	{
		is_huge = true;
		return;
	}

	seq = (char **)malloc(n * sizeof(int *));
	for (int i = 0; i < n; i++)
		seq[i] = (char *)malloc((tot/n) * sizeof(char));

}

void read_fasta() {
	
	char system_call_string[500];
	sprintf(system_call_string, "./split_fasta %s", input_fasta_filename);
	int ret_system = system(system_call_string);

	sprintf(system_call_string, "grep \">\" %s > isolate_names_%s.txt", input_fasta_filename, input_fasta_filename);
        ret_system = system(system_call_string);

	if (is_huge) {
		// process fasta file
		char isolates_name_file[500];
		sprintf(isolates_name_file, "isolate_names_%s.txt", input_fasta_filename);	
		FILE* fpi = fopen(isolates_name_file, "r");
		char isolate_name[1005];
		while(fgets(isolate_name, sizeof(isolate_name), fpi)) {
			if (isolate_name[0] == '>')
				names_taxa.push_back(isolate_name);
		}
		fclose(fpi);
		
		C = 0;
		read_huge_sequence(0, 0);		
		while (huge_seq[0][C])
			C++;

		for (int i = 0; i < n; i++)
			flag.push_back(1);

		return;  
	}

	FILE* fpi = fopen(input_fasta_filename, "r");
	names_taxa.clear();
	flag.clear();
	n = -1;
	int filled;
	while(fgets(buff, sizeof(buff), fpi)) {
		if (buff[0] == '>') {
			names_taxa.push_back((string)buff);
			if (n >= 0) {
				seq[n][filled] = 0;
				C = filled;
			}
			flag.push_back(1);
			n++;
			filled = 0;
			continue;
		}
		for (int i = 0; buff[i]; i++)
			if (isalpha(buff[i])) {
				seq[n][filled] = buff[i];
				filled++;
			}
	}
	seq[n][filled] = 0;
	n++;
	fclose(fpi);
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

char group_filename[1005];

void print_groups() {
	strcpy(group_filename, input_fasta_filename);
	int v = strlen(group_filename);
	while (v > 0 && group_filename[v] != '.')
		v--;
	group_filename[v] = 0;
	strcat(group_filename, "_groups.txt");
	FILE* fpo = fopen(group_filename, "w");
	 for (int i = 0; i < groups.size(); i++)
                if (groups[i].second.size() + 1 >= 4) {
			for (int j = 0; j < names_taxa[groups[i].first].size(); j++)
				if (names_taxa[groups[i].first][j] != '>' && names_taxa[groups[i].first][j] > 32)
					fprintf(fpo, "%c", names_taxa[groups[i].first][j]);
			for (int k = 0; k < groups[i].second.size(); k++) {
				fprintf(fpo, ",");
				int idx = groups[i].second[k];
				for (int j = 0; j < names_taxa[idx].size(); j++)
					if ( names_taxa[idx][j] != '>' &&  names_taxa[idx][j] > 32)
						fprintf(fpo, "%c", names_taxa[idx][j]);
			}
			fprintf(fpo, "\n");
		}
	fclose(fpo);
}

// will group isolates which are scale times closer than the most distant isolates
int select_isolates(int lim_rand, double scale = 50.0) {
	int i, j;
	// now we find the approximate maximum distance
	int mx_dist = 0;
	while (lim_rand--) {		
		int idx = rand() % n;
		for (i = 0; i < n; i++) {
			int v = is_huge? huge_seq_distance(i, idx) : seq_distance(i, idx);
			if (v > mx_dist)
				mx_dist = v;
		}
	}

	//double scale = 40.0;	//found by simulation	//bin_search_scale(mx_dist);
	//fprintf(stderr, "scale found %lf\n", scale);
	printf("partition of %d isolates\n", n);

	groups.push_back(make_pair(0, vector<int>()));
	for (i = 1; i < n; i++) {
		for (j = 0; j < groups.size() && flag[i]; j++) {
			int v = is_huge? huge_seq_distance(i, groups[j].first) : seq_distance(i, groups[j].first);
			if (v * scale < mx_dist) {
				flag[i] = 0;
				groups[j].second.push_back(i);
			}
		}
		if (flag[i])
			groups.push_back(make_pair(i, vector<int>()));
	}
	int cnt_grp = 0;
	for (i = 0; i < groups.size(); i++) 
		if (groups[i].second.size() + 1 < 4) {
			for (j = 0; j < groups[i].second.size(); j++)
				flag[groups[i].second[j]] = 1;
		}
		else {
			cnt_grp++;
			represent[groups[i].first] = cnt_grp;
		}
	print_groups();
	
	int high_res = 0;
	for (i = 0; i < n; i++) if (flag[i]) {
		high_res++;
	}
	// rule found from simulation, can be improved by more simulations: more optimal lengths for different sample sizes
	if (high_res <= 30)
		sample_size = 8000;
	else if (high_res <= 60)
		sample_size = 20000;
	else
		sample_size = 30000;
	if (sample_size < user_provided_sample_size)
		sample_size = user_provided_sample_size;
	return cnt_grp;
}

// here we have used reservoir sampling: equal probability of selecting columns

vector<int> reservoir() {
	int i;
	vector<int> ret(sample_size,0);
	for (i = 0; i < sample_size && i < C; i++)
		ret[i] = i;	// 1st sample_size columns are chosen initially
	for (; i < C; i++) {
		int idx = gene_idx(i+1);
		if (idx < sample_size)
			ret[idx] = i;	// if randomly picked index is less than sample_size, it is replaced
	}
	return ret;
}

void build_one_sample(int idx) {
	
	srand ( (unsigned)  time(0)  );
	vector<int> now = reservoir();

	char out_file_name[200];
	sprintf(out_file_name, "alignment_sampled_%.2d_%s", idx, input_fasta_filename);
	FILE* fpo = fopen(out_file_name, "w");

	for (int i = 0; i < n; i++) if (flag[i]) {
		if (represent[i])
			fprintf(fpo, ">group%.2d\n", represent[i]);
		else
			fprintf(fpo, "%s", names_taxa[i].c_str());
		if (is_huge) {
			read_huge_sequence(0, i);
		}
		for (int j = 0; j < now.size(); j++) {
			fprintf(fpo, "%c", is_huge? huge_seq[0][now[j]] : seq[i][now[j]]);
		}
		fprintf(fpo, "\n");
	}
	fclose(fpo);
}

void print_fasta_with_groupname() {
	strcpy(fasta_filename_with_groups, input_fasta_filename);
	int v = strlen(fasta_filename_with_groups);
        while (v > 0 && fasta_filename_with_groups[v] != '.')
                v--;
        fasta_filename_with_groups[v] = 0;
        strcat(fasta_filename_with_groups, "_with_groups_high_res.fasta");
	FILE* fpo = fopen(fasta_filename_with_groups, "w");
	for (int i = 0; i < n; i++) if (flag[i]) {
                if (represent[i])
                        fprintf(fpo, ">group%.2d\n", represent[i]);
                else
                        fprintf(fpo, "%s", names_taxa[i].c_str());
		if (is_huge) {            
                	read_huge_sequence(0, i);
		}
		for (int j = 0; j < C; j++) {
			if (j && j%70 == 0)
				fprintf(fpo, "\n");
			fprintf(fpo, "%c", is_huge? huge_seq[0][j] : seq[i][j]);
		}
		fprintf(fpo, "\n");
	}
	fclose(fpo);
}

char tree_cons_prog;
char tree_cons_prog_path[1005];

bool is_branch_length(string branch_len) {
	int count = 0;
	for (int i = 0; i < branch_len.size(); i++) {
		if (!isdigit(branch_len[i]) && branch_len[i] != '.')
			return false;
		if (branch_len[i] == '.')
			count++;
	}
	return count == 1;
}

int count_too_small_branch(string tree_newick, double scale) {
	int i;
	for (i = 0; i < tree_newick.size(); i++)
		if (tree_newick[i] == ',' || tree_newick[i] == ':' || tree_newick[i] == '(' || tree_newick[i] == ')')
			tree_newick[i] = ' ';
	
	char nows[5000];
	sprintf(nows, "%s", tree_newick.c_str());
	int bs = 0, db;
	char br_len[105];
	vector<double> branch_lengths;
	double max_len = 0;
	while (sscanf(nows + bs, "%s%n", br_len, &db) == 1) {
		bs += db;
		if (is_branch_length((string)br_len)) {
			double now_v;
			sscanf(br_len, "%lf", &now_v);
			if (now_v > max_len)
				max_len = now_v;
			branch_lengths.push_back(now_v);
		}	
	}
	
	int ret = 0;
	for (i = 0; i < branch_lengths.size(); i++) {
		if (max_len / branch_lengths[i] > scale)
			ret++;
	}
	return ret;
} 

void build_samples(int K, double scale) {
	print_fasta_with_groupname();	// whole alignment of high resolution isolates
	
	int count_prev = 0;
	int iter_limit = 3;
	char system_call_string[2005];
	while (iter_limit--) {

		for (int i = 0; i < K; i++) {
			build_one_sample(i);

			// compute maximum likelihood tree - according to user's choice
			if (tree_cons_prog == 'R') { // run RAxML
				sprintf(system_call_string,"%s -f a -x 45079 -p 35217 -# %d -m GTRGAMMA -s alignment_sampled_%.2d_%s -n tree_from_sample_%.2d_%s -T %d",
					full_path_raxml, num_bootstrap, i, input_fasta_filename, i, input_fasta_filename, num_threads);
			} else if (tree_cons_prog == 'I') { // run IQtree
				sprintf(system_call_string, "%s -s alignment_sampled_%.2d_%s -m GTR+G4 -nt %d", tree_cons_prog_path, i, input_fasta_filename, num_threads);
			} else {	// run FastTree
				sprintf(system_call_string, "%s -nt -gtr < alignment_sampled_%.2d_%s > alignment_sampled_%.2d_%s.nwk", tree_cons_prog_path, i, input_fasta_filename, i, input_fasta_filename);
			}
        		int ret_system = system(system_call_string);
		}
	
		int mx_idx = 0;
        	if (tree_cons_prog == 'R')
                	sprintf(system_call_string, "cat RAxML_bestTree.tree_from_sample_%.2d_%s > best_high_resolution_tree_%s.nwk", mx_idx, input_fasta_filename, input_fasta_filename);
        	else if (tree_cons_prog == 'I')
                	sprintf(system_call_string, "cat alignment_sampled_%.2d_%s.treefile > best_high_resolution_tree_%s.nwk", mx_idx, input_fasta_filename, input_fasta_filename);
        	else
                	sprintf(system_call_string, "cat alignment_sampled_%.2d_%s.nwk > best_high_resolution_tree_%s.nwk", mx_idx, input_fasta_filename, input_fasta_filename);
        	int ret_system = system(system_call_string);
	
		char tree_newick[5000];
		char tree_filename[500];
		sprintf(tree_filename, "best_high_resolution_tree_%s.nwk", input_fasta_filename);
		FILE* fpi_t = fopen(tree_filename, "r");
		fgets(tree_newick, sizeof(tree_newick), fpi_t);
		fclose(fpi_t);

		int count_now = count_too_small_branch((string)tree_newick, scale);
		if (sample_size >= C || count_now <= count_prev + 1)
			break;
		sample_size *= 2;
	}
	if (is_huge == false) {
		for (int j = 0; j < n; j++)
                	free(seq[j]);
        	free(seq);
	}	
	/* this part is excluded: TODO: statistical analysis later to see if the likelihood increases or RF ditance decreases.
	if (tree_cons_prog == 'R') 
		sprintf(system_call_string, "cat RAxML_bestTree.tree_from_sample_* > trees_all_sample.txt");
	else if (tree_cons_prog == 'I')
		sprintf(system_call_string, "cat alignment_sampled_*.fasta.treefile > trees_all_sample.txt");
	else
		sprintf(system_call_string, "cat alignment_sampled_*.nwk > trees_all_sample.txt");
	int ret_system = system(system_call_string);


	// computing likelihood of this ML tree w.r.t. all columns
        // so that we can find the tree with maximum value of likelihood
	sprintf(system_call_string,"%s -f N -z trees_all_sample.txt -s %s -m GTRGAMMA -n all_samples_lh -T %d",
                        full_path_raxml, fasta_filename_with_groups, num_threads);	
	ret_system = system(system_call_string);

	FILE* fpi = fopen("RAxML_info.all_samples_lh", "r");
	char nowb[5000];
	int idx, mx_idx = -1;
	double mx_lh, lh;
	while(fgets(nowb, sizeof(nowb), fpi)) {
		if(nowb[0] < '0' || nowb[0] > '9')
			continue;
		sscanf(nowb, "%d%lf", &idx, &lh);
		if (mx_idx == -1)
			mx_idx = idx, mx_lh = lh;
		else if (mx_lh < lh)
			mx_idx = idx, mx_lh = lh;
	}
	fclose(fpi);
	*/
} 	

void construct_low_resolution_trees(int num_groups) {
	char system_call_string[2005];	
	
	sprintf(system_call_string, "./construct_next_level_alignments %s %s", input_fasta_filename, group_filename);
	int ret_system = system(system_call_string);

	for (int i = 1; i <= num_groups; i++) {
		sprintf(system_call_string, "./HPC_new merged_alignment_group_%d_%s -t %d -p %c", i, input_fasta_filename, num_threads, tree_cons_prog);
		ret_system = system(system_call_string);
	}
}

// 1st argument is Fasta filename, 2nd argument is K -> K samples are created
// optional argument is list of isolates who are selected -> not implemented
int main(int argc, char* argv[]) {
	int no_arg = 2;
	if(argc < no_arg)
	{
		fprintf(stderr,"%s inputFile1\n\n", argv[0]);
		fprintf(stderr,"  inputFile1: alignment file containing the sequences for all isolates in fasta format\n");
		fprintf(stderr,"  options:\n");

		fprintf(stderr,"  -s N specify number of randomly chosen isolates to compute maximum sequence dissimilarity <= 20 (5) \n");
		//fprintf(stderr,"  -c N specify the number of columns to be sampled >= 100 (defaul values are found by simulation for different values of high resolution isolates) \n");
		fprintf(stderr,"  -b N specify the number of bootstrap iterations for trees generated by RAxML >= 20 (40) \n");
		//fprintf(stderr,"  -k N specify the number of samples and corresponding trees we want to generate >= 1 (5)\n");
		fprintf(stderr,"  -t N specify number of threads to be requested >= 2 (10)\n");
		fprintf(stderr,"  -p C choices are 'I' for IQtree or 'F' for FastTree followed by full path of the executible ('R' for RAxML)\n\n");
		exit(1);
	}
	user_provided_sample_size = 500;
	num_threads = 10;
	num_bootstrap = 40;
	tree_cons_prog = 'R';
	strcpy(full_path_raxml, "./raxmlHPC-PTHREADS-SSE3");
	strcpy(tree_cons_prog_path, full_path_raxml);
	int num_rand_iter = 1;
	int i, number;
	int K = 1;
	if(argc > no_arg)
	{
		for(i = no_arg; i < argc - 1; i += 2)
		{
			if(argv[i][0] != '-')
			{
				fprintf(stderr,"Each option must begin with a dash (-)\n");
				exit(1);
			}
			switch(argv[i][1])
			{
				case 'b':
					sscanf(argv[i+1], "%d", &number);
					if(number < 20)
					{
						fprintf(stderr, "Number of bootstrap iterations should be at least 20\n");
						exit(1);
					}
					num_bootstrap = number;
					break;
				/*case 'c':
					sscanf(argv[i+1], "%d", &number);
					if(number < 100)
					{
						fprintf(stderr, "Number of randomly selected columns should be at least 100\n");
						exit(1);
					}
					user_provided_sample_size = number;
					break;*/
				/*case 'k':
                    			sscanf(argv[i+1], "%d", &number);
                    			if(number < 1)
                    			{
                            			fprintf(stderr, "Number of samples should be at least 1\n");
                            			exit(1);
                    			}
                    			K = number;
                    		break;*/
				case 's':
                                        sscanf(argv[i+1], "%d", &number);
                                        if(number > 20 || number < 1)
                                        {
                                                fprintf(stderr, "Number of randomly chosen isolates are limited from 1 to 20 \n");
                                                exit(1);
                                        }
                                        num_rand_iter = number;
                                break;
				case 't':
					sscanf(argv[i+1],"%d",&number);
					if(number < 2)
					{
						fprintf(stderr, "Value of number of threads to be assigned to raxml should be at least 2\n");
						exit(1);
					}
					num_threads = number;
					break;
				case 'p':
					char temp_s[30];
					sscanf(argv[i+1], "%s", temp_s);
					if (temp_s[0] != 'I' && temp_s[0] != 'F' && temp_s[0] != 'R')
					{
						fprintf(stderr, "Wrong option letter for the choice of program for maximum likelihood tree construction\n");
						exit(1);
					}
					tree_cons_prog = temp_s[0];	
					//sscanf(argv[i+2], "%s", tree_cons_prog_path);
					//if (tree_cons_prog_path[0] == '!') {
						if (tree_cons_prog == 'R')
							strcpy(tree_cons_prog_path, "./raxmlHPC-PTHREADS-SSE3");
						else if (tree_cons_prog == 'I') 
							strcpy(tree_cons_prog_path, "./iqtree");
						else
							strcpy(tree_cons_prog_path, "./FastTreeMP");
					//}
					//i++;
					break; 
				default: 
					fprintf(stderr, "Wrong option letter\n");
					exit(1);
					break;
			}
		}
	}

	//strcpy(full_path_raxml, "./raxmlHPC-PTHREADS-SSE3"); 
	strcpy(input_fasta_filename, argv[1]);
	
	initialize();

	printf(" ->> initialization finished\n");

	// read the fasta file, i.e. alignment
	read_fasta();

	printf(" ->> reading sequence has completed\n");

	char system_call_string[500];
	sprintf(system_call_string, "./compute_dist_btwn_isolates_analysis %d", n);
	int ret_system = system(system_call_string);
	
	double scale;
	FILE* fpi = fopen("ratio_output.txt", "r");
	fscanf(fpi, "%lf", &scale);
	fclose(fpi);
	
	printf("====>> scale found %lf\n", scale);
	// check to do : multiple levels...
	
	if (scale < 25) {
		scale = 30;	// default value of alpha
	}	
	else if (scale > 100) {
		scale = 100.0;
	}
	// select candidate isolates by marking the flag array
	int num_groups = select_isolates(num_rand_iter, scale);

	printf(" ->> selecting isolates has finished\n");
	

	// build K samples and trees	

	build_samples(K, scale);

	printf(" ->> building alignment by sampling has finished\n");

	// build low resolution alignments and trees
	construct_low_resolution_trees(num_groups);
	return 0;
}
