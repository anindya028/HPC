#include<stdio.h>
#include<string.h>
#include<algorithm>
#include<string>
#include<set>

using namespace std;

char filename[505];
char contig_name[2000][205];

set<string> set_contigs;
char buf[3000005];

int main(int argc,char *argv[])
{
	int no_arg = 4;
	if(argc < no_arg)
	{
		fprintf(stderr,"Usage: %s fasta_filename contig_names output_filename\n\n",argv[0]);
		fprintf(stderr,"  fasta_filename: name of input file which contains all the contigs in fasta format\n");
		fprintf(stderr,"  contig_names: comma separated list of isolate names\n");
		fprintf(stderr,"  ouput_filename: output filename of fasta file containing only alignment of selected isolates\n");
	}
	strcpy(filename,argv[1]);

	strcpy(buf,argv[2]);
	int i;
	for(i = 0; buf[i]; i++)
		if(buf[i] == ',')
			buf[i] = ' ';
	int bs = 0, db;
	i = 0;
	set_contigs.clear();
	while(sscanf(buf+bs,"%s%n",contig_name[i],&db)==1)
	{
		bs += db;
		set_contigs.insert(contig_name[i]);
		i++;
	}
//	freopen(filename,"r",stdin);
	FILE* fpi = fopen(filename,"r");
	FILE* fpo = fopen(argv[3],"w");
	int on = 0;
	while(fgets(buf,sizeof(buf),fpi))
	{
		if(buf[0] == '>')
		{
			char name[505];
			buf[0] = ' ';
			sscanf(buf,"%s",name);
			if(set_contigs.find((string)name) != set_contigs.end())
				on = 1;
			else
				on = 0;
			buf[0] = '>';
		}
		if(on)
			fprintf(fpo,"%s",buf);
	}
	return 0;
}
