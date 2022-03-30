#include <ctype.h>
#include <math.h>
#include <cmath> // the library used to calculate the CDF function
#include <algorithm>
#include <Python.h> // call the python functions
#include <deque>
#include <time.h>
#include <iomanip>
#include <limits.h>
#include <numeric> // include the function to calculate the sum of a vector
#include "Markov.cc"

using namespace std;

#define TRIVIAL_LEN 6
#define DEL_LEN 6
#define OVERLAP 8
#define PRECISE 4
#define MULTIPLE 100
#define MIN(a, b) ((a)<(b)?(a):(b))
#define ALPHA 0.01 // the p-value cutoff (without correction) to select significant k-mers
#define BETA 0.05 // the p-value cutoff (without correction) to select sub-significant k-mers

typedef struct str_int_int
{
	string header;
	char strand;
	int begin;
	int end;
} st;

typedef vector<int> int_vector;
typedef vector<float> float_vector;
typedef vector<char> char_vector;

typedef struct kmer_occurr
{
	vector<string> kmer;
	int_vector occurr;
	float_vector cover1;
	float_vector cover2;
	float_vector score;
	vector<int_vector> site;
} kmer_set;

typedef struct kmer_lmer
{
	string kmer;
	int occurr;
	float cover1;
	float cover2;
	float score;
	vector<int> site;
	vector<float_vector> left;
	vector<float_vector> right;
} kmer_all;

typedef struct psm_occurr
{
	vector<int> set;
	vector<float_vector> mat;
	int occurr;
	string cons;
	float score;
	float cover1;
	float cover2;
} psm;

typedef struct left_psm_right
{
	vector<int> set;
	vector<float_vector> mat;
	int occurr;
	vector<float_vector> left;
	vector<float_vector> right;
	vector<int> site;
	int begin;
	int end;
	float cover1;
	float cover2;
	float score;
	double pvalue; // chi square p-value
} pre_mtf;

typedef struct set_mat_alength_nsites_cons
{
	vector<int> set;
	vector<float_vector> mat;
	int alength;
	int nsites;
	string cons;
	string rev_cons;
	string deg;
	string rev_deg;
	float score;
	vector<st> site;
	float cover1;
	float cover2;
	double pvalue; // chi square p-value
} mtf;

typedef struct mat_vec_mat
{
	vector<float_vector> pfm;
	vector<float> ic;
	vector<float_vector> pssm;
} spic;

int f_in_id = -1;
int f_bg_id = -1;
int f_out_id = -1;
int top_num;

int kmer_length = 8;
int lmer_length = 6;
int num_deg = 1;
int hd_thr = 1;
int redundant_thr = 2;
int num_mtf = -1;
int str_flag = 2;
int help_flag = 0;
int num_iter = 100;
float sw_thr = 1.8;

float sum_exp = 0;
float sum_bg = 0;

float thr1 = 8;
float thr2 = 1.96;
float thr3 = 4.5;

char alphabet[4] = {'A', 'C', 'G', 'T'};
map<char, int> r_alphabet;
float nt_freq[4];
char wild_card[14], reverse_wild[14], begin_pkg[64], end_pkg[64];
time_t begin_t, end_t;

/*************************************************************************************************
 *
 * FUnction free_pt used to free the saving space of pointer
 *
 ************************************************************************************************/

template <typename T>

void free_pt(T* tmp_pt)
{
	if(tmp_pt != NULL)
	{
		delete [] tmp_pt;
		tmp_pt = NULL;
	}
}

/*************************************************************************************************
 *
 * Function reverse_sort used to reverse the sort order
 *
 ************************************************************************************************/

void reverse_sort(int id[], int length)
{
	int *temp_id;
	temp_id = new int [length];
	int i;

	for (i = 0; i < length; i ++)
	{
		temp_id[i] = id[length - 1 - i];
	}
	for (i = 0; i < length; i ++)
	{
		id[i] = temp_id[i];
	}

	free_pt(temp_id);
}

/*************************************************************************************************
 *
 * Function usage used to print the usage of this program
 *
 ************************************************************************************************/

void usage()
{
	cout << "\n=================================================================================================================================================================\n";
	cout << "Usage: ./ProSampler [options]\n";
	cout << "\nParameters:\n";
	cout << "-i: Name of the input file in FASTA format\n";
	cout << "-b: Name of the background file in FASTA format or order of the Markov " 
		<< "model to generate background sequences (default: 3; 3rd order Markov model)\n";
	cout << "-d: The cutoff of Hamming Distance between any two k-mers in a PWM (default: 1)\n";
	cout << "-o: Prefix of the names of output files\n";
	cout << "-m: Number of motifs to be output (default: All)\n";
	cout << "-f: Number of cycles of Gibbs Sampling to find  each preliminary motif (default: 100)\n";
	cout << "-k: Length of preliminary motifs (default: 8)\n";
	cout << "-l: Length of the flanking l-mers (default: 6)\n";
	cout << "-c: Cutoff of Hamming distance to merge similar k-mers (default: 1)\n";
	cout << "-r: Cutoff of Hamming distance to delete redundant motifs basedn on consensus (default: 1)\n";
	cout << "-p: Number(1 or 2) of strands to be considered(default: 2)\n";
	cout << "-t: Cutoff of z-value to choose significant k-mers (default: 8.00)\n";
	cout << "-w: Cutoff of z-value to choose sub-significant k-mers (default: 4.50)\n";
	cout << "-z: Cutoff of z-value to extend preliniary motifs(default: 1.96)\n";
	cout << "-s: Cutoff of SW score to construct graph (default: 1.80)\n";
        cout << "-h: Print this message (default: 0)\n";
	cout << "=================================================================================================================================================================\n";
}

/*************************************************************************************************
 *
 * Function parse_opt used to parse the input parameters to this program
 *
 ************************************************************************************************/

 void parse_opt(int argc, char** argv)
 {
	 int i;

	 if (argc < 3)
	 {
		 cout << "Error: The number of input parameters is wrong. Please check the input parameters again.\n";
		 usage();
		 exit(0);
	 }

	 for (i = 1; i < argc - 1; i ++)
	 {
		 if (strcmp(argv[i], "-i") == 0) {f_in_id = i + 1; i ++;}
		 else if (strcmp(argv[i], "-b") == 0) {f_bg_id = i + 1; i ++;}
		 else if (strcmp(argv[i], "-o") == 0) {f_out_id = i + 1; i ++;}
		 else if (strcmp(argv[i], "-d") == 0) {num_deg = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-k") == 0) {kmer_length = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-l") == 0) {lmer_length = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-m") == 0) {num_mtf = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-c") == 0) {hd_thr = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-r") == 0) {redundant_thr = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-p") == 0) {str_flag = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-f") == 0) {num_iter = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-t") == 0) {thr1 = atof(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-w") == 0) {thr3 = atof(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-z") == 0) {thr2 = atof(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-h") == 0) {help_flag = atoi(argv[i + 1]); i ++;}
		 else if (strcmp(argv[i], "-s") == 0) {sw_thr = atof(argv[i + 1]); i ++;}
	 }

	 cout << "Finished parsing input parameters to the program.\n";
 }

/*************************************************************************************************
 *
 * Function de_lower used to substitue lowercase letters with uppercase letters
 *
 ************************************************************************************************/

void de_lower(sequence& seq)
{
	int i, j;

	for (i = 0; i < seq.strand.size(); i ++)
	{
		for (j = 0; j < seq.strand[i].length(); j ++)
		{
			if (islower(seq.strand[i][j]))
			{
				seq.strand[i][j] = toupper(seq.strand[i][j]);
			}
		}
	}
}

/*************************************************************************************************
 *
 * Function: de_other
 * Used to get rid of the unknown letters
 *
 ************************************************************************************************/

void de_other(sequence& seq)
{
	for (int i = 0; i < seq.strand.size(); i ++)
	{
		for (int j = 0; j < seq.strand[i].length(); j ++)
		{
			if (seq.strand[i][j] != 'A' && seq.strand[i][j] != 'C' 
					&& seq.strand[i][j] != 'G' && 
					seq.strand[i][j] != 'T')
			{
				seq.strand[i][j] = 'N';
			}
		}
	}

	cout << "Finished getting rid of the unknown letters in FASTA file.\n";
}

/*************************************************************************************************
 *
 * Function: copy_str
 * Used to concatenate strings
 *
 ************************************************************************************************/

void copy_str(string& small, string& big)
{
	for (int i = 0; i < small.length(); i ++)
	{
		if (small[i] != ' ')
		{
			big += small[i];
		}
	}
}

/*************************************************************************************************
 *
 * Function seq_reverse used to compute the reverse strand of the input sequences
 *
 ************************************************************************************************/

void seq_reverse(sequence& seq, vector<string>& strand)
{
	int i, j;
	string str;
	strand.clear();

	for (i = 0;i < seq.strand.size(); i ++)
	{
		str.assign(seq.strand[i].rbegin(), seq.strand[i].rend());
		for (j = 0; j < str.length(); j ++)
		{
			switch(str[j])
			{
				case 'A':
				str[j]='T';
				break;
				case 'C':
				str[j]='G';
				break;
				case 'G':
				str[j]='C';
				break;
				case 'T':
				str[j]='A';
				break;
				default:
				break;
			}
		}
		strand.push_back(str);
	}
}

/*************************************************************************************************
 *
 * Function seq_generate to generate the final version of sequences
 *
 ************************************************************************************************/

void seq_generate(sequence& seq_pos, sequence& seq_final, 
		int str_flag)
{
	if (str_flag == 1)
	{
		seq_final = seq_pos;
	}
	else if (str_flag == 2)
	{
		vector<string> strand;
		seq_reverse(seq_pos, strand);

		if (seq_pos.strand.size() != strand.size())
		{
//			cout << seq_pos.strand.size() << "\t" << 
//				strand.size() << "\n";

			cerr << "Error: The number of positive strands of " 
			<< "the sequences is not equal to the number of negative strands.\n";
			exit(1);
		}
		seq_final.name = seq_pos.name;
		seq_final.strand.clear();

		for (int i = 0; i < seq_final.name.size(); i ++)
		{
			seq_final.strand.push_back(seq_pos.strand[i]);
			seq_final.strand.push_back(strand[i]);
		}

		if (seq_final.name.size()*2 != seq_final.strand.size())
		{
			cerr << "Error: The number of sequence names " 
			<< "mismatches that of sequence strands.\n";
			exit(1);
		}
	}
	else
	{
		cerr << "Error: The parameter str_flag should be chosen from 1 or 2.\n" 
		<< "Please check the parameters again.\n";
		exit(1);
	}
}

/*************************************************************************************************
 *
 * Function load_data used to read data into this program
 *
 ************************************************************************************************/

void load_data(char* file_name, sequence& seq)
{
	seq.name.clear();
	seq.strand.clear();

	ifstream file_op(file_name);
	if (!file_op)
	{
		cerr << "Error: Can't open file " << file_name << "!\n";
		exit(1);
	}

	string line;
	string merge_ln;

	while (!file_op.eof())
	{
		getline(file_op, line);
			
		if (line[0] == '>')
		{
			if(!merge_ln.empty())
			{
				seq.strand.push_back(merge_ln);
			}

			seq.name.push_back(line.substr(1));
			merge_ln.clear();
		}
		else if (!line.empty())
		{
			copy_str(line, merge_ln);
			// concatenate one string at the end of another string
		}
		else
		{
			continue;
		}
	}

	if(!merge_ln.empty())
	{
		seq.strand.push_back(merge_ln);
	}

	if(seq.name.size() != seq.strand.size())
	{
		cerr << "Error: The number of sequence name lines is not " 
		<< "equal to the number of sequence lines!\n";
		exit(1);
	}

	cout << "Finished loading FASTA file " << file_name 
	<< " into sequence variable.\n";
	cout << "There're altogether " << seq.name.size() 
	<< " FASTA sequences inputted into this program.\n";
}

/*************************************************************************************************
 *
 * Function free_seq used to free the space of sequences
 *
 ************************************************************************************************/

void free_seq(sequence& seq)
{
	for (int i = 0; i < seq.name.size(); i ++)
	{
		string().swap(seq.name[i]);
	}

	vector<string>().swap(seq.name);

	for (int i = 0; i < seq.strand.size(); i ++)
	{
		string().swap(seq.strand[i]);
	}

	vector<string>().swap(seq.strand);
}

/*************************************************************************************************
 *
 * Function nt_stat used to calculate the frequencies of nucleotides
 *
 ************************************************************************************************/

void nt_stat(sequence& seq_nondeg, float* nt_freq)
{
	int i, j;
	float sum = 0;

	for (i = 0; i < 4; i ++)
	{
		nt_freq[i] = 0;
	}

	for (i = 0; i < seq_nondeg.strand.size(); i ++)
	{
		for (j = 0; j < seq_nondeg.strand[i].length(); j ++)
		{
			switch (seq_nondeg.strand[i][j])
			{
				case 'A':
				nt_freq[0]++;
				break;
				case 'C':
				nt_freq[1]++;
				break;
				case 'G':
				nt_freq[2]++;
				break;
				case 'T':
				nt_freq[3]++;
				break;
				default:
				break;
			}
		}
	}

	for (i = 0; i < 4; i ++)
	{
		sum += nt_freq[i];
	}

	for (i = 0; i < 4; i ++)
	{
		nt_freq[i] /= sum;
	}

	cout << "Finished counting the nucleotide frequencies of the input FASTA file.\n";
}

/*************************************************************************************************
 *
 * Function kmer_count used to count the occurrence numbers of k-mers
 *
 ************************************************************************************************/

void kmer_count(map<string, int_vector>& kmer_site, 
		sequence& seq_final, int& lmer_length)
{
	int i, j;

	kmer_site.clear();

	if (str_flag == 1)
	{
		for (i = 0; i < seq_final.strand.size(); i ++)
		{
			if (seq_final.strand[i].length() < 
					kmer_length+2 * lmer_length)
			{
				continue;
			}

			for (j = lmer_length; j <= seq_final.strand[i].
					length() - kmer_length-lmer_length; j ++)
			{
				if (seq_final.strand[i].substr(j, kmer_length).
						find('N') != seq_final
						.strand[i].substr(j, kmer_length).npos)
				{
					continue;
				}

				kmer_site[seq_final.strand[i].substr(j, kmer_length)]
					.push_back(i);
				kmer_site[seq_final.strand[i].substr(j, kmer_length)]
					.push_back(j);
			}
		}
	}
	else
	{
		for (i = 0; i < seq_final.strand.size(); i += 2)
		{
			if (seq_final.strand[i].length() < kmer_length
					+ 2 * lmer_length)
			{
				continue;
			}

			for (j = lmer_length; j <= seq_final.strand[i].
					length() - kmer_length - lmer_length; j ++)
			{
				if (seq_final.strand[i].substr(j, kmer_length)
						.find('N') != seq_final.strand[i]
						.substr(j, kmer_length).npos)
				{
					continue;
				}

				kmer_site[seq_final.strand[i]
					.substr(j, kmer_length)].push_back(i);
				kmer_site[seq_final.strand[i]
					.substr(j, kmer_length)].push_back(j);
			}
		}
	}

	cout << "Finished counting " << kmer_length << "-mer occurrence numbers and sites.\n";
}

/*************************************************************************************************
 *
 * Function if_same used to get rid of the k-mers with low complexity (length: 5)
 *
 ************************************************************************************************/

int if_same(string str)
{
	int i, flag;
	flag=0;

	for (i = 1; i < str.length(); i ++)
	{
		if (str[0] != str[i])
		{
			flag = 1;
			break;
		}
	}

	return flag;
}

/*************************************************************************************************
 *
 * Function check_complex used to check complexity of a string
 *
 ************************************************************************************************/

int check_complex(string str)
{
	int flag = 0;

	map<char, int> map_str;

	for (int i = 0; i < str.length(); i ++)
	{
		map_str[str[i]] ++;
	}

	if (map_str.size() >= 2)
	{
		flag = 1;
	}

	return flag;
}

/*************************************************************************************************
 *
 * Function de_trivial used to get rid of the k-mers with low complexity (unfixed length)
 *
 ************************************************************************************************/

int de_trivial(string str)
{
	int flag = 1;

	for (int i = 0; i <= str.length() - TRIVIAL_LEN; i ++)
	{
		if (!if_same(str.substr(i, TRIVIAL_LEN)))
		{
			flag = 0;
			break;
		}
	}

	if(flag == 1)
	{
		flag = check_complex(str);

	}

	return flag;
}

/*************************************************************************************************
 *
 * Function last_trivial used to get rid of the k-mers with low complexity (unfixed length)
 *
 ************************************************************************************************/

int last_trivial(string str)
{
	int flag = 1;

	for (int i = 0; i <= str.length() - DEL_LEN; i ++)
	{
		if (!if_same(str.substr(i, DEL_LEN)))
		{
			flag = 0;
			break;
		}
	}

	if (flag == 1)
	{
		flag = check_complex(str);

	}

	return flag;
}

/*************************************************************************************************
 *
 * Function ztest used to perform two proportion z-test
 *
 ************************************************************************************************/

float ztest(float& cover1, float& cover2, float& sum1, 
		float& sum2)
{
	float p1 = cover1/sum1;
	float p2 = cover2/sum2;
	float p = (cover1 + cover2) / (sum1 + sum2);
	float up = p1 - p2;
	float down = sqrt(p * (1 - p) * (1 / sum1 + 1 / sum2));
	float res = up / down;
	return res;
}

/*************************************************************************************************
 *
 * Function init_array used to initiate the array
 *
 ************************************************************************************************/

void init_array(int* flag_ar, int len, int start)
{
	for(int i = 0; i < len; i ++)
	{
		flag_ar[i] = start;
	}
}

/*************************************************************************************************
 *
 * Function site2cover used to transform site information into coverage information
 *
 ************************************************************************************************/

float site2cover(vector<int>& site, int seq_num)
{
	int* seq_ar;
	seq_ar = new int [seq_num];
	init_array(seq_ar, seq_num, 0);
	float sum = 0;

	if (str_flag == 2)
	{
		for (int i = 0; i < site.size(); i += 2)
		{
			seq_ar[site[i] / 2] = 1;
		}
	}
	else
	{
		for (int i = 0; i < site.size(); i += 2)
		{
			seq_ar[site[i]] = 1;
		}
	}

	for (int i = 0; i < seq_num; i ++)
	{
		sum += seq_ar[i];
	}

	free_pt(seq_ar);

	return sum;
}

/*************************************************************************************************
 *
 * Function normalCDF used to calculate the CDF 
 *
 ************************************************************************************************/
 
 double normalCDF(double x)
 // x : the number of k-mers
 {
     PyObject* myModuleString = PyString_FromString((char*)"PvalueToZscore");
     PyObject* myModule = PyImport_Import(myModuleString);
     PyObject* myFunction = PyObject_GetAttrString(myModule,(char*)"pvalue_to_zscore");
     PyObject* args = PyTuple_Pack(1, PyFloat_FromDouble(x));
     PyObject* myResult = PyObject_CallObject(myFunction, args);
     return(double PyFloat_AsDouble(myResult));
 }

/*************************************************************************************************
 *
 * Function choose_kmer used to calculate significant k-mers
 *
 ************************************************************************************************/

void choose_kmer(map<string, int_vector>& map_exp, 
		map<string, int_vector>& map_bg, kmer_set& major_set, 
		int seq_num, map<string, int>& fa_flag, kmer_set& minor_set)
{
	major_set.kmer.clear();
	major_set.cover1.clear();
	major_set.cover2.clear();
	major_set.site.clear();
	major_set.score.clear();

	float tmp_cover1, tmp_cover2;
	float z_sc;
	sum_exp = 0;
	sum_bg = 0;
	thr1 = normalCDF(ALPHA / double(map_exp.size())); // find the z-score cutoff to select significant k-mers
	thr3 = normalCDF(BETA / double(map_exp.size())); // find the z-score cutoff to select significant k-mers

	for (map<string, int_vector>::iterator i = 
			map_exp.begin(); i != map_exp.end(); i ++)
	{
		if (fa_flag[i -> first] == 0)
		{
			continue;
		}

		sum_exp += (i -> second).size() / 2;
	}

	for (map<string, int_vector>::iterator i = map_bg.begin(); 
			i != map_bg.end(); i ++)
	{
		sum_bg += (i -> second).size() / 2;
	}

	for (map<string, int_vector>::iterator i = map_exp.begin(); 
			i != map_exp.end(); i ++)
	{
		if (fa_flag[i -> first] == 0)
		{
			continue;
		}

		tmp_cover1 = (i -> second).size() / 2;
		tmp_cover2 = map_bg[i -> first].size() / 2;
		z_sc = ztest(tmp_cover1, tmp_cover2, sum_exp, sum_bg); // two proportion z-test

		if (z_sc > thr1)
		{
			major_set.kmer.push_back(i -> first);
			major_set.occurr.push_back(tmp_cover1);
			major_set.cover1.push_back(tmp_cover1);
			major_set.cover2.push_back(tmp_cover2);
			major_set.site.push_back(i -> second);
			major_set.score.push_back(z_sc);
		}
		else if (z_sc > thr3)
		{
			minor_set.kmer.push_back(i -> first);
			minor_set.occurr.push_back(tmp_cover1);
			minor_set.cover1.push_back(tmp_cover1);
			minor_set.cover2.push_back(tmp_cover2);
			minor_set.site.push_back(i -> second);
			minor_set.score.push_back(z_sc);
		}
	}

	cout << "there are altogether " << major_set.kmer.size() << 
		" significant " << kmer_length << "-mers left.\n";
	cout << "there are altogether " << minor_set.kmer.size() << 
		" sub-significant " << kmer_length << "-mers left.\n";
}

/*************************************************************************************************
 *
 * Function get_rev used to get the reverse strand of a string
 *
 ************************************************************************************************/

string get_rev(string str)
{
	string rev_str;	// The reverse string

	for (int i = str.length() - 1; i >= 0; i --)
	{
		rev_str += alphabet[3 - r_alphabet[str[i]]];
	}

	return rev_str;
}

/*************************************************************************************************
 *
 * Function check_parlindrome used to check whether the string is palindrome
 *
 ************************************************************************************************/

int check_parlindrome(string str)
{
	string rev_str = get_rev(str);	// The reverse string

	int flag = 0;

	if (str == rev_str)
	{
		flag = 1;
	}

	return flag;
}

/*************************************************************************************************
 *
 * Function reverse_st used to reverse the occurring sites of reverse supplementary strands
 *
 ************************************************************************************************/

void reverse_st(vector<int>& positive, vector<int>& negative, 
		vector<string>& strand)
{
	negative = positive;

	for (int i = 0; i < negative.size() - 1; i += 2)
	{
		negative[i]++;
	}

	for(int i=1; i<negative.size(); i+=2)
	{
		negative[i]=strand[negative[i-1]].length()-negative[i]-kmer_length;
	}
}

/*************************************************************************************************
 *
 * Function comb_kmer used to combine reverse complementary k-mers
 *
 ************************************************************************************************/

void comb_kmer(map<string, int_vector>& kmer_site, map<string, int>& fa_flag, vector<string>& positive, map<string, int_vector>& kmer_bg, map<string, int>& bg_flag, vector<string>& negative)
{
	for(map<string, int_vector>::iterator it=kmer_site.begin(); it!=kmer_site.end(); it++)
	{
		fa_flag[it->first]=1;
	}

	for(map<string, int_vector>::iterator it=kmer_site.begin(); it!=kmer_site.end(); it++)
	{
		if(fa_flag[it->first]==0 || check_parlindrome(it->first))
		{
			continue;
		}

		string rev_str=get_rev(it->first);

		if(kmer_site.count(rev_str)>0)
		{
			fa_flag[rev_str]=0;
			vector<int> rev_site;
			rev_site.reserve(kmer_site[rev_str].size());
			reverse_st(kmer_site[rev_str], rev_site, positive);
			(it->second).insert((it->second).end(), rev_site.begin(), rev_site.end());
			kmer_site[rev_str].clear();
		}
	}

	bg_flag=fa_flag;

	for(map<string, int>::iterator it=bg_flag.begin(); it!=bg_flag.end(); it++)
	{
		if(it->second == 0)
		{
			continue;
		}

		string rev_str=get_rev(it->first);

		if(kmer_bg.count(rev_str)>0)
		{
			vector<int> rev_site;
			rev_site.reserve(kmer_bg[rev_str].size());
			reverse_st(kmer_bg[rev_str], rev_site, negative);
			kmer_bg[it->first].insert(kmer_bg[it->first].end(), rev_site.begin(), rev_site.end());
			kmer_bg[rev_str].clear();
		}
	}
}

/*************************************************************************************************
 *
 * Function kmer_sort used to sort the chosen k-mers
 *
 ************************************************************************************************/

void kmer_sort(kmer_set& chosen_kmers)
{
	float *kmer_score;
	int *kmer_id1;
	int num=chosen_kmers.score.size();
	kmer_score=new float [num];
	kmer_id1=new int [num];
	int i;

	for(i=0; i<num; i++)
	{
		kmer_score[i]=chosen_kmers.score[i];
		kmer_id1[i]=i;
	}

	quickSort(kmer_score, 0, num-1, kmer_id1);
	reverse_sort(kmer_id1, num);

	kmer_set temp_kmers;
	for(i=0; i<num; i++)
	{
		temp_kmers.kmer.push_back(chosen_kmers.kmer[kmer_id1[i]]);
		temp_kmers.score.push_back(chosen_kmers.score[kmer_id1[i]]);
		temp_kmers.cover1.push_back(chosen_kmers.cover1[kmer_id1[i]]);
		temp_kmers.cover2.push_back(chosen_kmers.cover2[kmer_id1[i]]);
		temp_kmers.occurr.push_back(chosen_kmers.occurr[kmer_id1[i]]);
		temp_kmers.site.push_back(chosen_kmers.site[kmer_id1[i]]);
	}

	chosen_kmers=temp_kmers;

	delete [] kmer_id1;
	delete [] kmer_score;

	cout<<"Finished sorting the significant "<<kmer_length<<"-mers based on z-scores.\n";
}

/*************************************************************************************************
 *
 * Function psm_sort used to sort all the PSMs and simultaneously delete the trivial PSMs
 *
 ************************************************************************************************/

void psm_sort(vector<psm>& psm_set)
{
	int *psm_id1;
	float *psm_scores;
	int num=psm_set.size();
	psm_scores=new float [num];
	psm_id1=new int [num];

	int i;

	for(i=0; i<num; i++)
	{
		psm_scores[i]=psm_set[i].score;
		psm_id1[i]=i;
	}

	quickSort(psm_scores, 0, num-1, psm_id1);
	reverse_sort(psm_id1, num);

	{
		vector<psm> tmp_psm;
		for(i=0; i<num; i++)
		{
			tmp_psm.push_back(psm_set[psm_id1[i]]);
		}

		psm_set.swap(tmp_psm);
	}

	delete [] psm_id1;
	delete [] psm_scores;

	cout<<"Finished sorting the PSMs in PSM set.\n";
	cout<<"Finished kicking out the trivial PSMs.\n"<<"There're "<<psm_set.size()<<" PSMs left."<<endl;
}

/*************************************************************************************************
 *
 * Function generate_flag used to generate flags
 *
 ************************************************************************************************/

void generate_flag(int length, vector<int>& flag)
{
	int i;

	flag.reserve(length);

	for(i=0; i<length; i++)
	{
		flag.push_back(0);
	}
}

/*************************************************************************************************
 *
 * Function mask_kmer used to mask the selected k-mers
 *
 ************************************************************************************************/

void mask_kmer(vector<int>& set, vector<int>& kmer_flag)
{
	int i;

	for(i=0; i<set.size(); i++)
	{
		kmer_flag[set[i]]=1;
	}
}

/*************************************************************************************************
 *
 * Function find_seed used to find the id of next seed for partitioning k-mers
 *
 ************************************************************************************************/

int find_seed(vector<int>& kmer_flag)
{
	int out_id=-1;
	int i;

	for(i=0; i<kmer_flag.size(); i++)
	{
		if(kmer_flag[i]==0)
		{
			out_id=i;
			break;
		}
	}

	return out_id;
}

/*************************************************************************************************
 *
 * Function comp_hd used to get the Hamming Distance between two strings
 *
 ************************************************************************************************/

int comp_hd(string str1, string str2)
{
	if(str1.length()!=str2.length())
	{
		cerr<<"Error: only two strings with the same length can be compared to get Hamming Distance.\n";
		exit(0);
	}

	int hd=0;

	for(int i=0; i<str1.length(); i++)
	{
		if(str1[i]!=str2[i])
		{
			hd++;
		}
	}

	return hd;
}

/*************************************************************************************************
 *
 * Function kmer2psm used to transform all significant k-mers into PSMs
 *
 ************************************************************************************************/

void kmer2psm(vector<kmer_all>& kmer_final, vector<psm>& psm_set, int& top_num)
{
	psm_set.clear();
	psm tmp_psm;

	for(int i=0; i<top_num; i++)
	{
		tmp_psm.set.clear();
		tmp_psm.set.push_back(i);
		psm_set.push_back(tmp_psm);
	}

	for(int i=1; i<top_num; i++)
	{
		for(int j=0; j<i; j++)
		{
			if(comp_hd(kmer_final[i].kmer, kmer_final[j].kmer)<=hd_thr)
			{
				psm_set[i].set.push_back(j);
				psm_set[j].set.push_back(i);
			}
		}

		for(int j=top_num; j<kmer_final.size(); j++)
		{
			if(comp_hd(kmer_final[i].kmer, kmer_final[j].kmer)<=hd_thr)
			{
				psm_set[i].set.push_back(j);
			}
		}
	}

	for(int i=0; i<psm_set.size(); i++)
	{
		psm_set[i].occurr=0;

		for(int j=0; j<psm_set[i].set.size(); j++)
		{
			psm_set[i].occurr+=kmer_final[psm_set[i].set[j]].occurr;
		}
	}

	cout<<"Finished constructing PSMs.\n";
	cout<<"There are altogether "<<psm_set.size()<<" PSMs constructed.\n";
}

/*************************************************************************************************
 *
 * Function test_psm used to compute the z-scores of all PSMs
 *
 ************************************************************************************************/

void test_psm(vector<psm>& psm_set, float& sum1, float& sum2)
{
	for(int i=0; i<psm_set.size(); i++)
	{
		psm_set[i].score=ztest(psm_set[i].cover1, psm_set[i].cover2, sum1, sum2);
	}


	cout<<"Finished computing the z-scores for all PSMs\n";
}

/*************************************************************************************************
 *
 * Function init_mat used to initiate the PSM matrix
 *
 ************************************************************************************************/

void init_mat(vector<float_vector>& psm_mat, int& length)
{
	for(int i=0; i<psm_mat.size(); i++)
	{
		psm_mat[i].clear();
	}
	psm_mat.clear();

	vector<float> tmp_vec;

	for(int i=0; i<4; i++)
	{
		tmp_vec.clear();
		for(int j=0; j<length; j++)
		{
			tmp_vec.push_back(0);
		}
		psm_mat.push_back(tmp_vec);
	}
}

/*************************************************************************************************
 *
 * Function comp_mat used to compute the PSM matrix based on k-mer set of PSM
 *
 ************************************************************************************************/

void comp_mat(vector<float_vector>& psm_mat, vector<int>& psm_set, vector<kmer_all>& kmer_final, int& length)
{
	int i, j;

	init_mat(psm_mat, length);

	for(i=0; i<psm_set.size(); i++)
	{
		for(j=0; j<length; j++)
		{
			psm_mat[r_alphabet[kmer_final[psm_set[i]].kmer[j]]][j]+=kmer_final[psm_set[i]].occurr;
		}
	}
}

/*************************************************************************************************
 *
 * Function norm_mat used to normalize the matrices of the PSMs
 *
 ************************************************************************************************/

void norm_mat(vector<float_vector>& tmp_mat, int& occurr, int& length)
{
	for(int i=0; i<4; i++)
	{
		for(int j=0; j<length; j++)
		{
			tmp_mat[i][j]/=occurr;
		}
	}
}

/*************************************************************************************************
 *
 * Function fill_psm used to compute the matrices for all PSMs
 *
 ************************************************************************************************/

void fill_psm(vector<psm>& psm_set, vector<kmer_all>& kmer_final, int& length)//, map<char, int>& r_alphabet)
{
	int i;

	for(i=0; i<psm_set.size(); i++)
	{
		comp_mat(psm_set[i].mat, psm_set[i].set, kmer_final, length);
		norm_mat(psm_set[i].mat, psm_set[i].occurr, length);
	}

	cout<<"Finished computing the matrices information of each PSM.\n";
}

/*************************************************************************************************
 *
 * Function init_str used to initiate the information of the consensus string
 *
 ************************************************************************************************/

void init_str(string& str, int length)
{
	str.assign(length, 'A');
}

/*************************************************************************************************
 *
 * Function max_col used to choose the element with the maximum value in this column
 *
 ************************************************************************************************/

int max_col(vector<float_vector>& mat, int col)
{
	float max_val=0;
	int max_id=0;

	for(int i=0; i<mat.size(); i++)
	{
		if(mat[i][col]>max_val)
		{
			max_val=mat[i][col];
			max_id=i;
		}
	}

	return max_id;
}

/*************************************************************************************************
 *
 * Function get_cons used to compute the consensus string for the current PSM
 *
 ************************************************************************************************/

void get_cons(string& cons, vector<float_vector>& mat)
{
	int i;

	init_str(cons, mat[0].size());

	for(i=0; i<cons.length(); i++)
	{
		cons[i]=alphabet[max_col(mat, i)];
	}
}

/*************************************************************************************************
 *
 * Function comp_cons used to compute the consensus strings of the PSMs
 *
 ************************************************************************************************/

void comp_cons(vector<psm>& psm_set)
{
	int i;

	for(i=0; i<psm_set.size(); i++)
	{
		get_cons(psm_set[i].cons, psm_set[i].mat);
	}

	cout<<"Finished computing the consensus strings for each PSM.\n";
}

/*************************************************************************************************
 *
 * Function add_str used to add one flanking segment to the matrix
 *
 ************************************************************************************************/

void add_str(vector<float_vector>& tmp_mat, string str)
{
	for(int i=0; i<str.size(); i++)
	{
		tmp_mat[r_alphabet[str[i]]][i]++;
	}
}

/*************************************************************************************************
 *
 * Function get_flank used to get the left and right flanking matrices
 *
 ************************************************************************************************/

void get_flank(vector<float_vector>& tmp_mat, vector<int>& site, sequence& seq_final, int direct, int id)
{
	init_mat(tmp_mat, lmer_length);

	for(int i=0; i<site.size(); i+=2)
	{
		if(direct==-1)
		{
			add_str(tmp_mat, seq_final.strand[site[i]].substr(site[i+1]-lmer_length, lmer_length));
		}
		else if(direct==1)
		{
			add_str(tmp_mat, seq_final.strand[site[i]].substr(site[i+1]+kmer_length, lmer_length));
		}
		else
		{
			cerr<<"Error: There's errors in function get_flank. Please check the parameters.\n";
			exit(0);
		}
	}
}

/*************************************************************************************************
 *
 * Function combine_major used to combine two sets
 *
 ************************************************************************************************/

int combine_major(kmer_set& major_set, kmer_set& minor_set)
{
	int res=major_set.kmer.size();
	major_set.kmer.insert(major_set.kmer.end(), minor_set.kmer.begin(), minor_set.kmer.end());
	major_set.occurr.insert(major_set.occurr.end(), minor_set.occurr.begin(), minor_set.occurr.end());
	major_set.cover1.insert(major_set.cover1.end(), minor_set.cover1.begin(), minor_set.cover1.end());
	major_set.cover2.insert(major_set.cover2.end(), minor_set.cover2.begin(), minor_set.cover2.end());
	major_set.score.insert(major_set.score.end(), minor_set.score.begin(), minor_set.score.end());
	major_set.site.insert(major_set.site.end(), minor_set.site.begin(), minor_set.site.end());

	return res;
}
	
/*************************************************************************************************
 *
 * Function fill_kmer used to obtain the k-mer information as well as their flanking segments
 *
 ************************************************************************************************/

void fill_kmer(kmer_set& major_set, vector<kmer_all>& kmer_final, int seq_num, sequence& seq_final)
{
	kmer_all tmp_kmer;

	for(int i=0; i<major_set.kmer.size(); i++)
	{
		tmp_kmer.kmer=major_set.kmer[i];
		tmp_kmer.occurr=major_set.occurr[i];
		tmp_kmer.cover1=major_set.cover1[i];
		tmp_kmer.cover2=major_set.cover2[i];
		tmp_kmer.site=major_set.site[i];
		tmp_kmer.score=major_set.score[i];
		get_flank(tmp_kmer.left, major_set.site[i], seq_final, -1, i);
		get_flank(tmp_kmer.right, major_set.site[i], seq_final, 1, i);
		kmer_final.push_back(tmp_kmer);
	}

	cout<<"Finished calculating the profile matrices by the flanking segments beside the "<<kmer_length<<"-mers.\n";
}

/*************************************************************************************************
 *
 * Function free_kmer used to free the memory of kmer_set
 *
 ************************************************************************************************/

void free_kmer(kmer_set& major_set)
{
	for(int i=0; i<major_set.kmer.size(); i++)
	{
		string().swap(major_set.kmer[i]);
		vector<int>().swap(major_set.site[i]);
	}
	vector<string>().swap(major_set.kmer);
	vector<int_vector>().swap(major_set.site);
}

/*************************************************************************************************
 *
 * Function comp_sw used to check whether two vectors have overlap
 *
 ************************************************************************************************/

float comp_sw(vector<float_vector>& mat1, vector<float_vector>& mat2)
{
	if(mat1[0].size() != mat2[0].size())
	{
		cerr<<"Error: The two matrices are not equal in width!\n";
		exit(1);
	}

	int length = mat1[0].size();
	float sw=0;

	for(int j = 0; j < length; j++)
	{
		for(int i=0; i<4; i++)
		{
			sw+=(mat1[i][j] - mat2[i][j]) * 
				(mat1[i][j] - mat2[i][j]);
		}
	}

	sw/=length;
	sw=2-sw;

	return sw;
}

/*************************************************************************************************
 *
 * Function comp_overlap used to check whether two PSMs have overlap
 *
 ************************************************************************************************/

float comp_overlap(vector<int>& set1, vector<int>& set2)
{
	int flag=0;

	for(int i=0; i<set1.size(); i++)
	{
		if(find(set2.begin(), set2.end(), set1[i])
				!=set2.end() && set1[i]<top_num)
		{
			flag=1;
			break;
		}
	}

	return flag;
}

/*************************************************************************************************
 *
 * Function build_graph used to build similarity graph of PSMs
 *
 ************************************************************************************************/

void build_graph(vector<psm>& psm_set, vector<int_vector>& psm_nbr)
{
	vector<int> tmp_vec;

	for(int i=0; i<psm_set.size(); i++)
	{
		psm_nbr.push_back(tmp_vec);
	}

	for(int i=1; i<psm_set.size(); i++)
	{
		for(int j=0; j<i; j++)
		{
			float tmp_sw=comp_sw(psm_set[i].mat, psm_set[j].mat);
			if(tmp_sw>=sw_thr)
			{
				psm_nbr[i].push_back(j);
				psm_nbr[j].push_back(i);
			}
		}
	}

	cout<<"Finished building PSM Adjacency Graph.\n";
	cout<<"There're altogether "<<psm_nbr.size()<<" nodes in this graph.\n";
}

/************************************************************************************************
 *
 * Function psm2kmer used to transform PSM set into k-mer set
 *
 ************************************************************************************************/

void psm2kmer(vector<int>& tmp_kmer, vector<int>& tmp_psm, vector<psm>& all_psm, int kmer_num)
{
	tmp_kmer.clear();
	vector<int> flag;
	generate_flag(kmer_num, flag);

	for(int i=0;i<tmp_psm.size(); i++)
	{
		mask_kmer(all_psm[tmp_psm[i]].set, flag);
	}

	for(int i=0; i<kmer_num; i++)
	{
		if(flag[i]==1)
		{
			tmp_kmer.push_back(i);
		}
	}
}

/*************************************************************************************************
 *
 * Function comp_cover used to calculate the approximate coverage of the current k-mer set
 *
 ************************************************************************************************/

float comp_cover(vector<kmer_all>& kmer_final, vector<int>& tmp_set, int flag)
{
	int tmp_cover=0;

	if(flag==1)
	{
		for(int i=0; i<tmp_set.size(); i++)
		{
			tmp_cover+=kmer_final[tmp_set[i]].cover1;
		}
	}
	else if(flag==-1)
	{
		for(int i=0; i<tmp_set.size(); i++)
		{
			tmp_cover+=kmer_final[tmp_set[i]].cover2;
		}
	}

	return tmp_cover;
}

/*************************************************************************************************
 *
 * Function cover_psm used to compute the coverage of PSMs
 *
 ************************************************************************************************/

void cover_psm(vector<psm>& psm_set, vector<kmer_all>& kmer_final)
{
	for(int i=0; i<psm_set.size(); i++)
	{
		psm_set[i].cover1=comp_cover(kmer_final, psm_set[i].set, 1);
		psm_set[i].cover2=comp_cover(kmer_final, psm_set[i].set, -1);
	}

	cout<<"Finished cumputing the coverage for all PSMs.\n";
}

/*************************************************************************************************
 *
 * Function comp_occurr used to calculate the occurrence of the current k-mer set
 *
 ************************************************************************************************/

int comp_occurr(vector<kmer_all>& kmer_final, vector<int>& tmp_set)
{
	int tmp_occurr=0;

	for(int i=0; i<tmp_set.size(); i++)
	{
		tmp_occurr+=kmer_final[tmp_set[i]].occurr;
	}

	return tmp_occurr;
}

/*************************************************************************************************
 *
 * Function sw_col used to calculate the SW similarity score between two columns
 *
 ************************************************************************************************/

float sw_col(vector<float_vector>& mat1, vector<float_vector>& mat2, int col_id)
{
	float result=2;

	for(int i=0; i<4; i++)
	{
		result-=(mat1[i][col_id]-mat2[i][col_id])*(mat1[i][col_id]-mat2[i][col_id]);
	}

	return result;
}

/*************************************************************************************************
 *
 * Function sw_sim used to calculate the SW similarity scores
 *
 ************************************************************************************************/

float sw_sim(vector<float_vector>& mat1, vector<float_vector>& mat2)
{
	float result=0;

	for(int j=0; j<mat1[0].size(); j++)
	{
		result+=sw_col(mat1, mat2, j);
	}

	return result;
}

/*************************************************************************************************
 *
 * Function get_rand used to get one random index based on one distribution
 *
 ************************************************************************************************/

int get_rand(vector<float>& sim_vec)
{
	float* p;
	p=new float [sim_vec.size()];

	for(int i=0; i<sim_vec.size(); i++)
	{
		p[i]=sim_vec[i];
	}

	int result=genmulone(p, sim_vec.size());

	free_pt(p);

	return result;
}

/*************************************************************************************************
 *
 * Function sam_node used to choose the new node to initiate the Gibbs Sampling
 *
 ************************************************************************************************/

int sam_node(vector<int>& tmp_pre, vector<float_vector>& tmp_mat, vector<psm>& psm_set)
{
	int result;

	vector<float> sim_vec;	// vector to hold the similarity scores
	sim_vec.reserve(tmp_pre.size());

	for(int i=0; i<tmp_pre.size(); i++)
	{
		sim_vec.push_back(sw_sim(tmp_mat, psm_set[tmp_pre[i]].mat));
	}

	result=0;
	float max_sim=sim_vec[0];
	for(int i=1; i<sim_vec.size(); i++)
	{
		if(max_sim<sim_vec[i])
		{
			max_sim=sim_vec[i];
			result=i;
		}
	}

	return tmp_pre[result];
}

/*************************************************************************************************
 *
 * Function choose_node used to choose the old node to initiate the Gibbs Sampling
 *
 ************************************************************************************************/

int choose_node(vector<int>& tmp_pre)
{
	float* sam_vec;
	sam_vec=new float [tmp_pre.size()];

	for(int i=0; i<tmp_pre.size(); i++)
	{
		sam_vec[i]=1/(float)tmp_pre.size();
	}

	int res;
	res=genmulone(sam_vec, tmp_pre.size());
	res=tmp_pre[res];

	free_pt(sam_vec);
	return res;
}

/*************************************************************************************************
 *
 * Function get_nbr used to get the neighbor list of old node
 *
 ************************************************************************************************/

void get_nbr(int& old_node, vector<int>& psm_flag, vector<int>& tmp_nbr, vector<int_vector>& psm_nbr, vector<int>& tmp_pre)
{
	tmp_nbr.clear();

	for(int i=0; i<psm_nbr[old_node].size(); i++)
	{
		if(psm_flag[psm_nbr[old_node][i]]==0 && find(tmp_pre.begin(), tmp_pre.end(), psm_nbr[old_node][i])!=tmp_pre.end())
		{
			tmp_nbr.push_back(psm_nbr[old_node][i]);
		}
	}
}

/*************************************************************************************************
 *
 * Function comp_sc used to calculate the motif score of the current preliminary motif
 *
 ************************************************************************************************/

float comp_sc(vector<float_vector> mat, int& occurr, int& length)
{
	float sc=0;

	for(int i=0; i<4; i++)
	{
		for(int j=0; j<length; j++)
		{
			if(mat[i][j]==0)
			{
				mat[i][j]=0.00001;
			}

			mat[i][j]=mat[i][j]*log(mat[i][j]/nt_freq[i]);
			sc+=mat[i][j];
		}
	}

	sc/=length;
	sc=exp(sc);
	sc=sc*occurr;

	return sc;
}

/*************************************************************************************************
 *
 * Function del_element used to delete element in vector
 *
 ************************************************************************************************/

void del_element(vector<int>& tmp_vec, int& val)
{
	for(vector<int>::iterator i=tmp_vec.begin(); i!=tmp_vec.end(); i++)
	{
		if(*i==val)
		{
			tmp_vec.erase(i);
			break;
		}
	}
}

/*************************************************************************************************
 *
 * Function if_accept used to check whether we add one node to the current preliminary motif
 *
 ************************************************************************************************/

int if_accept(vector<int>& tmp_pre, int& tmp_occurr, vector<float_vector>& tmp_mat, int& old_node, int& new_node, vector<kmer_all>& kmer_final, vector<psm>& psm_set, int seq_num)
{
	vector<int> old_pre=tmp_pre;
	old_pre.push_back(old_node);
	vector<int> old_kmer;
	psm2kmer(old_kmer, old_pre, psm_set, kmer_final.size());
	vector<float_vector> old_mat;
	int old_occurr=comp_occurr(kmer_final, old_kmer);
	comp_mat(old_mat, old_kmer, kmer_final, kmer_length);
	norm_mat(old_mat, old_occurr, kmer_length);
	float old_score=comp_sc(old_mat, old_occurr, kmer_length);

	vector<int> new_pre=old_pre;
	new_pre.push_back(old_node);
	new_pre.push_back(new_node);
	vector<int> new_kmer;
	psm2kmer(new_kmer, new_pre, psm_set, kmer_final.size());
	int new_occurr=comp_occurr(kmer_final, new_kmer);
	vector<float_vector> new_mat;
	comp_mat(new_mat, new_kmer, kmer_final, kmer_length);
	norm_mat(new_mat, new_occurr, kmer_length);
	float new_score=comp_sc(new_mat, new_occurr, kmer_length);

	vector<int> del_pre=old_pre;
	del_pre.push_back(new_node);
	vector<int> del_kmer;
	psm2kmer(del_kmer, del_pre, psm_set, kmer_final.size());
	int del_occurr=comp_occurr(kmer_final, del_kmer);
	vector<float_vector> del_mat;
	comp_mat(del_mat, del_kmer, kmer_final, kmer_length);
	norm_mat(del_mat, del_occurr, kmer_length);
	float del_score=comp_sc(del_mat, del_occurr, kmer_length);

    vector<double> v_score(old_score, del_score, new_score);
    double sum = accumulate(v_score.begin(), v_score.end(), 0); // get the sum of the vector
    
    for (int i = 0; i < v_score.size(); i ++)
    {
        v_score[i] /= sum;
    }
    
    return(genmulone(v_score)); // randomly sample an id according to the probability distributation
    // calculated based on the motif score
}

/*************************************************************************************************
 *
 * Function psm2pre used to transform PSMs into preliminary motifs
 *
 ************************************************************************************************/

void psm2pre(vector<int_vector>& pre_out, vector<psm>& psm_set)
{
	pre_out.clear();

	for(int i=0; i<psm_set.size(); i++)
	{
		pre_out.push_back(psm_set[i].set);
	}
}

/*************************************************************************************************
 *
 * Function gibbs_sam used to perform Gibbs Sampling
 *
 ************************************************************************************************/

void gibbs_sam(vector<int_vector>& pre_out, vector<psm>& psm_set, 
		vector<kmer_all>& kmer_final, int& num_pre, int& num_iter, 
		vector<int_vector>& psm_nbr, int seq_num)
{
	vector<int> psm_flag;
	generate_flag(psm_set.size(), psm_flag);
	int seed_psm;
	vector<int> tmp_pre, tmp_kmer, tmp_nbr;
	vector<float_vector> tmp_mat;
	int tmp_occurr, old_node, new_node, del_node;
	float tmp_score, tmp_cover1, tmp_cover2;
	int sam_flag;

	cout<<"Begin Gibbs Sampling ...\n";

	for(int i=1; i<=num_pre; i++)
	{
		seed_psm=find_seed(psm_flag);

		if(seed_psm==-1)
		{
			break;
		}

		tmp_pre.clear();
		tmp_pre.push_back(seed_psm);
		for(int j=0; j<psm_nbr[seed_psm].size(); j++)
		{
			if(psm_flag[psm_nbr[seed_psm][j]]==0)
			{
				tmp_pre.push_back(psm_nbr[seed_psm][j]);
			}
		}
		mask_kmer(tmp_pre, psm_flag);
		psm2kmer(tmp_kmer, tmp_pre, psm_set, kmer_final.size());

		if(tmp_pre.size()==1)
		{
			pre_out.push_back(tmp_kmer);
			continue;
		}

		for(int j=1; j<=num_iter; j++)
		{
			old_node=choose_node(tmp_pre);
			del_element(tmp_pre, old_node);
			psm_flag[old_node]=0;
			psm2kmer(tmp_kmer, tmp_pre, psm_set, kmer_final.size());
			tmp_occurr=comp_occurr(kmer_final, tmp_kmer);
			comp_mat(tmp_mat, tmp_kmer, kmer_final, kmer_length);
			norm_mat(tmp_mat, tmp_occurr, kmer_length);
			tmp_score=comp_sc(tmp_mat, tmp_occurr, kmer_length);
			get_nbr(old_node, psm_flag, tmp_nbr, psm_nbr, tmp_pre);

			if(tmp_nbr.size() == 0)
			{
				psm_flag[old_node]=1;
				tmp_pre.push_back(old_node);
				continue;
			}

			new_node = sam_node(tmp_nbr, tmp_mat, psm_set);

			sam_flag = if_accept(tmp_pre, tmp_occurr, tmp_mat, old_node, 
					new_node, kmer_final, psm_set, seq_num);
			
			if(sam_flag == 2)
			{
				tmp_pre.push_back(old_node);
				tmp_pre.push_back(new_node);
			}
			else if(sam_flag == 1)
			{
				tmp_pre.push_back(new_node);
			}
			else
			{
				tmp_pre.push_back(old_node);
			}

			mask_kmer(tmp_pre, psm_flag);
		}

		pre_out.push_back(tmp_kmer);
	}

	cout<<"Finished Gibbs Sampling.\n"<<"We have obtained altogether "<<pre_out.size()<<" preliminary motifs."<<endl;
}

/*************************************************************************************************
 *
 * Function build_flank used to build flanking matrices
 *
 ************************************************************************************************/

void build_flank(vector<int>& tmp_set, vector<kmer_all>& kmer_final, vector<float_vector>& flank, int direct, int length)
{
	init_mat(flank, length);

	for(int i=0; i<tmp_set.size(); i++)
	{
		if(direct==-1)
		{
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<length; k++)
				{
					flank[j][k]+=kmer_final[tmp_set[i]].left[j][k];
				}
			}
		}
		else
		{
			for(int j=0; j<4; j++)
			{
				for(int k=0; k<length; k++)
				{
					flank[j][k]+=kmer_final[tmp_set[i]].right[j][k];
				}
			}
		}
	}
}

/*************************************************************************************************
 *
 * Function site_pre used to combine the site information of all k-mers contained
 *
 ************************************************************************************************/

void site_pre(pre_mtf& tmp_pre, vector<kmer_all>& kmer_final)
{
	tmp_pre.site.clear();
	pair<int, int> tmp_pair;
	map<pair<int, int>, int> site_map;

	for(int i=0; i<tmp_pre.set.size(); i++)
	{
		for(int j=0; j<kmer_final[tmp_pre.set[i]].site.size(); j+=2)
		{
			tmp_pair=make_pair(kmer_final[tmp_pre.set[i]].site[j], 
					kmer_final[tmp_pre.set[i]].site[j+1]);
			site_map[tmp_pair]=1;
		}
	}

	for(map<pair<int, int>, int>::iterator it=site_map.begin(); it!=site_map.end(); it++)
	{
		tmp_pre.site.push_back((it->first).first);
		tmp_pre.site.push_back((it->first).second);
	}
}

/*************************************************************************************************
 *
 * Function get_pre used to get the set of preliminary motifs based on preliminary output from Gibbs Sampling
 *
 ************************************************************************************************/

void get_pre(vector<pre_mtf>& all_pre, vector<int_vector>& pre_out, 
		vector<kmer_all>& kmer_final, int seq_num)
{
	for(int i=0; i<pre_out.size(); i++)
	{
		pre_mtf tmp_pre;
		
		tmp_pre.set=pre_out[i];
		site_pre(tmp_pre, kmer_final);
		tmp_pre.occurr=comp_occurr(kmer_final, tmp_pre.set);
		comp_mat(tmp_pre.mat, tmp_pre.set, kmer_final, kmer_length);
		norm_mat(tmp_pre.mat, tmp_pre.occurr, kmer_length);
		build_flank(tmp_pre.set, kmer_final, tmp_pre.left, -1, lmer_length);
		build_flank(tmp_pre.set, kmer_final, tmp_pre.right, 1, lmer_length);
		norm_mat(tmp_pre.left, tmp_pre.occurr, lmer_length);
		norm_mat(tmp_pre.right, tmp_pre.occurr, lmer_length);
		tmp_pre.cover1=comp_cover(kmer_final, tmp_pre.set, 1);
		tmp_pre.cover2=comp_cover(kmer_final, tmp_pre.set, -1);
		tmp_pre.score=ztest(tmp_pre.cover1, tmp_pre.cover2, sum_exp, sum_bg);

		all_pre.push_back(tmp_pre);
	}
}

/*************************************************************************************************
 *
 * Function test_cell used to check one cell in one column (Two Proportion z-test)
 *
 ************************************************************************************************/

float test_cell(float this_val, float that_val, int occurr)
{
	float pos_occurr=occurr*this_val;
	float neg_occurr=occurr*that_val;
	float p=(pos_occurr+neg_occurr)/(2*occurr);
	float z=fabs((this_val-that_val)/sqrt(p*(1-p)*(2/(float)occurr)));
	
	return z;
}

/*************************************************************************************************
 *
 * Function test_col used to check one column (Two Proportion z-test)
 *
 ************************************************************************************************/

float test_col(vector<float_vector>& flank, int col_id, int occurr)
{
	float max_z=0;
	float tmp_z;

	for(int i=0; i<4; i++)
	{
		tmp_z=test_cell(flank[i][col_id], nt_freq[i], occurr);

		if(tmp_z>max_z)
		{
			max_z=tmp_z;
		}
	}

	return max_z;
}

/*************************************************************************************************
 *
 * Function refine_one used to refine one preliminary motif
 *
 ************************************************************************************************/

void refine_one(pre_mtf& one_pre)
{
	if(lmer_length <= 0)
	{
		one_pre.begin=0;

		return;
	}

	one_pre.begin=lmer_length-1;

	for(int i=lmer_length-1; i>=0; i--)
	{
		if(test_col(one_pre.left, i, one_pre.occurr)<thr2)
		{
			one_pre.begin=i+1;
			break;
		}
	}

//	cout << one_pre.begin << "\n";

	one_pre.end=0;

	for(int i=0; i<lmer_length; i++)
	{
		if(test_col(one_pre.right, i, one_pre.occurr)<thr2)
		{
			one_pre.end=i-1;
			break;
		}
	}

//	cout << one_pre.begin << " " << one_pre.end << "\n";
}

/*************************************************************************************************
 *
 * Function refine_mtf used to refine the preliminary motifs
 *
 ************************************************************************************************/

void refine_mtf(vector<pre_mtf>& all_pre)
{
	for(int i=0; i<all_pre.size(); i++)
	{
		refine_one(all_pre[i]);

//		cout << all_pre[i].begin << " " << 
//			all_pre[i].end << "\n";
	}
}

/*************************************************************************************************
 *
 * Function sift_pvalue used to sift preliminary motifs using chi square p-values
 *
 ************************************************************************************************/

vector<pre_mtf> sift_pvalue(vector<pre_mtf>& all_pre, vector<kmer_all>& kmer_final, 
    trans_mat& T)
{
    vector<pre_mtf> out_pre; // declare a vector to hold retained preliminary motifs
    out_pre.reserve(all_pre.size()); // reserve space
    
    for (int i = 0; i < all_pre.size(); i ++) // for each preliminary motif
    {
        int n = 0; // initialize the counting number
        // calculate the chi square p-value
        double chi_val = 0; // initialize the chi-square static
        
        for (int j = 0; j < all_pre[i].set.size(); j ++) // for each k-mer in the preliminary motif
        {
            n ++; // update the number of k-mers
            P = T.A[kmer_final[all_pre[i].set[j]][0]] * T.B[kmer_final[all_pre[i].set[j]].substr(0, 2)] * 
                T.C[kmer_final[all_pre[i].set[j]].substr(0, 3)]; // probability of the first three elements
                // in the string
            
            for (k = 3; k < all_pre[i].set[j]].length(); k ++) // for the remaining elements
            {
                P *= T.D[kmer_final[all_pre[i].set[j]].substr(k - 3, k)]; // probabilities of 
                // the remaining elements
                // maybe I can logrithmize the p-value for acceleration
            }
            
            chi_val += pow((kmer_final[all_pre[i].set[j]].occurr / all_pre[i].occurr - P), 2) / P;
            // calculate the chi-square statistic
            
        }
    }
}

/*************************************************************************************************
 *
 * Function cut_mat used to transform matrices
 *
 ************************************************************************************************/

void cut_mat(pre_mtf& tmp_pre, vector<float_vector>& mat)
{
	mat.clear();
	deque<float> tmp_ln;
	vector<float> tmp_vec;

	for (int i = 0; i < 4; i ++)
	{
		tmp_ln.clear();
		tmp_vec.clear();

		for (int j = 0; j < kmer_length; j ++)
		{
			tmp_ln.push_back
				(tmp_pre.mat[i][j]);
		}

		if (lmer_length > 0)
		{
			if (tmp_pre.end >= 0)
			{
				for (int j = 0; j <= tmp_pre.end; j ++)
				{
					tmp_ln.push_back
						(tmp_pre.right[i][j]);
				}
			}
	
			if (tmp_pre.begin < lmer_length)
			{
				for (int j = lmer_length - 1; 
					j >= tmp_pre.begin; j --)
				{
					tmp_ln.push_front
					(tmp_pre.left[i][j]);
				}
			}
		}

		for (deque<float>::iterator i = tmp_ln.begin(); 
				i != tmp_ln.end(); i ++)
		{
			tmp_vec.push_back(*i);
		}

		mat.push_back(tmp_vec);
	}
}

/*************************************************************************************************
 *
 * Function make_wild used to make hash table of wild cards
 *
 ************************************************************************************************/

void make_wild(char* wild_card, char* reverse_wild)
{
	wild_card[0]='N';
	reverse_wild[0]='N';
	wild_card[1]='A';
	reverse_wild[1]='T';
	wild_card[2]='C';
	reverse_wild[2]='G';
	wild_card[3]='M';	// A, C
	reverse_wild[3]='K';
	wild_card[4]='G';
	reverse_wild[4]='C';
	wild_card[5]='R';	// A, G
	reverse_wild[5]='Y';
	wild_card[6]='S';	// C, G
	reverse_wild[6]='S';
	wild_card[7]='V';	// A, C, G
	reverse_wild[7]='B';
	wild_card[8]='T';
	reverse_wild[8]='A';
	wild_card[9]='W';	// A, T
	reverse_wild[9]='W';
	wild_card[10]='Y';	// C, T
	reverse_wild[10]='R';
	wild_card[11]='H';	// A, C, T
	reverse_wild[11]='D';
	wild_card[12]='K';	// G, T
	reverse_wild[12]='M';
	wild_card[13]='D';	// A, G, T
	reverse_wild[13]='H';
	wild_card[14]='B';	// G, C, T
	reverse_wild[14]='V';

	cout<<"Finished constructing wild card tables.\n";
}

/*************************************************************************************************
 *
 * Function col_deg used to get the degenerate consensus string 
 *
 ************************************************************************************************/

void col_deg(string& deg, string& rev_deg, vector<float_vector>& mat, int& col_id)
{
	int sum=0;

	for(int i=0; i<4; i++)
	{
		if(mat[i][col_id]>0.25)
		{
			sum+=(int)(rint(pow(2, i)));
		}
	}

	deg[col_id]=wild_card[sum];
	rev_deg[rev_deg.length()-1-col_id]=reverse_wild[sum];
}

/*************************************************************************************************
 *
 * Function get_deg used to get the degenerate consensus strings
 *
 ************************************************************************************************/

void get_deg(string& deg, string& rev_deg, vector<float_vector>& mat)
{
	init_str(deg, mat[0].size());
	init_str(rev_deg, mat[0].size());

	for(int i=0; i<deg.length(); i++)
	{
		col_deg(deg, rev_deg, mat, i);
	}
}

/*************************************************************************************************
 *
 * Function build_rev used to get the reverse consensus string
 *
 ************************************************************************************************/

void build_rev(string& rev_cons, string& cons)
{
	init_str(rev_cons, cons.length());

	for(int i=0; i<cons.length(); i++)
	{
		rev_cons[cons.length()-1-i]=alphabet[3-r_alphabet[cons[i]]];
	}
}

/*************************************************************************************************
 *
 * Function build_mtf used to transform one preliminary motif into one motif
 *
 ************************************************************************************************/

void build_mtf(pre_mtf& tmp_pre, mtf& tmp_mtf)
{
	tmp_mtf.set=tmp_pre.set;
	tmp_mtf.nsites=tmp_pre.occurr;
	cut_mat(tmp_pre, tmp_mtf.mat);
	tmp_mtf.alength=tmp_mtf.mat[0].size();
	get_cons(tmp_mtf.cons, tmp_mtf.mat);
	get_deg(tmp_mtf.deg, tmp_mtf.rev_deg, tmp_mtf.mat);
	build_rev(tmp_mtf.rev_cons, tmp_mtf.cons);
	tmp_mtf.score=tmp_pre.score;
}

/*************************************************************************************************
 *
 * Function get_site used to get the sites of motif
 *
 ************************************************************************************************/

void get_site(pre_mtf& tmp_pre, mtf& tmp_mtf, sequence& seq_final)
{
	st tmp_st;

	if(str_flag==1)
	{
		for(int i=0; i<tmp_pre.site.size(); i+=2)
		{
			tmp_st.header=seq_final.name[tmp_pre.site[i]];
			tmp_st.strand='+';
			tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin+1;
			tmp_st.end=tmp_pre.site[i+1]+kmer_length+tmp_pre.end+1;

			tmp_mtf.site.push_back(tmp_st);
		}
	}
	else
	{
		for(int i=0; i<tmp_pre.site.size(); i+=2)
		{
			tmp_st.header=seq_final.name[tmp_pre.site[i]/2];
			
			if(tmp_pre.site[i]%2==0)
			{
				tmp_st.strand='+';
//				cout << tmp_pre.site[i+1] << "\n";
				tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin+1;
				tmp_st.end=tmp_pre.site[i+1]+kmer_length+tmp_pre.end+1;
//				cout << tmp_st.begin << " " << lmer_length 
//					<< " " << tmp_pre.begin << "\n";
			}
			else
			{
				tmp_st.strand='-';
				tmp_st.begin=tmp_pre.site[i+1]-lmer_length+tmp_pre.begin;
				tmp_st.end=tmp_pre.site[i+1]+kmer_length-1+tmp_pre.end;
				int tmp_val=tmp_st.begin;
				tmp_st.begin=seq_final.strand[tmp_pre.site[i]].
					size()-tmp_st.end-1;
				tmp_st.end=seq_final.strand[tmp_pre.site[i]].size()-tmp_val;
			}

			tmp_mtf.site.push_back(tmp_st);
		}
	}
}

/*************************************************************************************************
 *
 * Function pre2mtf used to transform preliminary motifs into motifs
 *
 ************************************************************************************************/

void pre2mtf(vector<pre_mtf>& all_pre, vector<mtf>& all_mtf, sequence& seq_final)
{
	for(int i=0; i<all_pre.size(); i++)
	{
		mtf tmp_mtf;
		build_mtf(all_pre[i], tmp_mtf);
		get_site(all_pre[i], tmp_mtf, seq_final);
		all_mtf.push_back(tmp_mtf);
	}

	cout<<"Finished transforming preliminary motifs into motifs.\n";
}

/*************************************************************************************************
 *
 * Function sort_mtf used to sort all the motifs
 *
 ************************************************************************************************/

void sort_mtf(vector<mtf>& mtf_set)
{
		float *mtf_scores;
		int *mtf_id1;
		int num=mtf_set.size();
		mtf_scores=new float [num];
		mtf_id1=new int [num];

		int i;

		for(i=0; i<num; i++)
		{
				mtf_scores[i]=mtf_set[i].score;
				mtf_id1[i]=i;
		}

		quickSort(mtf_scores, 0, num-1, mtf_id1);
		reverse_sort(mtf_id1, num);

		vector<mtf> tmp_mtf;
		for(i=0; i<num; i++)
		{
			if(!last_trivial(mtf_set[mtf_id1[i]].cons))
			{
				continue;
			}

			tmp_mtf.push_back(mtf_set[mtf_id1[i]]);
		}

		mtf_set=tmp_mtf;

		delete [] mtf_id1;
		delete [] mtf_scores;

		cout << "Finished sorting the motifs according to motif scores.\n";
		cout << "Finished kicking out the trivial motifs.\nThere're " 
			<< mtf_set.size() << " motifs left.\n";
}

/*************************************************************************************************
 *
 * Function str2vec used to transform strings into vectors
 *
 ************************************************************************************************/

void str2vec(string& str, vector<char_vector>& vec)
{
	vec.clear();
	vector<char> single_vec;

	for(int i=0; i<str.length(); i++)
	{
		single_vec.clear();

		if(str[i]=='A')
		{
			single_vec.push_back('A');
		}
		else if(str[i]=='C')
		{
			single_vec.push_back('C');
		}
		else if(str[i]=='M')
		{
			single_vec.push_back('A');
			single_vec.push_back('C');
		}
		else if(str[i]=='G')
		{
			single_vec.push_back('G');
		}
		else if(str[i]=='R')
		{
			single_vec.push_back('A');
			single_vec.push_back('G');
		}
		else if(str[i]=='S')
		{
			single_vec.push_back('C');
			single_vec.push_back('G');
		}
		else if(str[i]=='V')
		{
			single_vec.push_back('A');
			single_vec.push_back('C');
			single_vec.push_back('G');
		}
		else if(str[i]=='T')
		{
			single_vec.push_back('T');
		}
		else if(str[i]=='W')
		{
			single_vec.push_back('A');
			single_vec.push_back('T');
		}
		else if(str[i]=='Y')
		{
			single_vec.push_back('C');
			single_vec.push_back('T');
		}
		else if(str[i]=='H')
		{
			single_vec.push_back('A');
			single_vec.push_back('C');
			single_vec.push_back('T');
		}
		else if(str[i]=='K')
		{
			single_vec.push_back('G');
			single_vec.push_back('T');
		}
		else if(str[i]=='D')
		{
			single_vec.push_back('A');
			single_vec.push_back('G');
			single_vec.push_back('T');
		}
		else if(str[i]=='B')
		{
			single_vec.push_back('G');
			single_vec.push_back('C');
			single_vec.push_back('T');
		}
		else
		{
			single_vec.push_back('A');
			single_vec.push_back('C');
			single_vec.push_back('G');
			single_vec.push_back('T');
		}

		vec.push_back(single_vec);
	}
}

/*************************************************************************************************
 *
 * Function eval_hd used to calculate similarity between motifs of the same length
 *
 ************************************************************************************************/

int eval_hd(string str1, string str2)
{
	if(str1.length()!=str2.length())
	{
		cerr<<"Error: The lengths of the two strings are not equal. Please check the details.\n";
		exit(1);
	}

	vector<char_vector> str_vec1, str_vec2;
	str2vec(str1, str_vec1);
	str2vec(str2, str_vec2);

	int dist=0;

	for(int i=0; i<str_vec1.size(); i++)
	{
		int dt=1;

		for(int j=0; j<str_vec1[i].size(); j++)
		{
			if(find(str_vec2[i].begin(), str_vec2[i].end(), str_vec1[i][j])!=str_vec2[i].end())
			{
				dt=0;
				break;
			}
		}

		dist+=dt;
	}

	return dist;
}

/*************************************************************************************************
 *
 * Function cmp_mtf used to determine the similarity between two motifs in slack conditions
 *
 ************************************************************************************************/

int cmp_mtf(mtf& mtf1, mtf& mtf2)
{
	int length;
	if(mtf1.alength >= mtf2.alength)
	{
		length=mtf2.alength;
	}
	else
	{
		length=mtf1.alength;
	}
	int min_val=length;
	int flag;
	int tmp_val;

	for(int i=0; i<=mtf1.alength-length; i++)
	{
		for(int j=0; j<=mtf2.alength-length; j++)
		{
			tmp_val=eval_hd(mtf1.deg.substr(i, length), mtf2.deg.substr(j, length));

			if(tmp_val<min_val)
			{
				min_val=tmp_val;
			}

			tmp_val=eval_hd(mtf1.deg.substr(i, length), mtf2.rev_deg.substr(j, length));

			if(tmp_val<min_val)
			{
				min_val=tmp_val;
			}
		}
	}

	if(min_val>=redundant_thr)
	{
		flag=0;
	}
	else
	{
		flag=1;
	}

	return flag;
}



/*************************************************************************************************
 *
 * Function sift_mtf used to delete the similar motifs with lower motif scores
 *
 ************************************************************************************************/

void sift_mtf(vector<mtf>& all_mtf)
{
	vector<int> flag;
	generate_flag(all_mtf.size(), flag);

	if(all_mtf.size() <= 0)
	{
		cerr << "0 motifs were identified!\n";
			exit(0);
	}

	for(int i=0; i<all_mtf.size()-1; i++)
	{
		for(int j=i+1; j<all_mtf.size(); j++)
		{
			if(flag[j]==0 && cmp_mtf(all_mtf[i], all_mtf[j])==1)
			{
				flag[j]=1;
			}
		}
	}

	vector<mtf> tmp_all;

	for(int i=0; i<all_mtf.size(); i++)
	{
		if(flag[i]==0)
		{
			tmp_all.push_back(all_mtf[i]);
		}
	}

	all_mtf=tmp_all;

	cout<<"Finished deleting similar motifs with lower motif scores.\n"<<"There are altogether "<<all_mtf.size()<<" motifs left."<<endl;
}

/*************************************************************************************************
 *
 * Function print_mtf used to print the result into one file
 *
 ************************************************************************************************/

void print_mtf(char *in_file, string& bg_file, string& out_file, vector<mtf>& all_mtf)
{
	string tmp_out=out_file;
	tmp_out+=".meme";
	
	ofstream f_out_op(tmp_out.c_str());

	if(!f_out_op)
	{
		cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
		exit(1);
	}

	f_out_op<<"# MEME 4.0.0\n"<<"# Command: ./ProSampler -i "<<in_file
		<<" -b "<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf
		<<" -f "<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "
		<<redundant_thr<<" -p "<<str_flag<<" -t "<<thr1<<" -c "<<hd_thr<<" -z "
		<<thr2<<" -h "<<help_flag<<endl;
	f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
	f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
	f_out_op<<"MEME Version 4\n"<<endl;
	f_out_op<<"ALPHABET= ACGT\n"<<endl;

	if(str_flag==1)
	{
		f_out_op<<"Strands: +\n"<<endl;
	}
	else
	{
		f_out_op<<"Strands: + -\n"<<endl;
	}

	f_out_op<<"Background letter frequencies (from dataset):\n";
	f_out_op.setf(ios::fixed);

	for(int i=0; i<3; i++)
	{
		f_out_op<<alphabet[i]<<" "<<fixed<<setprecision(PRECISE)<<nt_freq[i]<<" ";
	}
	f_out_op<<alphabet[3]<<" "<<fixed<<setprecision(PRECISE)<<nt_freq[3]<<endl<<endl;

	if(num_mtf>all_mtf.size())
	{
		num_mtf=all_mtf.size();
	}

	for(int i=0; i<num_mtf; i++)
	{
		f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "
			<<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
		f_out_op<<"letter-probability matrix: alength= 4 w= "
			<<all_mtf[i].alength<<" nsites= "<<all_mtf[i].nsites
			<<" score= "<<all_mtf[i].score<<endl;

		for(int j=0; j<all_mtf[i].alength; j++)
		{
			f_out_op.setf(ios::fixed);
			for(int k=0; k<3; k++)
			{
				f_out_op<<fixed<<setprecision(PRECISE)<<all_mtf[i].mat[k][j]<<" ";
			}
			f_out_op<<fixed<<setprecision(PRECISE)<<all_mtf[i].mat[3][j]<<endl;
		}

		f_out_op<<endl<<endl;
	}

	f_out_op<<endl;
	f_out_op.unsetf(ios::fixed);
	f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
	f_out_op.close();

	cout<<"Finished generating the motif file "<<tmp_out<<".\n";
	cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
}

/*************************************************************************************************
 *
 * Function print_site used to print the result into one file
 *
 ************************************************************************************************/

void print_site(char *in_file, string& bg_file, string& out_file, vector<mtf>& all_mtf)
{
	string tmp_out=out_file;
	tmp_out+=".site";
	ofstream f_out_op(tmp_out.c_str());

	if(!f_out_op)
	{
		cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
		exit(1);
	}

	f_out_op<<"# ProSampler 1.0.0\n"<<"# Command: ./ProSampler -i "<<in_file<<" -b "
		<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf<<" -f "
		<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "<<redundant_thr
		<<" -p "<<str_flag<<" -t "<<thr1<<" -c "<<hd_thr<<" -z "<<thr2<<" -h "<<help_flag<<endl;
	f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
	f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
	f_out_op<<"ProSampler Version 1.0.0\n"<<endl;
	f_out_op<<"ALPHABET= ACGT\n"<<endl;

	if(str_flag==1)
	{
		f_out_op<<"Strands: +\n"<<endl;
	}
	else
	{
		f_out_op<<"Strands: + -\n"<<endl;
	}

	f_out_op<<"Background letter frequencies (from dataset):\n";
	
	if(num_mtf>all_mtf.size())
	{
		num_mtf=all_mtf.size();
	}

	for(int i=0; i<num_mtf; i++)
	{
		f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "<<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
		f_out_op<<"letter-probability matrix: alength= 4 w= "<<all_mtf[i].alength<<" nsites= "<<all_mtf[i].nsites<<" score= "<<all_mtf[i].score<<endl<<endl<<"********************\n"<<endl;

		for(int j=0; j<all_mtf[i].site.size(); j++)
		{
			f_out_op<<all_mtf[i].site[j].header<<"\t"
				<<all_mtf[i].site[j].strand<<"\t"<<all_mtf[i].site[j].begin
				<<"\t"<<all_mtf[i].site[j].end<<endl;
		}

		f_out_op<<endl<<endl;
	}

	f_out_op<<endl;

	f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
	f_out_op.close();

	cout<<"Finished generating the motif file "<<tmp_out<<".\n";
	cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
	cout<<"Thank you for using ProSampler!\n";
}

/*************************************************************************************************
 *
 * Function get_spic used to get SPIC format output file from motif
 *
 ************************************************************************************************/

void get_spic(vector<mtf>& all_mtf, vector<spic>& all_spic)
{
	spic tmp_spic;

	for(int i=0; i<all_mtf.size(); i++)
	{
		tmp_spic.pfm=all_mtf[i].mat;
		tmp_spic.pssm=all_mtf[i].mat;
		
		for(int j=0; j<4; j++)
		{
			for(int k=0; k<all_mtf[i].alength; k++)
			{
				tmp_spic.pfm[j][k]*=all_mtf[i].nsites;

				if(tmp_spic.pssm[j][k]==0)
				{
					tmp_spic.pssm[j][k]=0.00001;
				}

				tmp_spic.pssm[j][k]/=nt_freq[j];
				tmp_spic.pssm[j][k]=log(tmp_spic.pssm[j][k]);
			}
		}

		tmp_spic.ic.clear();

		for(int j=0; j<all_mtf[i].alength; j++)
		{
			float tmp_sum=0;

			for(int k=0; k<4; k++)
			{
				tmp_sum+=all_mtf[i].mat[k][j]*tmp_spic.pssm[k][j];
			}

			tmp_spic.ic.push_back(tmp_sum);
		}

		all_spic.push_back(tmp_spic);
	}

	cout<<"Finished transforming motifs into SPIC format.\n";
}

/*************************************************************************************************
 *
 * Function print_spic used to print the result into one file
 *
 ************************************************************************************************/

void print_spic(char *in_file, string& bg_file, string& out_file, vector<spic>& all_spic, 
		vector<mtf>& all_mtf)
{
	string tmp_out=out_file;
	tmp_out+=".spic";
	ofstream f_out_op(tmp_out.c_str());

	if(!f_out_op)
	{
		cerr<<"Error: Can't open file "<<tmp_out.c_str()<<" for output!\n";
		exit(1);
	}

	f_out_op<<"# ProSampler 1.0.0\n"<<"# Command: ./ProSampler -i "<<in_file
		<<" -b "<<bg_file.c_str()<<" -o "<<out_file.c_str()<<" -d "<<num_deg<<" -m "<<num_mtf
		<<" -f "<<num_iter<<" -k "<<kmer_length<<" -l "<<lmer_length<<" -r "
		<<redundant_thr<<" -p "<<str_flag<<" -t "<<thr1<<" -c "<<hd_thr<<" -z "
		<<thr2<<" -h "<<help_flag<<endl;
	f_out_op<<endl<<"# Begin: "<<begin_pkg<<endl;
	f_out_op<<"#   End: "<<end_pkg<<endl<<endl;
	f_out_op<<"ProSampler Version 1.0.0\n"<<endl;
	f_out_op<<"ALPHABET= ACGT\n"<<endl;

	if(str_flag==1)
	{
		f_out_op<<"Strands: +\n"<<endl;
	}
	else
	{
		f_out_op<<"Strands: + -\n"<<endl;
	}

	f_out_op<<"Background letter frequencies (from dataset):\n";
	
	if(num_mtf>all_spic.size())
	{
		num_mtf=all_spic.size();
	}

	for(int i=0; i<num_mtf; i++)
	{
		f_out_op<<"MOTIF "<<all_mtf[i].deg<<" "<<all_mtf[i].rev_deg<<" ProSampler\n"<<endl;
		
		for(int j=0; j<4; j++)
		{
			f_out_op<<alphabet[j]<<"\t";
			
			for(int k=0; k<all_mtf[i].alength-1; k++)
			{
				f_out_op<<setiosflags(ios::fixed)<<setprecision(PRECISE);
				f_out_op<<setw(6)<<all_spic[i].pssm[j][k]<<"\t";
			}
			f_out_op<<setiosflags(ios::fixed)<<setprecision(PRECISE);
			f_out_op<<setw(6)<<all_spic[i].pssm[j][all_mtf[i].alength-1]<<endl;
		}
		
		f_out_op<<"I"<<"\t";
		f_out_op.setf(ios::fixed);

		for(int k=0; k<all_mtf[i].alength-1; k++)
		{
			f_out_op<<fixed<<setprecision(PRECISE);
			f_out_op<<setw(6)<<all_spic[i].ic[k]<<"\t";
		}
		f_out_op<<fixed<<setprecision(PRECISE);
		f_out_op<<setw(6)<<all_spic[i].ic[all_mtf[i].alength-1]<<endl;

		f_out_op.unsetf(ios::fixed);

		for(int j=0; j<4; j++)
		{
			char tmp_char;
			tmp_char=tolower(alphabet[j]);
			f_out_op<<tmp_char<<"\t";
			for(int k=0; k<all_mtf[i].alength-1; k++)
			{
				f_out_op<<int(all_spic[i].pfm[j][k])<<"\t";
			}
			f_out_op<<int(all_spic[i].pfm[j][all_mtf[i].alength-1])<<endl;
		}

		f_out_op<<endl;
	}

	f_out_op<<endl;

	f_out_op<<"Time "<<difftime(end_t, begin_t)<<" secs.\n";
	f_out_op.close();

	cout<<"Finished generating the motif file "<<tmp_out<<".\n";
	cout<<"There're altogether "<<num_mtf<<" motifs output.\n";
	cout<<"Thank you for using ProSampler!\n";
}

/*************************************************************************************************
 *
 * Main function
 *
 ************************************************************************************************/

 int main(int argc, char **argv)
 {
 	 begin_t=time(0);
	strftime(begin_pkg, sizeof(begin_pkg), 
			 "%b %d %Y %a %X %Z", localtime(&begin_t));

	parse_opt(argc, argv);

	if(help_flag==1)
	{
		usage();
	}

	r_alphabet['A']=0;
	r_alphabet['C']=1;
	r_alphabet['G']=2;
	r_alphabet['T']=3;

	sequence seq_nondeg, bg_nondeg;	// sequence with degenerate positions
	load_data(argv[f_in_id], seq_nondeg);
	int markov_flag;
	string bg_file = "order Markov chain model";
	string OutFile = argv[f_in_id];
	map_string kmer_B, kmer_T; // the 3-mer frequencies and transition matrix

	if(f_out_id != -1)
	{
		OutFile = argv[f_out_id];
	}

	if(f_bg_id == -1)
	{
		markov_flag = 3;
		bg_file = "3 order Markov chain model";
		markov(seq_nondeg, bg_nondeg, markov_flag);
	}
	else if(!isdigit(argv[f_bg_id][0]))
	{
		bg_file = argv[f_bg_id];
		load_data(argv[f_bg_id], bg_nondeg);
	}
	else if(atoi(argv[f_bg_id]) < 3)
	{
		bg_file = argv[f_bg_id] + ' ' + bg_file;
		markov_flag = atoi(argv[f_bg_id]);
		markov(seq_nondeg, bg_nondeg, markov_flag);
	}
	else
	{
		bg_file = "3 order Markov chain model";
		cout << "Using 3-rd order Markov Chain model as background.\n";
		markov_flag = 3;
		kmer_package kmer_pkg = markov(seq_nondeg, bg_nondeg, markov_flag);
		kmer_B = kmer_pkg.B; // the 3-mer frequencies
		kmer_T = kmer_pkg.T; // the transition matrix
	}

	de_lower(seq_nondeg);
	de_other(seq_nondeg);
	nt_stat(seq_nondeg, nt_freq);
	sequence seq_final;
	sequence bg_final;
	seq_generate(seq_nondeg, seq_final, str_flag);
	seq_generate(bg_nondeg, bg_final, str_flag);
	free_seq(seq_nondeg);
	free_seq(bg_nondeg);

	map<string, int_vector> kmer_site;
	kmer_count(kmer_site, seq_final, lmer_length);

	map<string, int_vector> kmer_bg;	// hold the k-mers in background sequences
	kmer_count(kmer_bg, bg_final, lmer_length);
	kmer_set major_set, minor_set;	// k-mer sets
	map<string, int> fa_flag;	
	// The flag to represent the status of k-mers in ChIP-seq data
	map<string, int> bg_flag;	// Status in background data
	for(map<string, int_vector>::iterator it=kmer_site.begin(); it!=kmer_site.end(); it++)
	{
		fa_flag[it->first]=1;
	}

	if(str_flag != 1)
	{
		comb_kmer(kmer_site, fa_flag, seq_final.strand, 
				kmer_bg, bg_flag, bg_final.strand);
	}
	choose_kmer(kmer_site, kmer_bg, major_set, 
			 seq_final.name.size(), fa_flag, minor_set); // select the significant k-mers
	kmer_sort(major_set);
	kmer_sort(minor_set);
	kmer_site.clear();
	kmer_bg.clear();
	map<string, int_vector>().swap(kmer_site);
	map<string, int_vector>().swap(kmer_bg);
	free_seq(bg_nondeg);
	free_seq(bg_final);
	top_num=combine_major(major_set, minor_set);
	free_kmer(minor_set);
	vector<kmer_all> kmer_final;
	fill_kmer(major_set, kmer_final, seq_final.name.size(), seq_final);
	free_kmer(major_set);
	vector<psm> psm_set;
	kmer2psm(kmer_final, psm_set, top_num);
	fill_psm(psm_set, kmer_final, kmer_length);
	int seq_num=seq_final.name.size();
 	cover_psm(psm_set, kmer_final);
    test_psm(psm_set, sum_exp, sum_bg);
	comp_cons(psm_set);
	psm_sort(psm_set);
	vector<int_vector> pre_out;
	int num_pre;
	
	if(num_mtf != -1)
	{
		 num_pre=MULTIPLE*num_mtf;
	}
	else
	{
		 num_pre=psm_set.size();
	}

	vector<int_vector> psm_nbr;
	build_graph(psm_set, psm_nbr);
	gibbs_sam(pre_out, psm_set, kmer_final,
			 num_pre, num_iter, psm_nbr, seq_num);
	psm_set.clear();
	vector<psm>().swap(psm_set);
	psm_nbr.clear();
	vector<int_vector>().swap(psm_nbr);
	vector<pre_mtf> all_pre;
	get_pre(all_pre, pre_out, kmer_final, seq_num);
	vector<kmer_all>().swap(kmer_final);
	pre_out.clear();
	vector<int_vector>().swap(pre_out);
	refine_mtf(all_pre);
	sift_pvalue(all_pre, kmer_final); // sift preliminary motifs using chi square p-values
	// whether to perform correction?
	kmer_final.clear();
	make_wild(wild_card, reverse_wild);
	vector<mtf> all_mtf;
	pre2mtf(all_pre, all_mtf, seq_final);
	all_pre.clear();
	vector<pre_mtf>().swap(all_pre);
	sort_mtf(all_mtf);
	sift_mtf(all_mtf);

	end_t=time(0);
	strftime(end_pkg, sizeof(end_pkg), 
			 "%b %d %Y %a %X %Z", localtime(&end_t));
	print_mtf(argv[f_in_id], bg_file, OutFile, all_mtf);
	print_site(argv[f_in_id], bg_file, OutFile, all_mtf);
	vector<spic> all_spic;
	all_spic.reserve(all_mtf.size());
	get_spic(all_mtf, all_spic);
	print_spic(argv[f_in_id], bg_file, OutFile, all_spic, all_mtf);

	return 0;
 }
