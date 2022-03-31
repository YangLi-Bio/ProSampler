/****************************************************************************
 *                                                                          *
 *                                 Head files                               *
 *                                                                          *
 ****************************************************************************/


#include<iostream> // standard input and output stream
#include<fstream> // file input and output stream
#include<string> // C++ string class
#include<string.h> // C string class
#include<stdlib.h> // a C class containing many basic functions for variable transformation
#include<stdint.h> // a C class containing many macros regarding integers
#include<vector> // vector class
#include<map> // hash table class
#include<sstream> // stream to segment strings when reading files

using namespace std;


/****************************************************************************
 *                                                                          *
 *                      Macros and complex classes                          *
 *                                                                          *
 ****************************************************************************/


#define SIZE 4 // the size of alphabet, i.e., ACGT

typedef struct
{
	vector<string> name; // introduction information of sequencecs
	vector<string> strand; // the ACGT sequences of a dataset
} sequence; // DNA sequence struct

typedef map<string, float> map_string; // a hash table to map strings to float numbers, i.e., k-mers to frequencies

struct kmer_package {
    map_string AA; // one base frequencies
    map_string B; // two base frequencies
    map_string C; // three base frequencies
    map_string D; // four base frequencies
}; // transition martrices used to define Markov chain model


/****************************************************************************
 *                                                                          *
 *                                 Functions                                *
 *                                                                          *
 ****************************************************************************/


// Generate a random float number in a range uniformly
float genunf(float low, float high)
{
	  float random = ((float) rand()) / (float) RAND_MAX; // generate a random number without bounds
    float diff = high - low; // set the range of random numbers
    float r = random * diff; // zoom in the random number
    return low + r; // added to the lower bound
}


// Partition a vector of numbers into two sets without random sizes
// after partition, numbers in the first set are smaller than the second set
// this is a divide-and-conpueror approach for quickSort
template <class T> // T can be int, float, or double
int partition(T a[], int start, int stop, int id[])
{
        int temp_id, up = start, down = stop - 1;
        T temp_value, part = a[stop];
        if(stop <= start) return start;

        while (true)
        {    
                while(a[up] < part) up ++;
                while(part < a[down] && (up < down)) down --;

                if(up >= down) break;

                temp_value = a[up];  a[up] = a[down]; a[down] = temp_value;
                temp_id = id[up]; id[up] = id[down]; id[down] = temp_id;

                up ++; down --;
        }    

        temp_value = a[up]; a[up] = a[stop]; a[stop] = temp_value;
        temp_id = id[up]; id[up] = id[stop]; id[stop] = temp_id;
        return up;
}


// This is the fastest algorithm for sorting in implementation
template <class T>
void quickSort(T a[], int start, int stop, int id[])
{
        int i;
        if (stop <= start) return;

        i = partition(a, start, stop, id);
        quickSort(a,start,i - 1, id);
        quickSort(a, i + 1, stop, id);
}


// Generate a random number according to a distribution saved in an array
int genmulone(float *p_input, long ncat)
{
	int i;
	float *p = NULL;
	p = new float [ncat];
	for (i = 0; i < ncat; i ++) p[i] = p_input[i];
	
	int*q = NULL;
	q = new int [ncat];
	for (i = 0; i < ncat;i ++) q[i] = i;
	
	quickSort(p, 0, ncat - 1, q);
	
	float prob = genunf(0, 1);
	int outcome = ncat - 1;
	while (outcome >= 0  && prob > p[outcome])
	{
		prob -= p[outcome];
		outcome --;
	}
	
	return q[outcome];
}


// Genrate Markov Chain model
// input the experiment sequences, file path to save the background sequences, and the order 
// of Markov Chain model
// output the Markov Chain transition matrices saved in kmer_package object
struct kmer_package markov(sequence& f_in)
{
  // Construct the hash table for various k-mers (k = 1, 2, 3, 4)
	char* alphabet;	// the alphabet for DNA nucleotides A, C, G, T
	alphabet = new char [SIZE];

	alphabet[0]='A';
	alphabet[1]='C';
	alphabet[2]='G';
	alphabet[3]='T';

	map_string hash_onebase, hash_twobase, hash_threebase, hash_fourbase;
	// hash tables to correspond k-mers (k = 1, 2, 3, and 4) to their frequencies

	
	// For one-kmer
	string onebase;	// the string to contain A, C, G, T
	
	// Initialization of the hash table by zeros
	for (int i = 0; i < SIZE; i ++)
	{
	  onebase.clear();
	  onebase.push_back(alphabet[i]);
	  hash_onebase[onebase] = 0;
	}	
	
	
	// For two-mer
	string twobase;
	
	// initialization for the hash table of 2-mers
	for (int i = 0; i < SIZE; i ++)
	{
	  for (int j = 0; j < SIZE; j ++)
	  {
	    twobase.clear();
	    twobase.push_back(alphabet[i]);
	    twobase.push_back(alphabet[j]);
	    hash_twobase[twobase] = 0;
	  }
	}
	
	
	// For three-mer
	string threebase;
	
	// initialization for the hash table of 3-mers
	for (int i = 0; i < SIZE; i ++)
	{
	  for (int j = 0; j < SIZE; j ++)
	  {
	    for (int k = 0; k < SIZE; k ++)
	    {
	      threebase.clear();
	      threebase.push_back(alphabet[i]);
	      threebase.push_back(alphabet[j]);
	      threebase.push_back(alphabet[k]);
	      hash_threebase[threebase] = 0;
	    }
	  }
	}
	
	
	// For four-mer
	string fourbase;
	
	// Initialization for the hash table of 4-mers
	for (int i = 0; i < SIZE; i ++)
	{
	  for (int j = 0; j < SIZE; j ++)
	  {
	    for (int k = 0; k < SIZE; k ++)
	    {
	      for (int l = 0; l < SIZE; l ++)
	      {
	        fourbase.clear();
	        fourbase.push_back(alphabet[i]);
	        fourbase.push_back(alphabet[j]);
	        fourbase.push_back(alphabet[k]);
	        fourbase.push_back(alphabet[l]);
	        hash_fourbase[fourbase] = 0;
	      }
	    }
	  }
	}
	

	// Count the frequencies of k-mers (k = 0, 1, 2, and 3)
	
	for (int i = 0; i < f_in.strand.size(); i ++)
	{
	  for (int j = 0; j < f_in.strand[i].size(); j ++)
	  {
	    if (f_in.strand[i][j] != 'N') // not 'N'
	    {
	      hash_onebase[f_in.strand[i].substr(j, 1)] ++; // 1-mers
	    }
	    if (j < f_in.strand[i].size() - 1 && f_in.strand[i].substr(j, 2).find('N') 
           == f_in.strand[i].substr(j, 2).npos) // no 'N' is included
	    {
	      hash_twobase[f_in.strand[i].substr(j, 2)] ++; // 2-mers
	    }
	    
	    if (j < f_in.strand[i].size() - 2 && f_in.strand[i].substr(j, 3).find('N') 
           == f_in.strand[i].substr(j, 3).npos) // no 'N' is included
	    {
	      hash_threebase[f_in.strand[i].substr(j,3)] ++; // 3-mers
	    }
	    
	    if (j < f_in.strand[i].size() - 3 && f_in.strand[i].substr(j, 4).find('N') 
           == f_in.strand[i].substr(j, 4).npos) // no 'N' is included
	    {
	      hash_fourbase[f_in.strand[i].substr(j, 4)] ++; // 4-mers
	    }
	  }
	}

	cout << "k-mer counting has been completed.\n";


	// Build the transition matrices
	map<string, float>::iterator iter; // an iterator to traverse a hash table
	float sum = 0;	// the sum of all nucleotides in sequences

	map_string threebase_freq = hash_threebase; // make a copy of the hash table of 3-mers
	
	
	// Normalize the 4-mers by dividing their prefix 3-mers
	for (iter = hash_fourbase.begin(); iter != hash_fourbase.end(); iter ++)
	{
	  iter->second /= hash_threebase[iter->first.substr(0, 3)]; // divide the occurrence of 
	  // a 4-mer by the that of its 3-prefix
	}
	
	
	// Normalize the 3-mers by dividing their prefix 2-mers
	float threebase_sum = 0; // the sum of frequencies of all 3-mers
	
	for (iter = hash_threebase.begin(); iter != hash_threebase.end(); iter ++)
	{
	  iter->second /= hash_twobase[iter->first.substr(0, 2)];
	  threebase_sum += iter->second; // add the frequency of each 3-mer
	}
	
	for (iter = threebase_freq.begin(); iter != threebase_freq.end(); iter ++)
	{
	  iter->second /= threebase_sum;
	} // divide the occurrence of 3-mers by the sum of occurrences of all 3-mers
	
	
	// Normalize the 2-mers by dividing their prefix 1-mers
	for (iter = hash_twobase.begin(); iter != hash_twobase.end(); iter ++)
	{
	  iter->second /= hash_onebase[iter->first.substr(0, 1)];
	}
	
	
	// Normalize the 1-mers by dividing the total number of nucleotides
	for (iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
	{
	  sum += iter->second;// sum of frequencies of all nucleotides
	}
	
	for (iter = hash_onebase.begin(); iter != hash_onebase.end(); iter ++)
	{
	  iter->second /= sum; // divide the frequencies of all nucleotides by
	  // the sum of frequencies of all nucleotides
	}
	

	cout << "Finished claculating the k-mer frequencies.\n";

	
	struct kmer_package kmer_pkg; // the 3-mer frequencies and transition matrix
	kmer_pkg.AA = hash_onebase; 
	kmer_pkg.B = hash_twobase;
	kmer_pkg.C = hash_threebase; // the 3-mer frequencies
	kmer_pkg.D = hash_fourbase; // the hash table to count the 4-mer frequencies,
	// which can serve as the transition matrix B
	
	return(kmer_pkg);
}


// Count the k-mers in background based on the probability model
void kmer_in_bg(map<string, float>& kmer_bg, struct kmer_package kmer_pkg, 
                map<string, vector<int> > kmer_site, int nkmer)
{
  // kmer_bg : the struct to save the k-mer frequencies
  // kmer_pkg : the probability model of background sequences
  
  // int nkmer = 0;
  // 
  // for (map<string, int_vector>::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
  // {
  //   nkmer += (it -> second).size(); // Add the occurrence of this k-mer
  // }
  
  cout << "There are in total " << nkmer << " " << (kmer_site.begin() -> first).length() << 
    "-mers.\n";
  
  for (map<string, vector<int> >::iterator it = kmer_site.begin(); it != kmer_site.end(); it ++)
  {
    string str = it -> first; // get the string from positive sequences
    kmer_bg[str] = nkmer * kmer_pkg.AA[str.substr(0, 1)] * kmer_pkg.B[str.substr(0, 2)] * 
      kmer_pkg.C[str.substr(0, 3)]; // The first three nucleotides
    for (int i = 3; i < str.length(); i ++)
    {
      // cout << kmer_pkg.D[str.substr(i - 3, 4)] << "\n";
      kmer_bg[str] *= kmer_pkg.D[str.substr(i - 3, 4)];
    }
    // cout << kmer_bg[it -> first] << "\n";
    // exit(0);
    // kmer_bg[it -> first] *= nkmer;
    // cout << str << " " << kmer_bg[str] << "\n";
    // exit(0);
  }
  
  // cout << kmer_bg.size() << "\n";
  // exit(0);
  
  cout << "Finish calculating the frequencies for " << nkmer << " " << (kmer_site.begin() -> first).length() << 
    "-mers.\n";
}