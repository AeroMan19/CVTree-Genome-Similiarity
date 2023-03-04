#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include <algorithm>

int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010
#define NUM_THREADS		4
void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);   number of possible 4-mers
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);	 number of possible 5-mers
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);       number of possible 6-mers
}

class Bacteria
{
private:
	std::unordered_map<long, double> vector;
	double* second;
	double one_l[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;

	void InitVectors()
	{
		second = new double[M1];
		memset(second, 0, M1 * sizeof(double));
		memset(one_l, 0, AA_NUMBER * sizeof(double));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i = 0; i < LEN - 1; i++)
		{
			short enc = encode(buffer[i]);
			one_l[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		second[indexs]++;
	}

	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		one_l[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		vector[index]++;
		total++;
		indexs = (indexs % M2) * AA_NUMBER + enc;
		second[indexs]++;
	}

public:
	long count;
	std::vector<std::pair<long, double>> tvVector; // sort can sort by first item in pair

	Bacteria(char* filename)
	{
		FILE* bacteria_file;
		errno_t OK = fopen_s(&bacteria_file, filename, "r");

		if (OK != 0)
		{
			fprintf(stderr, "Error: failed to open file %s\n", filename);
			exit(1);
		}

		InitVectors();

		char ch;
		while ((ch = fgetc(bacteria_file)) != EOF)
		{
			if (ch == '>')
			{
				while (fgetc(bacteria_file) != '\n'); // skip rest of line

				char buffer[LEN - 1];
				fread(buffer, sizeof(char), LEN - 1, bacteria_file);
				init_buffer(buffer);
			}
			else if (ch != '\n')
				cont_buffer(ch);
		}
		long total_plus_complement = total + complement;
		double total_div_2 = total * 0.5;
		int i_mod_aa_number = 0;
		int i_div_aa_number = 0;
		long i_mod_M1 = 0;
		long i_div_M1 = 0;

		for (int i = 0; i < AA_NUMBER; i++)
			one_l[i] = (double)one_l[i] / total_l;

		for (int i = 0; i < M1; i++)
			second[i] = (double)second[i] / total_plus_complement;

		count = 0;

		for (long i = 0; i < M; i++)
		{
			double p1 = second[i_div_aa_number];
			double p2 = one_l[i_mod_aa_number];
			double p3 = second[i_mod_M1];
			double p4 = one_l[i_div_M1];
			double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

			if (i_mod_aa_number == AA_NUMBER - 1)
			{
				i_mod_aa_number = 0;
				i_div_aa_number++;
			}
			else
				i_mod_aa_number++;

			if (i_mod_M1 == M1 - 1)
			{
				i_mod_M1 = 0;
				i_div_M1++;
			}
			else
				i_mod_M1++;

			if (stochastic > EPSILON)
			{
				double freq_i = 0.0;
				auto it = vector.find(i);
				if (it != vector.end()) freq_i = it->second;
				count++;
				tvVector.push_back({ i, (freq_i - stochastic) / stochastic });
			}
		}

		delete second;
		fclose(bacteria_file);
	}
};

// proccesses list.txt, and creates vector of file names 
void ReadInputFile(const char* input_name)
{
	FILE* input_file;
	errno_t OK = fopen_s(&input_file, input_name, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
		exit(1);
	}

	fscanf_s(input_file, "%d", &number_bacteria);
	bacteria_name = new char* [number_bacteria];

	for (long i = 0; i < number_bacteria; i++)
	{
		char name[10];
		fscanf_s(input_file, "%s", name, 10);
		bacteria_name[i] = new char[20];
		sprintf_s(bacteria_name[i], 20, "data/%s.faa", name);
	}
	fclose(input_file);
}

double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1 = 0;
	double vector_len2 = 0;
	long p1 = 0;
	long p2 = 0;
	while (p1 < b1->count && p2 < b2->count) // stop whenever a bacteria runs out of info
	{
		long n1 = b1->tvVector[p1].first;
		long n2 = b2->tvVector[p2].first;
		if (n1 < n2)
		{
			double t1 = b1->tvVector[p1].second;
			vector_len1 += (t1 * t1);
			p1++;
		}
		else if (n2 < n1)
		{
			double t2 = b2->tvVector[p2].second;
			p2++;
			vector_len2 += (t2 * t2);
		}
		else
		{
			double t1 = b1->tvVector[p1++].second;
			double t2 = b2->tvVector[p2++].second;
			vector_len1 += (t1 * t1);
			vector_len2 += (t2 * t2);
			correlation += t1 * t2;
		}
	}
	while (p1 < b1->count)
	{
		double t1 = b1->tvVector[p1++].second;
		vector_len1 += (t1 * t1);
	}
	while (p2 < b2->count)
	{
		double t2 = b2->tvVector[p2++].second;
		vector_len2 += (t2 * t2);
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void CompareAllBacteriaPar()
{
	Bacteria** b = new Bacteria * [number_bacteria];

	int i;
	#pragma omp parallel for schedule(dynamic) 
	for (i = 0; i < number_bacteria; i++)
	{
		b[i] = new Bacteria(bacteria_name[i]);
	}

	int j;
	//std::ofstream results("test.csv");
	#pragma omp parallel for private(j) schedule(dynamic) 
	for (i = 0; i < number_bacteria - 1; i++)
	{
		for (j = i + 1; j < number_bacteria; j++)
		{
			double ans = CompareBacteria(b[i], b[j]);
			#pragma omp critical
			{
				printf("%2d %2d -> %.20lf\n", i, j, ans);
				//results << i << "," << j << "," << ans << "\n";
			}

		}
	}
	//results.close();
}



void CompareAllBacteriaSeq()
{
	// read all bacteria info into buffer so we don't have to 
	// re-read it everytime a particular bacteria is required
	// Store the bacteria's info in a Bacteria Object.
	Bacteria** b = new Bacteria * [number_bacteria];

	for (int i = 0; i < number_bacteria; i++)
	{
		b[i] = new Bacteria(bacteria_name[i]);
	}

	for (int i = 0; i < number_bacteria - 1; i++)
		for (int j = i + 1; j < number_bacteria; j++)
		{
			double correlation = CompareBacteria(b[i], b[j]);
			printf("%2d %2d -> %.20lf\n", i, j, correlation);
			//results << i << "," << j << "," << correlation << "\n";
		}
	//results.close();
}


int main(int argc, char* argv[])
{
	// start timer
	auto start_time = std::chrono::high_resolution_clock::now();

	omp_set_dynamic(0);
	omp_set_num_threads(12);
	Init();
	ReadInputFile("list.txt");
	CompareAllBacteriaPar();

	// end timer
	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end_time - start_time;
	printf("time elapsed: %f seconds\n", duration.count());
	return 0;
}






