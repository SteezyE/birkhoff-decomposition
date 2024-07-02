#include <bits/stdc++.h>

using namespace std;
using edge = pair<int,int>;

void output_mtx(vector<edge> &edges, vector<double> &values, int n, string filename)
{
	ofstream file;
	file.open(filename+"_skewed.mtx");
	file << "%%MatrixMarket matrix coordinate real general\n";
	file << n << " " << n << " " << edges.size() << '\n';
	int k=0;
	for(const auto &[i, j] : edges)
		file << i+1 << " " << j+1 << " " << values[k++] << "\n";
	file.close();
}

void generate_mtx(int n, int pb, int ps, double s)
{
	assert(n > 0 && pb > 0 && ps > 0 && s >= 0.0 && s <= 1.0);
	int p = pb + ps;
	double factors[p];
	int * permutations = (int *) malloc(p * n * sizeof(int));
	double bf = s / pb;
	double bs = (1.0-s) / ps;
	for(int j=0; j<p; ++j) 
	{
		factors[j] = j < pb ? bf : bs;
		for(int k=0; k<n; ++k)
			permutations[j*n+k] = k; 
	}
	// Fisher-Yates shuffle
	for(int j=0; j<p; ++j)
	{
		for(int k=0; k<n-1; ++k)
		{
			int h = k + rand() % (n - k);
			swap(permutations[j*n+k], permutations[j*n+h]);	
		}
	}
	int offset = 0;
	int * positions = (int *) calloc(n, sizeof(int));
	int * r_positions = (int *) calloc(n, sizeof(int));
	bool * column_indices = (bool *) calloc(n, sizeof(bool));
	vector<edge> edges;
	vector<double> values;
	for(int k=0; k<n; ++k)
	{
		int h = 0;
		offset = values.size();
		for(int j=0; j<p; ++j)
		{
			int index = permutations[j*n+k];
			if(!column_indices[index])
			{
				positions[index] = h;
				r_positions[h++] = index;
				column_indices[index] = 1;
				values.push_back(factors[j]);
				edges.push_back({k, index});
			}
			else
				values[offset + positions[index]] += factors[j];
		}
		while(h--)
			column_indices[r_positions[h]] = 0;
	}
	output_mtx(edges, values, n, to_string(n)+"_"+
                                     to_string(pb)+"_"+
                                     to_string(ps)+"_"+
                                     to_string(int(s*100)));
	free(positions); free(r_positions); free(column_indices); free(permutations);
}

int main(int argc, char **argv)
{
	if(argc != 5)	
	{
		printf("USAGE: mm_generator_skewed [N] [B] [S] [V]\n\n"
                       "N := matrix size\n"
                       "B := big flow count\n"
                       "S := small flow count\n"
                       "V := skew value\n");
		return 1;
	}
	generate_mtx(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atof(argv[4]));
	return 0;
}
