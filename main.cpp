#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include <thread>
#include <mutex>
#include <set>
#include <memory>
#include <limits>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <stdlib.h>

#include<omp.h>

#include <seal/seal.h>

using namespace std;
using namespace seal;

int Find(vector<int> &root, int index) {
	return root[index] == index ? index : root[index] = Find(root, root[index]);
}

void Union(vector<int> &root, int lhs, int rhs) {
	int fl = Find(root, lhs);
	int fr = Find(root, rhs);
	if (fl != fr) {
		root[fl] = fr;
	}
}

template<class T>
vector<double> standardization(const vector<T> &input) {
	double avg = 0;
	for (int i = 0; i < input.size(); i++) {
		avg += input[i];
	}
	avg /= input.size();
	double s = 0;
	for (int i = 0; i < input.size(); i++) {
		s += (input[i] - avg) * (input[i] - avg);
	}
	s /= input.size();
	s = sqrt(s);
	vector<double> result(input.size());
	for (int i = 0; i < input.size(); i++) {
		result[i] = (input[i] - avg) / s;
	}
	return result;
}

template<class T>
vector<double> standardization_sq(const vector<T> &input) {
	double avg = 0;
	for (int i = 0; i < input.size(); i++) {
		avg += input[i];
	}
	avg /= input.size();
	double s = 0;
	for (int i = 0; i < input.size(); i++) {
		s += (input[i] - avg) * (input[i] - avg);
	}
	s /= input.size();
	s = sqrt(s);
	vector<double> result(input.size());
	for (int i = 0; i < input.size(); i++) {
		result[i] = input[i] / s;
	}
	return result;
}


template<class T>
vector<double> standardization_minmax(const vector<T> &input) {
	T ele_min = *min_element(input.begin(), input.end());
	T ele_max = *max_element(input.begin(), input.end());
	vector<double> result(input.size());
	for (int i = 0; i < input.size(); i++) {
		result[i] = 1.0 * (input[i] - ele_min) / (ele_max - ele_min);
	}
	return result;
}
/*
Helper function: Prints the name of the example in a fancy banner.
*/
void print_example_banner(string title)
{
	if (!title.empty())
	{
		size_t title_length = title.length();
		size_t banner_length = title_length + 2 + 2 * 10;
		string banner_top(banner_length, '*');
		string banner_middle = string(10, '*') + " " + title + " " + string(10, '*');

		cout << endl
			<< banner_top << endl
			<< banner_middle << endl
			<< banner_top << endl
			<< endl;
	}
}

/*
Helper function: Prints the parameters in a SEALContext.
*/
void print_parameters(shared_ptr<SEALContext> context)
{
	// Verify parameters
	if (!context)
	{
		throw invalid_argument("context is not set");
	}
	auto &context_data = *context->context_data();

	/*
	Which scheme are we using?
	*/
	string scheme_name;
	switch (context_data.parms().scheme())
	{
	case scheme_type::BFV:
		scheme_name = "BFV";
		break;
	case scheme_type::CKKS:
		scheme_name = "CKKS";
		break;
	default:
		throw invalid_argument("unsupported scheme");
	}

	cout << "/ Encryption parameters:" << endl;
	cout << "| scheme: " << scheme_name << endl;
	cout << "| poly_modulus_degree: " <<
		context_data.parms().poly_modulus_degree() << endl;

	/*
	Print the size of the true (product) coefficient modulus.
	*/
	cout << "| coeff_modulus size: " << context_data.
		total_coeff_modulus_bit_count() << " bits" << endl;

	/*
	For the BFV scheme print the plain_modulus parameter.
	*/
	if (context_data.parms().scheme() == scheme_type::BFV)
	{
		cout << "| plain_modulus: " << context_data.
			parms().plain_modulus().value() << endl;
	}

	cout << "\\ noise_standard_deviation: " << context_data.
		parms().noise_standard_deviation() << endl;
	cout << endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
ostream &operator <<(ostream &stream, parms_id_type parms_id)
{
	stream << hex << parms_id[0] << " " << parms_id[1] << " "
		<< parms_id[2] << " " << parms_id[3] << dec;
	return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template<typename T>
void print_vector(vector<T> vec, size_t print_size = 4, int prec = 3)
{
	/*
	Save the formatting information for std::cout.
	*/
	ios old_fmt(nullptr);
	old_fmt.copyfmt(cout);

	size_t slot_count = vec.size();

	cout << fixed << setprecision(prec) << endl;
	if (slot_count <= 2 * print_size)
	{
		cout << "    [";
		for (size_t i = 0; i < slot_count; i++)
		{
			cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
		}
	}
	else
	{
		vec.resize(max(vec.size(), 2 * print_size));
		cout << "    [";
		for (size_t i = 0; i < print_size; i++)
		{
			cout << " " << vec[i] << ",";
		}
		if (vec.size() > 2 * print_size)
		{
			cout << " ...,";
		}
		for (size_t i = slot_count - print_size; i < slot_count; i++)
		{
			cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
		}
	}
	cout << endl;

	/*
	Restore the old std::cout formatting.
	*/
	cout.copyfmt(old_fmt);
}

struct Data {
	int index;
	int type;
	vector<int64_t> data;

	int64_t& operator[] (int index) {
		return data[index];
	}
	size_t size() {
		return data.size();
	}
};

Data parser(const vector<int64_t> &vec) {
	if (vec.size() < 2) {
		return Data();
	}
	Data result;
	result.index = vec[0];
	result.type = vec[1];
	for (int i = 2; i < vec.size(); i++) {
		result.data.push_back(vec[i]);
	}
	return result;
}

vector<Data> read_data(const string &filename) {
	ifstream fin;
	fin.open(filename);
	if (!fin.is_open()) {
		cout << "Could not open this file!!!" << endl;
		exit(EXIT_FAILURE);
	}
	char buffer[1024];
	
	vector<Data> result;
	while (fin.getline(buffer, 1024)) {
		stringstream ss(buffer);
		vector<int64_t> vec;
		int64_t temp;
		while (ss >> temp) {
			vec.push_back(temp);
		}
		if (vec.size() < 1) {
			break;
		}
		result.push_back(parser(vec));
		if (result.size() > 1 && result[0].size() != result[result.size() - 1].size()) {
			cout << "Data Length Inconsistent in index:" << result[result.size() - 1].index << endl;
			exit(EXIT_FAILURE);
		}
	}
	fin.close();
	return result;
}

void norm(const vector<int64_t> &input) {

}

int main() {
#ifdef SEAL_VERSION
	cout << "SEAL version: " << SEAL_VERSION << endl;
#endif

	int cluster_num = 2;
	// Step 0
	// Init SEAL
	EncryptionParameters parms(scheme_type::BFV);

	parms.set_poly_modulus_degree(64);
	parms.set_coeff_modulus(coeff_modulus_128(8192));

	parms.set_plain_modulus(9100000000032769);

	auto context = SEALContext::Create(parms);
	print_parameters(context);

	auto qualifiers = context->context_data()->qualifiers();
	cout << "Batching enabled: " << boolalpha << qualifiers.using_batching << endl;

	KeyGenerator keygen(context);
	auto public_key = keygen.public_key();
	auto secret_key = keygen.secret_key();

	auto gal_keys = keygen.galois_keys(30);

	auto relin_keys = keygen.relin_keys(30);

	Encryptor encryptor(context, public_key);
	Evaluator evaluator(context);
	Decryptor decryptor(context, secret_key);

	BatchEncoder batch_encoder(context);

	size_t slot_count = batch_encoder.slot_count();
	size_t row_size = slot_count / 2;
	cout << "Plaintext matrix row size: " << row_size << endl;


	// Step 0
	// Load Dataset
	string input;
	vector<Data> dataset;
	dataset = read_data("data.txt");

	const int M = dataset.size();
	const int len = dataset[0].size();

	// Helper function to flatten
	auto pos = [M](int row, int col) {
		return ((2 * M - 1 - row) * row / 2 + col - row - 1);
	};
	const int flatten_size = (M + 1) * (M - 2) / 2 + 1;

	cout << "The training matrix has " << M << " rows and " << dataset[0].size() << " columns." << endl;



	// Step 0
	// Standardize
	cout << "Need standardize ? (Y/N)" << endl;
	cin >> input;
	if (toupper(input[0]) == 'Y') {
		const int scale = 1000;
		for (int index = 0; index < dataset[0].size(); index++) {
			double sum = 0, sqsum = 0;
			for (int i = 0; i < M; i++) {
				sum += dataset[i][index];
			}
			double avg = sum / M;
			for (int i = 0; i < M; i++) {
				sqsum += (dataset[i][index] - avg) * (dataset[i][index] - avg);
			}
			double s = sqrt(sqsum / M);
			for (int i = 0; i < M; i++) {
				dataset[i][index] = (dataset[i][index] - avg) / s;
			}
		}
	}

	cout << "Clustering algorithm:" << endl;
	cout << "1.Plaintext  clustering algorithm" << endl;
	cout << "2.Ciphertext clustering algorithm" << endl;
	cout << "Your choose: ";
	cin >> input;


	chrono::high_resolution_clock::time_point time_start, time_end;

	chrono::microseconds time_distance_calculation_sum(0);
	chrono::microseconds time_sort_sum(0);


	// Step 1
	// Calculate the distance Dij

	vector<vector<int64_t>> dist(M, vector<int64_t>(M, 0));
	vector<int64_t> dist_flatten(flatten_size, 0);

	time_start = chrono::high_resolution_clock::now();
	if (input[0] == '1') {
		// Plaintext
		for (int i = 0; i < M - 1; i++) {
			for (int j = i + 1; j < M; j++) {
				int64_t sum = 0;
				for (int k = 0; k < len; k++) {
					int64_t diff = dataset[i][k] - dataset[j][k];
					sum += diff * diff;
				}
				dist[i][j] = dist[j][i] = dist_flatten[pos(i,j)] = sqrt(sum);
			}
		}
	}
	else if (input[0] == '2') {
		// Cipher
		vector<Ciphertext> cipher_data(M);
		vector<Ciphertext> cipher_dist(flatten_size);
		for (int i = 0; i < M; i++) {
			Plaintext plain;
			batch_encoder.encode(dataset[i].data, plain);
			encryptor.encrypt(plain, cipher_data[i]);
		}
		
		omp_set_nested(1);//设置支持嵌套并行
#pragma omp parallel for
		for (int i = 0; i < M - 1; i++) {
#pragma omp parallel for
			for (int j = i + 1; j < M; j++) {
				Ciphertext rotate_temp;
				Ciphertext cipher_sum;

				evaluator.sub(cipher_data[i], cipher_data[j], cipher_sum);
				evaluator.square_inplace(cipher_sum);
				evaluator.relinearize_inplace(cipher_sum, relin_keys);

				rotate_temp = cipher_sum;
				
				// original rotate sum
				//for (int k = 1; k < dataset[i].size(); k++) {
				//	evaluator.rotate_rows_inplace(rotate_temp, 1, gal_keys);
				//	evaluator.add_inplace(cipher_sum, rotate_temp);
				//}

				// logn rotate sum
				for (int k = 1; k < dataset[i].size(); k <<= 1) {
					evaluator.rotate_rows(cipher_sum, k, gal_keys, rotate_temp);
					evaluator.add_inplace(cipher_sum, rotate_temp);
				}
				cipher_dist[pos(i, j)] = cipher_sum; 
			}
		}

		vector<vector<int64_t>> result_dist(flatten_size, vector<int64_t>(slot_count, 0));
		for (int i = 0; i < flatten_size; i++) {
			Plaintext plain;
			decryptor.decrypt(cipher_dist[i], plain);
			batch_encoder.decode(plain, result_dist[i]);
		}
		ofstream decrypted_distance("dec_distances.txt");
		for (int i = 0; i < M - 1; i++)
		{
			for (int j = i + 1; j < M; j++)
			{
				dist[i][j] = dist[j][i] = dist_flatten[pos(i, j)] = sqrt(result_dist[pos(i, j)][0]);
				
				decrypted_distance << i << "-" << j << ":" << dist_flatten[pos(i, j)] << endl;
			}
		}
		decrypted_distance.close();
	}

	cout << "The distance between all nodes is calculated!" << endl;
	time_end = chrono::high_resolution_clock::now();
	time_distance_calculation_sum += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	cout << "Distance calculation time: " << time_distance_calculation_sum.count() << " microseconds" << endl;



	// Step 2
	// Ascending the distance between all nodes. Take a cutoff distance dc.
	time_start = chrono::high_resolution_clock::now();
	sort(dist_flatten.begin(), dist_flatten.end());
	time_end = chrono::high_resolution_clock::now();
	time_sort_sum += chrono::duration_cast<chrono::microseconds>(time_end - time_start);
	cout << "Sorting time: " << time_sort_sum.count() << " microseconds" << endl;

	ofstream ascending_distance("ascending_distance.txt");
	for (int i = 0; i < flatten_size; i++)
	{
		ascending_distance << dist_flatten[i] << endl;
	}
	ascending_distance.close();

	double percent = 0.02;
	int percent_index = flatten_size * percent;
	int64_t dc = dist_flatten[percent_index];
	cout << "Percent_num: " << percent_index << endl;
	cout << "Cut-off distance: " << dc << endl;

	// Step 3
	// Calculate the density ρi of each node i
	vector<int> density(M);
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			if (dist[i][j] <= dc) {
				density[i]++;
			}
		}
	}

	ofstream Density("density.txt");
	for (int i = 0; i < M; i++)
	{
		Density << density[i] << endl;
	}
	Density.close();



	// Step 4
	// For each node i, find all nodes j that are denser than the node i, and select the smallest dij, denoted as δi.
	// In particular, for the node with the highest density (ρ), δi is the maximum distance from all nodes to node i


	// Older Density count
	//vector<int64_t> si(M, 0);
	//vector<int> max_count(M, 0);
	//for (int i = 0; i < M; i++) {
	//	for (int j = 0; j < M; j++) {
	//		if (density[j] > density[i]) {
	//			if (dist[i][j] < si[i]) {
	//				si[i] = dist[i][j];
	//			}
	//			max_count[i]++;
	//		}
	//	}
	//}

	vector<int> density_index(M);
	for (int i = 0; i < M; i++) {
		density_index[i] = i;
	}
	auto density_comp = [&density] (const int lhs, const int rhs) {
		return density[lhs] > density[rhs];
	};
	sort(density_index.begin(), density_index.end(), density_comp);

	vector<int64_t> si(M, 0);
#pragma omp parallel for
	for (int i = 1; i < M; i++) {
		si[density_index[i]] = dist[density_index[0]][density_index[i]];
		if (density[density_index[i]] == density[density_index[0]])
			continue;
		for (int j = 0; j < i && density[density_index[j]] > density[density_index[i]]; j++) {
			si[density_index[i]] = min(si[density_index[i]],
				dist[density_index[j]][density_index[i]]);
		}
		
		//for (int j = i + 1; j < M && density[density_index[i]] == density[density_index[j]]; j++) {
		//	si[density_index[i]] = min(si[density_index[i]],
		//		dist[density_index[j]][density_index[i]]);
		//}
	}
	si[density_index[0]] = 0;
	for (int i = 1; i < M; i++) {
		si[density_index[0]] = max(si[density_index[0]],
			dist[density_index[0]][density_index[i]]);
	}
	

	ofstream Si("si.txt");
	for (int i = 0; i < M; i++) {
		Si << si[i] << endl;
	}
	Si.close();

	// Step 5
	// Identify cluster centers
	vector<int64_t> mul(M);
	vector<double> density_std = standardization_sq(density);
	vector<double> si_std = standardization_sq(si);
	vector<int> mul_index(M);
	for (int i = 0; i < M; i++) {
		mul[i] = density[i] * si[i];
		mul_index[i] = i;
	}
	auto mul_comp = [&mul](const int lhs, const int rhs) {
		return mul[lhs] > mul[rhs];
	};
	sort(mul_index.begin(), mul_index.end(), mul_comp);

	set<int> cluster_center;
	for (int i = 0; i < cluster_num; i++) {
		cluster_center.insert(mul_index[i]);
		cout << "The " << i + 1 << " cluster center!" << endl;
		cout << "Index = " << mul_index[i]
			<< " , label = " << dataset[mul_index[i]].type
			<< " , mul = " << mul[mul_index[i]]
			<< " , density = " << density[mul_index[i]]
			<< " , si = " << si[mul_index[i]]
			<< " , den_std = " << density_std[mul_index[i]]
			<< " , si_std = " << si_std[mul_index[i]] << endl;
	}

	// Step 6
	// Allocate the remaining points (non-central points)
	
	// Union-Set Data structure
	vector<int> fa(M);
	for (int i = 0; i < M; i++) {
		fa[i] = i;
	}
	for (int i = 0; i < M; i++) {
		if (cluster_center.find(i) != cluster_center.end()) {
			// i is cluster center
			continue;
		}
		bool pre = false;
		for (auto &selected : cluster_center) {
			if (density[i] == density[selected] && dist[i][selected] <= dc) {
				Union(fa, i, selected);
				pre = true;
				break;
			}
		}
		if (pre) {
			continue;
		}
		int min_index = -1;
		for (int j = 0; j < M; j++) {
			if (density[j] > density[i]) {
				if (min_index == -1 || dist[i][min_index] > dist[i][j]) {
					min_index = j;
				}
			}
		}
		if (min_index != -1) {
			Union(fa, i, min_index);
		}
		else {
			cout << "fail" << endl;
		}
	}

	ofstream Label("label.txt");
	map<int, int> map_center;
	vector<int> label(M, -1);
	for (auto &center : cluster_center) {
		map_center[Find(fa, center)] = center;
	}
	for (int i = 0; i < M; i++)
	{
		label[i] = map_center[Find(fa, i)];
		Label << i << "-" << label[i] << endl;
		//Label << i << " -fa- " << Find(fa, i) << endl;
	}
	Label.close();


	// Step 7
	// Boundary and isolated point judgement


	// Step 8
	// Count each node tag and calculate clustering accuracy
	int node_labeled = M;
	int acc_count = 0;
	for (int i = 0; i < M; i++) {
		if (label[i] == -1) {
			node_labeled--;
		}
		else if (dataset[label[i]].type == dataset[i].type) {
			acc_count++;
		}
	}
	double acc = 1.0 * acc_count / node_labeled;
	cout << "acc = " << acc << endl;


	//cout << "Debug:" << endl;
	//map<int, int> type_count;
	//for (auto &ele : dataset) {
	//	type_count[ele.type]++;
	//}
	//for (auto &type : type_count) {
	//	cout << "Type " << type.first << " has " << type.second << " element" << endl;
	//}
	//for (int i = 0; i < M; i++) {
	//	if (i == 0 || density[density_index[i]] == density[density_index[0]]) {
	//		cout << "Data " << density_index[i] << " type: " << dataset[density_index[i]].type
	//			<< "  density = " << density[density_index[i]] << endl;
	//	}
	//}


	// S_Dbw

	// map the index of root node
	//map<int, vector<double>> c_s(cluster_num + 1, vector<double>(dataset[0].size()));
	//map<int, vector<double>> c_avg(cluster_num + 1, vector<double>(dataset[0].size()));
	//map<int, double> c_std(cluster_num + 1);
	//map<int, int> c_count(cluster_num + 1);
	map<int, vector<double>> c_s, c_avg;
	map<int, double> c_std;
	map<int, int> c_count;
	double stdev = 0;
	for (int i = 0; i < M; i++) {
		if (label[i] < 0) {
			continue;
		}
		if (c_avg.count(label[i]) == 0) {
			c_avg[label[i]].resize(dataset[0].size());
		}
		for (int k = 0; k < dataset[0].size(); k++) {
			c_avg[label[i]][k] += dataset[i][k];
		}
		c_count[label[i]]++;
	}
	//for (int i = 0; i <= cluster_num; i++) {
	for (int i : cluster_center) {
		if (c_count[i] > 0) {
			for (int k = 0; k < dataset[0].size(); k++) {
				c_avg[i][k] /= c_count[i];
			}
		}
	}

	for (int i = 0; i < M; i++) {
		if (label[i] < 0) {
			continue;
		}
		if (c_s.count(label[i]) == 0) {
			c_s[label[i]].resize(dataset[0].size());
		}
		for (int k = 0; k < dataset[0].size(); k++) {
			c_s[label[i]][k] += (c_avg[label[i]][k] - dataset[i][k]) * (c_avg[label[i]][k] - dataset[i][k]);
		}
	}
	//for (int i = 0; i <= cluster_num; i++) {
	for (int i : cluster_center) {
		if (c_count[i] > 0) {
			for (int k = 0; k < dataset[0].size(); k++) {
				c_s[i][k] /= c_count[i];
				c_std[i] += c_s[i][k];
			}
			stdev += sqrt(c_std[i]);
		}
	}
	stdev /= cluster_num;
	vector<int> cluster_center_vec(cluster_center.begin(), cluster_center.end());
	auto density_count = [&dataset, &cluster_center_vec, M, &label, stdev](const vector<int> &cluster) {
		// cluster: index of cluster_center_vec , range [0, cluster_center_vec.size() - 1]
		vector<int64_t> center(dataset[0].size());
		if (cluster.size() == 2) {
			int &lhs = cluster_center_vec[cluster[0]];
			int &rhs = cluster_center_vec[cluster[1]];
			for (int k = 0; k < center.size(); k++) {
				center[k] = (dataset[lhs][k] + dataset[rhs][k]) / 2;
			}
		}
		else {
			int &idx = cluster_center_vec[cluster[0]];
			for (int k = 0; k < center.size(); k++) {
				center[k] = dataset[idx][k];
			}
		}
		int res = 0;
		for (int i = 0; i < M; i++) {
			for (int vindex : cluster) {
				if (label[i] != cluster_center_vec[vindex])
					continue;
				// distance
				int64_t dist = 0;
				for (int k = 0; k < dataset[0].size(); k++) {
					dist += (dataset[i][k] - center[k]) * (dataset[i][k] - center[k]);
				}
				if (sqrt(dist) < stdev)
					res++;
			}
		}
		return res;
	};
	vector<int> density_list(cluster_center_vec.size());
	for (int i = 0; i < density_list.size(); i++) {
		density_list[i] = density_count({ i });
	}
	double Dens_bw = 0;
	for (int i = 0; i < density_list.size(); i++) {
		for (int j = 0; j < density_list.size(); j++) {
			if (i == j)
				continue;
			Dens_bw += double(density_count({ i, j })) / max(density_list[i], density_list[j]);
		}
	}
	Dens_bw /= density_list.size() * (density_list.size() - 1);

	
	vector<double> col_s(dataset[0].size()), col_avg(dataset[0].size());
	int line_count = 0;
	for (int i = 0; i < M; i++) {
		// TODO(sstan): should exclude points that do not belong to any label?
		for (int k = 0; k < dataset[0].size(); k++)
			col_avg[k] += dataset[i][k];
		line_count++;
	}
	for (int k = 0; k < dataset[0].size(); k++)
		col_avg[k] /= line_count;
	for (int i = 0; i < M; i++) {
		for (int k = 0; k < dataset[0].size(); k++)
			col_s[k] += (dataset[i][k] - col_avg[k]) * (dataset[i][k] - col_avg[k]);
	}
	double sigma_s_norm = 0;
	for (int k = 0; k < dataset[0].size(); k++) {
		col_s[k] /= line_count;
		sigma_s_norm += col_s[k];
	}
	sigma_s_norm = sqrt(sigma_s_norm);

	double sum_sigma_norm = 0;
	for (int i : cluster_center) {
		sum_sigma_norm += sqrt(c_std[i]);
	}
	
	double Scat = sum_sigma_norm / (sigma_s_norm * cluster_center.size());

	double S_Dbw = Dens_bw + Scat;
	cout << "S_Dbw = " << S_Dbw << endl;
	return 0;
}