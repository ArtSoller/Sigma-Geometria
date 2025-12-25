#include <iostream>
#include <vector>
#include <fstream>
#include <array>

using namespace std;

#include <time.h>
#include <omp.h>
#include <cmath>

#include <locale.h>
#include <Windows.h>
#include <algorithm>
#include <iomanip>

using namespace std;

int maxiter = 100; // максимальное количество итераций прямая задача
int max_iter = 100; // максимальное количество итераций обратная задача
double eps = 1e-15; // величина требуемой относительной невзяки
int n = 0; // количество конечных элементов
int nv = 0; // количество узлов в сетке
int ne = 0; // количество краевых условий
int nn = 0; // количество главных условий
double r_min = 0; // наименьшая координата по r
double r_max = 0; // наибольшая координата по r
double z_min = 0; // наименьшая координата по z
double z_max = 0; // наибольшая координата по z
int nr_full = 0; // количество точек по r
int nz_full = 0; // количество точек по z
int n1 = 0; // количество 1 краевых
int n2 = 0; // количество 2 краевых
const double VEL_LEN = 100.0;

vector<vector<double>> G(4); // матрица жесткости
vector<vector<double>> M(4); // матрица массы
vector<vector<double>> localA(4); // локальная матрица
vector<double> localb(4); // локальный вектор
vector<vector<double>> grid; // матрица координат узлов
vector<vector<int>> num_elem; // матрица глобальных номеров узлов
vector<vector<int>> edge1; // 1 краевые условия (узел1, узел2, формула)
vector<vector<int>> edge2; // 2 краевые условия (узел1, узел2, формула)
vector<double> localb2(2); // для учета 2 и 3 краевых
vector<vector<double>> locala2(2); // для учета 3 краевых

vector<double> gridr;
vector<double> gridz;
vector<vector<double>> sloy;
vector<vector<double>> sloy_dop;


vector<int> ig; // строки
vector<int> jg; // столбцы
vector<double> al; // нижний треугольник 
vector<double> au; // верхний треугольник
vector<double> di; // диагональ
vector<double> b; // вектор правой части
// для ЛОС
vector<double> L; // нижний треугольник 
vector<double> U; // верхний треугольник
vector<double> D; // диагональ
vector<double> x0;
vector<double> z;
vector<double> r;
vector<double> y;
vector<double> p;
vector<double> t;
// обратная задача
vector<double> true_eps(5); // истинное решение
vector<double> tmp_eps(5); // решение на текущей итерации прямой задачи
vector<double> next_eps(5); // решение на следующей итерации прямой задачи
vector<double> delta_tmp_eps(5); // решение на текущей итерации прямой задачи с прираще-нием
double start_u = 1e-6; // начальное приближение (u0)
double tmp_u; // решение на текущей итерации обратной задачи
double delta_u; // приращение на текущей итерации
double w = 1; // веса
double alpha = 0;
double betta = 1;
double gamma = 1e-4;
double tok = 1; // ток

vector<double> pointz; // контрольные точки по z
vector<double> nz; // количество интервалов по z
vector<double> kz; // коэффициент разрядки по z
vector<int> kz_route; // направление разрядки

vector<double> pointr; // контрольные точки по r
vector<double> nr; // количество интервалов по r
vector<double> kr; // коэффициент разрядки по r
vector<int> kr_route; // направление разрядки

ifstream input("input.txt");
vector<double> qn;
vector<double> qv;
vector<double> q_dop;
vector<vector<double>> grid_n;
vector<vector<int>> num_elem_n;
vector<vector<double>> grid_dop;
vector<vector<int>> num_elem_dop;

double xA = 0;
double yA = 0;
double zA = 0;
double xB = 0;
double yB = 100.0;
double zB = -VEL_LEN;


vector<double_t> initial_params = { -200.0 /*y1*/ };


struct anomaly {
	double x0 = 0.0;
	double x1 = 0.0;
	double y0 = 0.0;
	double y1 = 0.0;
} synthetic_anomaly, generated_anomaly;

struct mesh_buffer {
	size_t r_amount = 0; vector<double> r{};
	vector<size_t> nr{}; vector<double> kr{};
	vector<int8_t> r_direction{};

	size_t z_amount = 0; vector<double> z{};
	vector<size_t> nz{}; vector<double> kz{};
	vector<int8_t> z_direction{};
} main_mesh, dop_mesh;

// подсчет в заданной точке
double 
result_q(double r, double z, vector<double> q, vector<vector<double>> grid, vector<vector<int>> num_elem)
{
	double res = 0;
	// определить, какому элементу принадлежит
	for (int i = 0; i < num_elem.size(); i++)
	{
		double r0 = grid[num_elem[i][0]][0];
		double r1 = grid[num_elem[i][1]][0];
		double z0 = grid[num_elem[i][0]][1];
		double z1 = grid[num_elem[i][2]][1];
		double hr = r1 - r0;
		double hz = z1 - z0;

		if (r >= r0 && r <= r1 && z >= z0 && z <= z1)
		{// вычисление соответствующих базисных функций

			vector<double> psi(4);
			psi[0] = (r1 - r) / hr * (z1 - z) / hz;
			psi[1] = (r - r0) / hr * (z1 - z) / hz;
			psi[2] = (r1 - r) / hr * (z - z0) / hz;
			psi[3] = (r - r0) / hr * (z - z0) / hz;

			res = psi[0] * q[num_elem[i][0]] + psi[1] * q[num_elem[i][1]] + psi[2] * q[num_elem[i][2]] + psi[3] * q[num_elem[i][3]]; // значение V
			return res;
		}
	}
	return res;
}

//коэффициент 
//коэффициент 
double sigma(int num, vector<vector<double>> grid, vector<vector<int>> num_elem)
{
	for (int i = 0; i < sloy.size(); i++)
	{
		// Получаем координаты вершин элемента
		double z0 = grid[num_elem[num][0]][1];  // нижняя левая
		double z1 = grid[num_elem[num][2]][1];  // верхняя левая
		double r0 = grid[num_elem[num][0]][0];  // нижняя левая
		double r1 = grid[num_elem[num][1]][0];  // нижняя правая

		// Координаты центра элемента
		double center_z = (z0 + z1) / 2.0;
		double center_r = (r0 + r1) / 2.0;

		// Получаем границы слоя
		double layer_r_min = sloy[i][0];
		double layer_r_max = sloy[i][1];
		double layer_z_min = sloy[i][2];
		double layer_z_max = sloy[i][3];
		double layer_sigma = sloy[i][4];

		// Проверяем, попадает ли центр элемента в слой
		bool in_z_range = (center_z >= layer_z_min && center_z <= layer_z_max);
		bool in_r_range = (center_r >= layer_r_min && center_r <= layer_r_max);

		if (in_z_range && in_r_range) {
			return layer_sigma;
		}
	}
	return 0;
}

double sigma_dop(int num, vector<vector<double>> grid, vector<vector<int>> num_elem)
{
	for (int i = 0; i < sloy_dop.size(); i++)
	{
		// Получаем координаты вершин элемента
		double z0 = grid[num_elem[num][0]][1];
		double z1 = grid[num_elem[num][2]][1];
		double r0 = grid[num_elem[num][0]][0];
		double r1 = grid[num_elem[num][1]][0];

		// Координаты центра элемента
		double center_z = (z0 + z1) / 2.0;
		double center_r = (r0 + r1) / 2.0;

		// Получаем границы слоя
		double layer_r_min = sloy_dop[i][0];
		double layer_r_max = sloy_dop[i][1];
		double layer_z_min = sloy_dop[i][2];
		double layer_z_max = sloy_dop[i][3];
		double layer_sigma = sloy_dop[i][4];

		// Проверяем, попадает ли центр элемента в слой
		bool in_z_range = (center_z >= layer_z_min && center_z <= layer_z_max);
		bool in_r_range = (center_r >= layer_r_min && center_r <= layer_r_max);

		if (in_z_range && in_r_range) {
			return layer_sigma;
		}
	}
	return 0;
}



// подсчет матрицы жесткости
void GetLocalG(double rp, double zs, double hr, double hz)
{

	double a1 = (hz * rp) / (6 * hr);
	double a2 = (hz) / 12;
	double a3 = (hr * rp) / (6 * hz);
	double a4 = (hr * hr) / (12 * hz);

	G[0][0] = 2 * a1 + 2 * a2 + 2 * a3 + a4;
	G[0][1] = -2 * a1 - 2 * a2 + a3 + a4;
	G[0][2] = a1 + a2 - 2 * a3 - a4;
	G[0][3] = -a1 - a2 - a3 - a4;

	G[1][0] = -2 * a1 - 2 * a2 + a3 + a4;
	G[1][1] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
	G[1][2] = -a1 - a2 - a3 - a4;
	G[1][3] = a1 + a2 - 2 * a3 - 3 * a4;

	G[2][0] = a1 + a2 - 2 * a3 - a4;
	G[2][1] = -a1 - a2 - a3 - a4;
	G[2][2] = 2 * a1 + 2 * a2 + 2 * a3 + a4;
	G[2][3] = -2 * a1 - 2 * a2 + a3 + a4;

	G[3][0] = -a1 - a2 - a3 - a4;
	G[3][1] = a1 + a2 - 2 * a3 - 3 * a4;
	G[3][2] = -2 * a1 - 2 * a2 + a3 + a4;
	G[3][3] = 2 * a1 + 2 * a2 + 2 * a3 + 3 * a4;
}

// расчет начального шага
double step(double min, double max, int n, double k)
{
	if (k != 1)
		return (max - min) * (1 - k) / (1 - pow(k, n));
	else
		return (max - min) / n;
}

static void
sigma_read() {
	ifstream fin("sigma.txt");
	size_t n_sloy; // кол-во слоев
	fin >> n_sloy; sloy.resize(n_sloy);
	for (int i = 0; i < n_sloy; i++) {
		sloy[i].resize(5);
		fin >> sloy[i][0] >> sloy[i][1] >> sloy[i][2] >> sloy[i][3] >> sloy[i][4];
	}
	size_t n_sloy_dop; fin >> n_sloy_dop;
	sloy_dop.resize(n_sloy_dop);
	for (int i = 0; i < n_sloy_dop; i++) {
		sloy_dop[i].resize(5);
		fin >> sloy_dop[i][0] >> sloy_dop[i][1] >> sloy_dop[i][2] >> sloy_dop[i][3] >> sloy_dop[i][4];
	}
	synthetic_anomaly.x0 = sloy_dop[0][0];
	synthetic_anomaly.x1 = sloy_dop[0][1];
	synthetic_anomaly.y0 = sloy_dop[0][2];
	synthetic_anomaly.y1 = sloy_dop[0][3];
}

[[deprecated]]
void SigmaGeneration()
{
	ifstream input("sigma.txt");

	int n_sloy; // кол-во слоев

	input >> n_sloy;

	sloy.resize(n_sloy);

	for (int i = 0; i < n_sloy; i++)
	{
		sloy[i].resize(5);
		// левая и правая границы нижняя граница .. верхняя граница .. sigma
		input >> sloy[i][0] >> sloy[i][1] >> sloy[i][2] >> sloy[i][3] >> sloy[i][4];
	}
	int n_sloy_dop; // кол-во слоев
	input >> n_sloy_dop;
	sloy_dop.resize(n_sloy_dop);

	for (int i = 0; i < n_sloy_dop; i++)
	{
		sloy_dop[i].resize(5);
		// левая и правая границы нижняя граница .. верхняя граница .. sigma
		input >> sloy_dop[i][0] >> sloy_dop[i][1] >> sloy_dop[i][2] >> sloy_dop[i][3] >> sloy_dop[i][4];
	}
}

// Use it once !!! No need to read mesh iteratively.
static void 
read_mesh() {
	ifstream fin("input.txt");

	fin >> main_mesh.r_amount; main_mesh.r.resize(main_mesh.r_amount + 1);
	for (auto& ri : main_mesh.r) fin >> ri;
	main_mesh.nr.resize(main_mesh.r_amount); for (auto& n_i : main_mesh.nr) fin >> n_i;
	main_mesh.kr.resize(main_mesh.r_amount); for (auto& k_i : main_mesh.kr) fin >> k_i;
	main_mesh.r_direction.resize(main_mesh.r_amount); for (auto& d_i : main_mesh.r_direction) fin >> d_i;
	
	fin >> main_mesh.z_amount; main_mesh.z.resize(main_mesh.z_amount + 1);
	for (auto& ri : main_mesh.z) fin >> ri;
	main_mesh.nz.resize(main_mesh.z_amount); for (auto& n_i : main_mesh.nz) fin >> n_i;
	main_mesh.kz.resize(main_mesh.z_amount); for (auto& k_i : main_mesh.kz) fin >> k_i;
	main_mesh.z_direction.resize(main_mesh.z_amount); for (auto& d_i : main_mesh.z_direction) fin >> d_i;
}

static void 
read_dop_mesh() {
	ifstream fin("input_dop.txt");
	fin >> dop_mesh.r_amount; dop_mesh.r.resize(dop_mesh.r_amount + 1);
	for (auto& ri : dop_mesh.r) fin >> ri;
	dop_mesh.nr.resize(dop_mesh.r_amount); for (auto& n_i : dop_mesh.nr) fin >> n_i;
	dop_mesh.kr.resize(dop_mesh.r_amount); for (auto& k_i : dop_mesh.kr) fin >> k_i;
	dop_mesh.r_direction.resize(dop_mesh.r_amount); for (auto& d_i : dop_mesh.r_direction) fin >> d_i;

	fin >> dop_mesh.z_amount; dop_mesh.z.resize(dop_mesh.z_amount + 1);
	for (auto& ri : dop_mesh.z) fin >> ri;
	dop_mesh.nz.resize(dop_mesh.z_amount); for (auto& n_i : dop_mesh.nz) fin >> n_i;
	dop_mesh.kz.resize(dop_mesh.z_amount); for (auto& k_i : dop_mesh.kz) fin >> k_i;
	dop_mesh.z_direction.resize(dop_mesh.z_amount); for (auto& d_i : dop_mesh.z_direction) fin >> d_i;

	synthetic_anomaly.x0 = dop_mesh.r[1];
	synthetic_anomaly.x1 = dop_mesh.r[2];
	synthetic_anomaly.y0 = dop_mesh.z[1];
	synthetic_anomaly.y1 = dop_mesh.z[2];
}

static void 
main_grid_generation() {
	ofstream node("node.txt");
	ofstream elem("elem.txt");

	// построение r
	pointr = {}; pointr.resize(main_mesh.r.size()); copy(main_mesh.r.begin(), main_mesh.r.end(), pointr.begin());
	for (int i = 0; i < pointr.size() - 1; i++)
	{
		// начальный шаг
		double hr0 = step(pointr[i], pointr[i + 1], main_mesh.nr[i], main_mesh.kr[i]);

		if (main_mesh.r_direction[i] == 1) // возрастает шаг
		{
			double r = pointr[i];
			while (pointr[i + 1] - r > 1e-6)
			{
				gridr.push_back(r);
				r += hr0;
				hr0 *= main_mesh.kr[i];
			}
		}
		else // убывает
		{
			double r = pointr[i + 1] - hr0;
			while (r - pointr[i] > 1e-6)
			{
				gridr.push_back(r);
				hr0 *= main_mesh.kr[i];
				r -= hr0;
			}
			gridr.push_back(pointr[i]);
		}
	}
	//gridr.push_back(pointr[pointr.size() - 1]); // последнюю точку отдельно 
	gridr.push_back(pointr.back()); // последнюю точку отдельно 

	// построение z
	pointz = {}; pointz.resize(main_mesh.z.size()); copy(main_mesh.z.begin(), main_mesh.z.end(), pointz.begin());
	for (int i = 0; i < pointz.size() - 1; i++)
	{
		// начальный шаг
		double hz0 = step(pointz[i], pointz[i + 1], main_mesh.nz[i], main_mesh.kz[i]);

		if (main_mesh.z_direction[i] == 1) // возрастает шаг
		{
			double z = pointz[i];
			while (pointz[i + 1] - z > 1e-6)
			{
				gridz.push_back(z);
				z += hz0;
				hz0 *= main_mesh.kz[i];
			}
		}
		else // убывает
		{
			double z = pointz[i + 1] - hz0;
			while (z - pointz[i] > 1e-6)
			{
				gridz.push_back(z);
				hz0 *= main_mesh.kz[i];
				z -= hz0;
			}
			gridz.push_back(pointz[i]);
		}
	}
	//gridz.push_back(pointz[pointz.size() - 1]); // последнюю точку отдельно 
	gridz.push_back(pointz.back()); // последнюю точку отдельно 

	sort(gridz.begin(), gridz.end());
	sort(gridr.begin(), gridr.end());
	nv = gridr.size() * gridz.size();
	grid.resize(nv);

	z_min = gridz[0];
	r_max = gridr[gridz.size() - 1];

	node << nv << endl;
	int i = 0;

	for (int iy = 0; iy < gridz.size(); iy++)
		for (int ix = 0; ix < gridr.size(); ix++)
		{
			grid[i].resize(2);
			grid[i][0] = gridr[ix];
			grid[i][1] = gridz[iy];
			node << i << " " << grid[i][0] << " " << grid[i][1] << endl;
			i++;
		}

	nr_full = gridr.size() - 1;
	nz_full = gridz.size() - 1;
	n = nr_full * nz_full; // количество конечных элементов
	num_elem.resize(n);
	int tmp_num = 0;
	for (int j = 0; j < nz_full; j++)
		for (int i = 0; i < nr_full; i++)
		{
			num_elem[tmp_num].resize(4);
			num_elem[tmp_num][0] = i + j * (nr_full + 1); // 1 вершина
			num_elem[tmp_num][1] = i + j * (nr_full + 1) + 1; // 2 вершина
			num_elem[tmp_num][2] = i + (j + 1) * (nr_full + 1); // 3 вершина
			num_elem[tmp_num][3] = i + (j + 1) * (nr_full + 1) + 1; // 4 вершина
			// вывод в файл
			elem << tmp_num << " " << num_elem[tmp_num][0] << " " << num_elem[tmp_num][1] << " " << num_elem[tmp_num][2] << " " <<
				num_elem[tmp_num][3] << ' ' << sigma(tmp_num, grid, num_elem) << endl;
			tmp_num++;
		}
	node.close();
	elem.close();
}

static void 
dop_grid_generation() {
	ofstream fout_node_dop("node_dop.txt");
	ofstream fout_elem_dop("elem_dop.txt");

	// построение r
	pointr = {}; pointr.resize(dop_mesh.r.size()); copy(dop_mesh.r.begin(), dop_mesh.r.end(), pointr.begin());
	for (int i = 0; i < pointr.size() - 1; i++)
	{
		// начальный шаг
		double hr0 = step(pointr[i], pointr[i + 1], dop_mesh.nr[i], dop_mesh.kr[i]);

		if (dop_mesh.r_direction[i] == 1) // возрастает шаг
		{
			double r = pointr[i];
			while (pointr[i + 1] - r > 1e-6)
			{
				gridr.push_back(r);
				r += hr0;
				hr0 *= dop_mesh.kr[i];
			}
		}
		else // убывает
		{
			double r = pointr[i + 1] - hr0;
			while (r - pointr[i] > 1e-6)
			{
				gridr.push_back(r);
				hr0 *= dop_mesh.kr[i];
				r -= hr0;
			}
			gridr.push_back(pointr[i]);
		}
	}
	//gridr.push_back(pointr[pointr.size() - 1]); // последнюю точку отдельно 
	gridr.push_back(pointr.back()); // последнюю точку отдельно 

	// построение z
	pointz = {}; pointz.resize(dop_mesh.z.size()); copy(dop_mesh.z.begin(), dop_mesh.z.end(), pointz.begin());
	for (int i = 0; i < pointz.size() - 1; i++)
	{
		// начальный шаг
		double hz0 = step(pointz[i], pointz[i + 1], dop_mesh.nz[i], dop_mesh.kz[i]);

		if (dop_mesh.z_direction[i] == 1) // возрастает шаг
		{
			double z = pointz[i];
			while (pointz[i + 1] - z > 1e-6)
			{
				gridz.push_back(z);
				z += hz0;
				hz0 *= dop_mesh.kz[i];
			}
		}
		else // убывает
		{
			double z = pointz[i + 1] - hz0;
			while (z - pointz[i] > 1e-6)
			{
				gridz.push_back(z);
				hz0 *= dop_mesh.kz[i];
				z -= hz0;
			}
			gridz.push_back(pointz[i]);
		}
	}
	//gridz.push_back(pointz[pointz.size() - 1]); // последнюю точку отдельно 
	gridz.push_back(pointz.back()); // последнюю точку отдельно 

	sort(gridz.begin(), gridz.end());
	sort(gridr.begin(), gridr.end());
	nv = gridr.size() * gridz.size();
	grid.resize(nv);

	z_min = gridz[0];
	r_max = gridr[gridz.size() - 1];

	fout_node_dop << nv << endl;
	int i = 0;

	for (int iy = 0; iy < gridz.size(); iy++)
		for (int ix = 0; ix < gridr.size(); ix++)
		{
			grid[i].resize(2);
			grid[i][0] = gridr[ix];
			grid[i][1] = gridz[iy];
			fout_node_dop << i << " " << grid[i][0] << " " << grid[i][1] << endl;
			i++;
		}

	nr_full = gridr.size() - 1;
	nz_full = gridz.size() - 1;
	n = nr_full * nz_full; // количество конечных элементов
	num_elem.resize(n);
	int tmp_num = 0;
	for (int j = 0; j < nz_full; j++)
		for (int i = 0; i < nr_full; i++)
		{
			num_elem[tmp_num].resize(4);
			num_elem[tmp_num][0] = i + j * (nr_full + 1); // 1 вершина
			num_elem[tmp_num][1] = i + j * (nr_full + 1) + 1; // 2 вершина
			num_elem[tmp_num][2] = i + (j + 1) * (nr_full + 1); // 3 вершина
			num_elem[tmp_num][3] = i + (j + 1) * (nr_full + 1) + 1; // 4 вершина
			// вывод в файл
			fout_elem_dop << tmp_num << " " << num_elem[tmp_num][0] << " " << num_elem[tmp_num][1] << " " << num_elem[tmp_num][2] << " " <<
				num_elem[tmp_num][3] << ' ' << sigma(tmp_num, grid, num_elem) << endl;
			tmp_num++;
		}

	fout_node_dop.close();
	fout_elem_dop.close();
}

// cетка
[[deprecated]]
void GridGeneration()
{
	ofstream node("node.txt");
	ofstream elem("elem.txt");

	// построение r
	for (int i = 0; i < pointr.size() - 1; i++)
	{
		// начальный шаг
		double hr0 = step(pointr[i], pointr[i + 1], nr[i], kr[i]);

		if (kr_route[i] == 1) // возрастает шаг
		{
			double r = pointr[i];
			while (pointr[i + 1] - r > 1e-6)
			{
				gridr.push_back(r);
				r += hr0;
				hr0 *= kr[i];
			}
		}
		else // убывает
		{
			double r = pointr[i + 1] - hr0;
			while (r - pointr[i] > 1e-6)
			{
				gridr.push_back(r);
				hr0 *= kr[i];
				r -= hr0;
			}
			gridr.push_back(pointr[i]);
		}
	}
	gridr.push_back(pointr[pointr.size() - 1]); // последнюю точку отдельно 

	// построение z
	for (int i = 0; i < pointz.size() - 1; i++)
	{
		// начальный шаг
		double hz0 = step(pointz[i], pointz[i + 1], nz[i], kz[i]);

		if (kz_route[i] == 1) // возрастает шаг
		{
			double z = pointz[i];
			while (pointz[i + 1] - z > 1e-6)
			{
				gridz.push_back(z);
				z += hz0;
				hz0 *= kz[i];
			}
		}
		else // убывает
		{
			double z = pointz[i + 1] - hz0;
			while (z - pointz[i] > 1e-6)
			{
				gridz.push_back(z);
				hz0 *= kz[i];
				z -= hz0;
			}
			gridz.push_back(pointz[i]);
		}
	}
	gridz.push_back(pointz[pointz.size() - 1]); // последнюю точку отдельно 

	sort(gridz.begin(), gridz.end());
	sort(gridr.begin(), gridr.end());
	nv = gridr.size() * gridz.size();
	grid.resize(nv);

	z_min = gridz[0];
	r_max = gridr[gridz.size() - 1];

	node << nv << endl;
	int i = 0;

	for (int iy = 0; iy < gridz.size(); iy++)
		for (int ix = 0; ix < gridr.size(); ix++)
		{
			grid[i].resize(2);
			grid[i][0] = gridr[ix];
			grid[i][1] = gridz[iy];
			node << i << " " << grid[i][0] << " " << grid[i][1] << endl;
			i++;
		}

	nr_full = gridr.size() - 1;
	nz_full = gridz.size() - 1;
	n = nr_full * nz_full; // количество конечных элементов
	num_elem.resize(n);
	int tmp_num = 0;
	for (int j = 0; j < nz_full; j++)
		for (int i = 0; i < nr_full; i++)
		{
			num_elem[tmp_num].resize(4);
			num_elem[tmp_num][0] = i + j * (nr_full + 1); // 1 вершина
			num_elem[tmp_num][1] = i + j * (nr_full + 1) + 1; // 2 вершина
			num_elem[tmp_num][2] = i + (j + 1) * (nr_full + 1); // 3 вершина
			num_elem[tmp_num][3] = i + (j + 1) * (nr_full + 1) + 1; // 4 вершина
			// вывод в файл
			elem << tmp_num << " " << num_elem[tmp_num][0] << " " << num_elem[tmp_num][1] << " " << num_elem[tmp_num][2] << " " <<
				num_elem[tmp_num][3] << ' ' << sigma(tmp_num, grid, num_elem) << endl;
			tmp_num++;
		}
}

// портрет матрицы
void MatrixPortrait(vector<vector<int>> num_elem)
{
	vector<vector<int>> list(nv); // вспомогательный массив по количеству узлов
	vector<int> tmp(4); // вспомогательный вектор вершин одного элемента
	for (int i = 0; i < nv; i++)
		list[i];
	int number = 0; // количество элементов для массива ig
	for (int i = 0; i < num_elem.size(); i++) //идем по конечным элементам
	{
		for (int j = 0; j < 4; j++)
			tmp[j] = num_elem[i][j];
		reverse(tmp.begin(), tmp.end()); //сортируем по убывaнию, т.к. портрет сим-метричный (в нашем случае переворачиваем т.к. все изначально отсортировано)
		for (int j = 0; j < 4; j++)
			for (int k = j + 1; k < 4; k++)
			{
				int flag = 1;
				for (int p = 0; p < list[tmp[j]].size() && flag; p++)
					if (list[tmp[j]][p] == tmp[k]) flag = 0; // если элемент уже есть в списке смежности
				if (flag)
				{
					list[tmp[j]].push_back(tmp[k]); // если в списке еще нет такого элемента, то добавляем
					number++; // увеличиваем количество для jg
				}
			}
	}
	for (int i = 0; i < nv; i++)
		sort(list[i].begin(), list[i].end()); //сортируем по возрастанию 

	ig.resize(nv + 1);

	ig[0] = 0;
	ig[1] = 0;
	for (int i = 2; i < nv + 1; i++)
		ig[i] = ig[i - 1] + list[i - 1].size();

	for (int i = 0; i < nv; i++)
		for (int j = 0; j < list[i].size(); j++)
			jg.push_back(list[i][j]);

	//выделение памяти для работы с векторами
	al.resize(number);
	au.resize(number);
	di.resize(nv);
	b.resize(nv);
	//ЛОС
	L.resize(number);
	U.resize(number);
	D.resize(nv);
	x0.resize(nv);
	y.resize(nv);
	z.resize(nv);
	r.resize(nv);
	p.resize(nv);
	t.resize(nv);
	locala2[0].resize(2);
	locala2[1].resize(2);

	for (int i = 0; i < 4; i++)
	{
		G[i].resize(4);
		M[i].resize(4);
		localA[i].resize(4);
	}
}

// подсчет локальной матрицы
void LocalMatrix(int num, vector<vector<double>> grid, vector<vector<int>> num_elem)
{
	//вершины прямоугольника
	int a = num_elem[num][0];
	int b = num_elem[num][1];
	int c = num_elem[num][2];
	int d = num_elem[num][3];
	double rp, zs, hr, hz;
	rp = grid[a][0];
	zs = grid[a][1];
	hr = grid[b][0] - grid[a][0];
	hz = grid[c][1] - grid[b][1];

	GetLocalG(rp, zs, hr, hz);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			localA[i][j] = sigma(num, grid, num_elem) * G[i][j];
		}
}

// подсчет локальной матрицы
void LocalMatrix_dop(int num, vector<vector<double>> grid, vector<vector<int>> num_elem)
{
	for (int i = 0; i < 4; i++)
	{
		localb[i] = 0;
	}

	//вершины прямоугольника
	int a = num_elem[num][0];
	int b = num_elem[num][1];
	int c = num_elem[num][2];
	int d = num_elem[num][3];
	double rp, zs, hr, hz;
	rp = grid[a][0];
	zs = grid[a][1];
	hr = grid[b][0] - grid[a][0];
	hz = grid[c][1] - grid[b][1];

	vector<double> q; // правая часть задачи V_n
	q.push_back(result_q(grid[a][0], grid[a][1], qn, grid_n, num_elem_n));
	q.push_back(result_q(grid[b][0], grid[b][1], qn, grid_n, num_elem_n));
	q.push_back(result_q(grid[c][0], grid[c][1], qn, grid_n, num_elem_n));
	q.push_back(result_q(grid[d][0], grid[d][1], qn, grid_n, num_elem_n));

	GetLocalG(rp, zs, hr, hz);

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			localA[i][j] = sigma_dop(num, grid, num_elem) * G[i][j];
		}

	for (int i = 0; i < 4; i++)
		for (int j = 0; j < 4; j++)
		{
			localb[i] += (sigma(num, grid, num_elem) - sigma_dop(num, grid, num_elem)) * G[i][j] * q[j];
		}
}

// добавление локальной матрицы в глобальную
void AddInGlobal(int num, vector<vector<int>> num_elem)
{
	// диагональ
	for (int i = 0; i < 4; i++)
		di[num_elem[num][i]] += localA[i][i];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < i; j++)
		{
			int igg = num_elem[num][i];
			int jgg = num_elem[num][j];
			if (igg < jgg)
			{
				igg = jgg;
				jgg = num_elem[num][i];
			}
			int index = ig[igg];
			int flag = 1;
			for (; index < ig[igg + 1] && flag; index++)
				if (jg[index] == jgg) flag = 0;
			index--;
			al[index] += localA[i][j];
			au[index] += localA[j][i];
		}
	}
}

// добавление локальной матрицы в глобальную
void AddInGlobal_dop(int num, vector<vector<int>> num_elem)
{
	// диагональ
	for (int i = 0; i < 4; i++)
		di[num_elem[num][i]] += localA[i][i];

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < i; j++)
		{
			int igg = num_elem[num][i];
			int jgg = num_elem[num][j];
			if (igg < jgg)
			{
				igg = jgg;
				jgg = num_elem[num][i];
			}
			int index = ig[igg];
			int flag = 1;
			for (; index < ig[igg + 1] && flag; index++)
				if (jg[index] == jgg) flag = 0;
			index--;
			al[index] += localA[i][j];
			au[index] += localA[j][i];
		}
	}

	for (int i = 0; i < 4; i++)
		b[num_elem[num][i]] += localb[i];
}

// учет первых краевых условий
void First(vector<vector<double>> grid)
{
	for (int i = 0; i < grid.size(); i++) {
		// нижняя граница
		if (grid[i][1] == z_min) {
			di[i] = 1e+15;
			b[i] = 0;
		}
		// правая граница
		if (grid[i][0] == r_max) {
			di[i] = 1e+15;
			b[i] = 0;
		}
	}
}

#pragma region IMMUTABLE 
void LU()
{
	for (int i = 0; i < ig[nv]; i++)
	{
		L[i] = al[i];
		U[i] = au[i];
	}
	for (int i = 0; i < nv; i++)
		D[i] = di[i];

	for (int i = 0; i < nv; i++)
	{
		int ia0 = ig[i];
		int ia1 = ig[i + 1];
		double sd = 0;
		for (int k = ia0; k < ia1; k++)
		{
			double s1 = 0;
			double s2 = 0;
			int j = jg[k]; // номер столбца(строки) элемента al(au)
			int ja0 = ig[j]; // номер, с которого начинается элементы столб-ца(строки)
			int ja1 = ig[j + 1]; // номер, с которого начинаются элементы нового столбца(строки)
			int ki = ia0;
			int kj = ja0;
			while (ki < k && kj < ja1)
				if (jg[ki] == jg[kj])
				{
					s1 += L[ki] * U[kj];
					s2 += L[kj] * U[ki];
					ki++;
					kj++;
				}
				else if (jg[ki] > jg[kj])
					kj++;
				else
					ki++;
			L[k] = L[k] - s1;
			U[k] = U[k] - s2;
			U[k] = U[k] / D[jg[k]];
			sd += U[k] * L[k];
		}
		D[i] -= sd;
	}
}

void LSolve(vector<double>& x1, vector<double> x2)
{
	x1[0] = x2[0] / D[0];
	for (int i = 1; i < nv; i++)
	{
		double s = 0;
		for (int j = ig[i]; j < ig[i + 1]; j++)
			s += L[j] * x1[jg[j]];
		x1[i] = x2[i] - s;
		x1[i] = x1[i] / D[i];
	}
}

void USolve(vector<double>& x1, vector<double> x2)
{
	for (int i = 0; i < nv; i++)
		x1[i] = x2[i];
	x1[nv - 1] = x1[nv - 1] / D[nv - 1];
	for (int i = nv - 1; i >= 1; i--)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
			x1[jg[j]] -= U[j] * x1[i];
	}
}

// умножения вектора на матрицу A
void multiplication_matrix_on_vector(vector<double> a, vector<double>& b)
{
	for (int i = 0; i < nv; i++)
		b[i] = 0;
	for (int i = 0; i < nv; i++)
	{
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			b[i] += al[j] * a[jg[j]];
			b[jg[j]] += au[j] * a[i];
		}
		b[i] += di[i] * a[i];
	}
}

// умножение двух векторов
double vectors_multiplication(vector<double> v1, vector<double> v2)
{
	double s = 0;
	for (int i = 0; i < nv; i++)
		s += v1[i] * v2[i];
	return s;
}

// норма вектора
double norma(vector<double> vector)
{
	double s = 0;
	for (int i = 0; i < nv; i++)
		s += vector[i] * vector[i];
	return(sqrt(s));
}
// d = a + b * c
void summ(vector<double> a, double b, vector<double> c, vector<double>& d)
{
	for (int i = 0; i < nv; i++)
		d[i] = a[i] + b * c[i];
}

void LOS()
{
	double alpha, betta;
	double pk_1_rk_1, pk_1_pk_1;
	// нулевое начальное приближение
	for (int i = 0; i < nv; i++)
		x0[i] = 0;
	multiplication_matrix_on_vector(x0, y); // Ax0
	summ(b, -1, y, r); // r0 = f - Ax0
	LSolve(r, r); // r0 = (f - Ax0) / L
	USolve(z, r); // z0 = r0 / U
	multiplication_matrix_on_vector(z, y); // Az0
	LSolve(p, y); // p0 = Az0 / L
	for (int k = 1; k < maxiter; k++)
	{
		pk_1_rk_1 = vectors_multiplication(p, r); // (p_(k-1),r_(k-1))
		pk_1_pk_1 = vectors_multiplication(p, p); // (p_(k-1),p_(k-1))
		alpha = pk_1_rk_1 / pk_1_pk_1; // alpha_k = (p_(k-1),r_(k-1)) / (p_(k-1),p_(k-1))
		summ(x0, alpha, z, x0); // x_k = x_(k-1) + alpha_k * z_(k-1)
		summ(r, -alpha, p, r); // r_k = r_(k-1) - alpha_k * p_(k-1)

		USolve(t, r); // t = r_k/ U
		multiplication_matrix_on_vector(t, y); // y = A * r_k / U
		LSolve(t, y); // t = A * r_k / UL
		betta = -vectors_multiplication(p, t) / pk_1_pk_1; // betta_k = (p_k-1,L_1*A*U_1*r_k) / (p_(k-1),p_(k-1))
		USolve(y, r); // y = r_k / U
		summ(y, betta, z, z);//	z_k = r_k / U + betta_k * z_(k-1)
		summ(t, betta, p, p);//p_k = L_1*A*U_1* r_k + betta_k * p_(k-1)
		if (vectors_multiplication(r, r) < eps) // (r_k, r_k) < e
		{
			return;
		}
	}
}
#pragma endregion IMMUTABLE

// прямая задача
vector<double> direct_task()
{
	for (int i = 0; i < num_elem_n.size(); i++)
	{
		LocalMatrix(i, grid_n, num_elem_n);
		AddInGlobal(i, num_elem_n);
	}
	//b[0] = 1;


	// Make VEL
	//for (int i = 0; i < grid_n.size(); i++) {
	//	if (grid_n[i][0] == 0 && grid_n[i][1] == 0) {
	//		b[i] = tok;
	//	}
	//}

	if (num_elem.size() != 0) {
		for (int i(0); i < num_elem.size(); ++i) {
			if (grid_n[num_elem[i][0]][0] == 0.0 && (grid_n[num_elem[i][0]][1] <= zB && zB <= grid_n[num_elem[i][2]][1])) {
				b[num_elem[i][0]] = -tok * (grid_n[num_elem[i][2]][1] - zB) / (grid_n[num_elem[i][2]][1] - grid_n[num_elem[i][0]][1]);
				b[num_elem[i][2]] = -tok * (zB - grid_n[num_elem[i][0]][1]) / (grid_n[num_elem[i][2]][1] - grid_n[num_elem[i][0]][1]);
			}
			else if (grid_n[num_elem[i][0]][0] == 0.0 && (grid_n[num_elem[i][1]][1] <= zA && zA <= grid_n[num_elem[i][3]][1])) {
				b[num_elem[i][0]] = tok * (grid_n[num_elem[i][2]][1] - zA) / (grid_n[num_elem[i][2]][1] - grid_n[num_elem[i][0]][1]);
				b[num_elem[i][2]] = tok * (zA - grid_n[num_elem[i][0]][1]) / (grid_n[num_elem[i][2]][1] - grid_n[num_elem[i][0]][1]);
			}
		}
	}
	else {
		for (int i(0); i < num_elem_n.size(); ++i) {
			if (grid_n[num_elem_n[i][0]][0] == 0.0 && (grid_n[num_elem_n[i][0]][1] <= zB && zB <= grid_n[num_elem_n[i][2]][1])) {
				b[num_elem_n[i][0]] = -tok * (grid_n[num_elem_n[i][2]][1] - zB) / (grid_n[num_elem_n[i][2]][1] - grid_n[num_elem_n[i][0]][1]);
				b[num_elem_n[i][2]] = -tok * (zB - grid_n[num_elem_n[i][0]][1]) / (grid_n[num_elem_n[i][2]][1] - grid_n[num_elem_n[i][0]][1]);
			}
			else if (grid_n[num_elem_n[i][0]][0] == 0.0 && (grid_n[num_elem_n[i][1]][1] <= zA && zA <= grid_n[num_elem_n[i][3]][1])) {
				b[num_elem_n[i][0]] = tok * (grid_n[num_elem_n[i][2]][1] - zA) / (grid_n[num_elem_n[i][2]][1] - grid_n[num_elem_n[i][0]][1]);
				b[num_elem_n[i][2]] = tok * (zA - grid_n[num_elem_n[i][0]][1]) / (grid_n[num_elem_n[i][2]][1] - grid_n[num_elem_n[i][0]][1]);
			}
		}
	}

	First(grid_n);

	LU();
	LOS();
	//user_pardiso(ig, jg, di, al, au, b, x0);
	// меняем значения для нового слоя 
	for (int i = 0; i < grid_n.size(); i++)
	{
		di[i] = 0;
		b[i] = 0;
	}
	for (int i = 0; i < al.size(); i++)
	{
		al[i] = 0;
		au[i] = 0;
	}
	ofstream q("qn.txt");

	for (int i = 0; i < grid_n.size(); i++)
		q << scientific << setprecision(15) << x0[i] << endl;
	return x0;
}

// прямая задача
vector<double> direct_task_dop()
{
	for (int i = 0; i < num_elem_dop.size(); i++)
	{
		LocalMatrix_dop(i, grid_dop, num_elem_dop);
		AddInGlobal_dop(i, num_elem_dop);
	}
	First(grid_dop);
	LU();
	LOS();
	//user_pardiso(ig, jg, di, al, au, b, x0);
	// меняем значения для нового слоя 
	for (int i = 0; i < grid_dop.size(); i++)
	{
		di[i] = 0;
		b[i] = 0;
	}
	for (int i = 0; i < al.size(); i++)
	{
		al[i] = 0;
		au[i] = 0;
	}

	ofstream q("q+.txt");

	for (int i = 0; i < grid_dop.size(); i++)
		q << scientific << setprecision(15) << x0[i] << endl;

	return x0;
}

// SPORNO
double result_xyz_q(double x, double y, vector<double> q, vector<vector<double>> grid, vector<vector<int>> num_elem) {
	double r = sqrt(x * x + y * y);
	double z = 0.0;

	double r1 = sqrt((x - xA) * (x - xA) + (y - yA) * (y - yA));
	double r2 = sqrt((x - xB) * (x - xB) + (y - yB) * (y - yB));

	//double z1 = 0.0;
	//double z2 = -zB;
	return result_q(r, z, q, grid, num_elem); // -result_q(r2, z2, q, grid, num_elem);
	//return result_q(r1, 0.0, q, grid, num_elem) - result_q(r2, 0.0, q, grid, num_elem);
}

void result_function_q(vector<double> q, vector<vector<double>> grid, vector<vector<int>> num_elem) {
	cout << scientific << setprecision(8) << result_xyz_q(100, 0, q, grid, num_elem) << endl;
	cout << scientific << setprecision(8) << result_xyz_q(150, 0, q, grid, num_elem) << endl;
	cout << scientific << setprecision(8) << result_xyz_q(200, 0, q, grid, num_elem) << endl;
	cout << scientific << setprecision(8) << result_xyz_q(250, 0, q, grid, num_elem) << endl;
	cout << scientific << setprecision(8) << result_xyz_q(300, 0, q, grid, num_elem) << endl;
	cout << endl;
}

void result_function_q(vector<double>& vec, vector<double> q, vector<vector<double>> grid, vector<vector<int>> num_elem) {
	vec[0] = result_xyz_q(100, 0, q, grid, num_elem);
	vec[1] = result_xyz_q(150, 0, q, grid, num_elem);
	vec[2] = result_xyz_q(200, 0, q, grid, num_elem);
	vec[3] = result_xyz_q(250, 0, q, grid, num_elem);
	vec[4] = result_xyz_q(300, 0, q, grid, num_elem);
}

double derivative(double a, double b) {
	return (b - a) / (0.05 * tmp_u);
}

void clearAllVectors()
{
	// Основные матрицы и векторы
	G.clear();
	G.resize(4);

	M.clear();
	M.resize(4);

	localA.clear();
	localA.resize(4);

	localb.clear();
	localb.resize(4);

	grid.clear();
	num_elem.clear();
	edge1.clear();
	edge2.clear();

	localb2.clear();
	localb2.resize(2);

	locala2.clear();
	locala2.resize(2);

	// Векторы для СЛАУ
	ig.clear();
	jg.clear();
	al.clear();
	au.clear();
	di.clear();
	b.clear();

	// Векторы для ЛОС
	L.clear();
	U.clear();
	D.clear();
	x0.clear();
	z.clear();
	r.clear();
	y.clear();
	p.clear();
	t.clear();

	// Векторы для сетки
	gridr.clear();
	gridz.clear();
	pointz.clear();
	nz.clear();
	kz.clear();
	kz_route.clear();
	pointr.clear();
	nr.clear();
	kr.clear();
	kr_route.clear();
}

static void 
field_selection() {
	//SigmaGeneration();

	cout << "Synthetic data generation started" << endl;
	main_grid_generation();
	//GridGeneration();

	// сохранение сетки 
	for (int i = 0; i < grid.size(); i++) {
		grid_n.push_back(grid[i]);
	}
	for (int i = 0; i < num_elem.size(); i++) {
		num_elem_n.push_back(num_elem[i]);
	}

	MatrixPortrait(num_elem_n);
	qn = direct_task(); // прямая задача

	// обнуление параметров
	clearAllVectors();

	dop_grid_generation();
	//GridGeneration();
	// сохранение сетки 
	for (int i = 0; i < grid.size(); i++) {

		grid_dop.push_back(grid[i]);
	}
	for (int i = 0; i < num_elem.size(); i++) {

		num_elem_dop.push_back(num_elem[i]);
	}
	MatrixPortrait(num_elem_dop);
	q_dop = direct_task_dop(); // прямая задача дополнительная
	qv.resize(qn.size());

	// итоговый результат 
	for (int i = 0; i < grid_n.size(); i++)
	{
		qv[i] = result_q(grid_n[i][0], grid_n[i][1], qn, grid_n, num_elem_n) +
			result_q(grid_n[i][0], grid_n[i][1], q_dop, grid_dop, num_elem_dop);
	}
	cout << "Synthetic data generated" << endl;
}



void field_selection_direct_task()
{
	clearAllVectors();
	nv = grid_n.size();

	grid_n = {}; num_elem_n = {};
	main_grid_generation();
	//GridGeneration();

	// сохранение сетки 
	for (int i = 0; i < grid.size(); i++) {
		grid_n.push_back(grid[i]);
	}
	for (int i = 0; i < num_elem.size(); i++) {
		num_elem_n.push_back(num_elem[i]);
	}



	MatrixPortrait(num_elem_n);
	qn = direct_task(); // прямая задача
	// обнуление параметров
	clearAllVectors();
	nv = grid_dop.size();

	dop_grid_generation();
	//GridGeneration();
	// сохранение сетки
	grid_dop = {}; num_elem_dop = {};
	for (int i = 0; i < grid.size(); i++) {

		grid_dop.push_back(grid[i]);
	}
	for (int i = 0; i < num_elem.size(); i++) {

		num_elem_dop.push_back(num_elem[i]);
	}


	MatrixPortrait(num_elem_dop);
	q_dop = direct_task_dop(); // прямая задача дополнительная
	// итоговый результат 
	for (int i = 0; i < grid_n.size(); i++)
	{
		qv[i] = result_q(grid_n[i][0], grid_n[i][1], qn, grid_n, num_elem_n) +
			result_q(grid_n[i][0], grid_n[i][1], q_dop, grid_dop, num_elem_dop);
		//qv[i] = result_q(grid_n[i][0], grid_n[i][1], qn, grid_n, num_elem_n);
	}
}

void inverse_problem() {
	cout << "\tGenerated synthetic data on receivers:" << endl;
	result_function_q(true_eps, qv, grid_n, num_elem_n);
	result_function_q(qv, grid_n, num_elem_n);
	//true_eps[0] *= 1.05; // зашумление
	decltype(initial_params) template_params(initial_params.begin(), initial_params.end());
	//tmp_u = start_u;
	tmp_u = initial_params[0];
	for (int iter = 0; iter < max_iter; iter++) {
		double sum = 0; // левая часть уравнения
		double f = 0; // правая часть уравнения
		//tok = tmp_u;
		//sloy_dop[1][2] = tmp_u;


		// ACHTUNG !!!
		sloy_dop[0][3] = tmp_u;
		dop_mesh.z[2] = sloy_dop[0][3];


		// приращение
		field_selection_direct_task();
		result_function_q(tmp_eps, qv, grid_n, num_elem_n);
		result_function_q(qv, grid_n, num_elem_n);

		//tok = 1.05 * tmp_u;
		//sloy_dop[1][2] = 1.05 * tmp_u;
		sloy_dop[0][3] = 1.05 * tmp_u;
		dop_mesh.z[2] = sloy_dop[0][3];


		field_selection_direct_task();
		result_function_q(delta_tmp_eps, qv, grid_n, num_elem_n);
		result_function_q(qv, grid_n, num_elem_n);

		// слау
		for (int i = 0; i < 5; i++) {
			sum += w * w * pow(derivative(delta_tmp_eps[i], tmp_eps[i]), 2);
		}
		sum += alpha;
		for (int i = 0; i < 5; i++) {
			f -= w * w * derivative(delta_tmp_eps[i], tmp_eps[i]) * (true_eps[i] - tmp_eps[i]);
		}
		f -= alpha * (tmp_u - initial_params[0]);

		delta_u = f / sum;
		// итерация
		double J_prev = 0;
		double J_next = 0;
		for (int j = 0; j < 5; j++) {
			J_prev = 0;
			for (int i = 0; i < 5; i++) {
				//w = 1 / tmp_eps[i];
				J_prev += pow(w * (true_eps[i] - tmp_eps[i]), 2);
			}
			J_next = 0;
			//tok = delta_u * betta + tmp_u;
			//sloy_dop[1][2] = delta_u * betta + tmp_u;
			sloy_dop[0][3] = delta_u * betta + tmp_u;
			dop_mesh.z[2] = sloy_dop[0][3];


			field_selection_direct_task();
			result_function_q(next_eps, qv, grid_n, num_elem_n);
			result_function_q(qv, grid_n, num_elem_n);

			for (int i = 0; i < 5; i++) {
				//w = 1 / next_eps[i];
				J_next += pow(w * (true_eps[i] - next_eps[i]), 2);
			}
			if (J_next - J_prev > 0) {
				betta /= 2;
			}
			else {
				break;
			}
		}
		betta = 1;
		//cout << "iter = " << iter + 1 << setprecision(15) << " tok = " << tok << endl;
		//cout << "iter = " << iter + 1 << setprecision(15) << " sigma = " << sloy_dop[1][2] << endl;
		cout << "iter = " << iter + 1 << " x0 = " << sloy_dop[0][0] << endl;
		cout << "iter = " << iter + 1 << " x1 = " << sloy_dop[0][1] << endl;
		cout << "iter = " << iter + 1 << " y0 = " << sloy_dop[0][2] << endl;
		cout << "iter = " << iter + 1 << " y1 = " << sloy_dop[0][3] << endl;
		cout << delta_u << endl;

		dop_mesh.r[1] = sloy_dop[0][0];
		dop_mesh.r[2] = sloy_dop[0][1];
		dop_mesh.z[1] = sloy_dop[0][2];
		dop_mesh.z[2] = sloy_dop[0][3];


		if (abs(J_prev - J_next) < eps) {
			break;
		}
		else {
			//tmp_u = tok;
			//tmp_u = sloy_dop[1][2];
			tmp_u = sloy_dop[0][3];


		}
	}
}

int main()
{
	ofstream output("q.txt");
	// step 1. read data
	cout << "Reading data" << endl;
	read_mesh();
	read_dop_mesh();
	sigma_read();
	
	// step 2. generate synthetic
	cout << "Generate synthetic" << endl;
	field_selection();

	// step 3. solving inverse task.
	cout << "Solving inverse task" << endl;
	inverse_problem();
}
