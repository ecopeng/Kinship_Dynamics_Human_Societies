/*
THE HYBRID VERSION OF THE NEWLY DEVELOPED MODEL OF KINSHIP DYNAMICS
*/

#include <iostream>
#include <random>
#include <vector>
#include <numeric>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <Eigen/Dense>
#include <mpi.h>

std::vector<int> SAMPLE(int n, const std::vector<int>& s, const std::vector<double>& p) {
  std::random_device rd;
  std::mt19937 generator(rd());
  std::discrete_distribution<int> dist(p.begin(), p.end());
  std::vector<int> sample;
  for(int i = 0; i < n; ++i) {sample.push_back(s[dist(generator)]);}
  return sample;
}

int SAMPLE(double p) {
  std::vector<double> pq = {1 - p, p};
  std::random_device rd;
  std::mt19937 generator(rd());
  std::discrete_distribution<int> dist(pq.begin(), pq.end());
  return dist(generator);
}

std::vector<int> FATE(const std::vector<int>& ac, const std::vector<double>& SUR) {
  std::vector<int> vec;
  for(int i = 0; i < ac.size(); ++i) {vec.push_back(SAMPLE(SUR[ac[i] - 1]));}
  return vec;
}

std::vector<double> READ_VECTOR_DOU(const std::string& path_fn) {
  std::ifstream i_file(path_fn);
  if(!i_file.is_open()) {std::cout << "\n/* --- cannot open file! --- */\n"; std::exit(EXIT_FAILURE);}
  std::string str;
  std::vector<double> vec;
  while(std::getline(i_file, str)) {vec.push_back(std::stod(str));}
  return vec;
}

std::vector<double> STABLE_AC_FREQ(const std::string& path_fn) {
  std::ifstream i_file(path_fn);
  if(!i_file.is_open()) {std::cout << "\n/* --- cannot open file! --- */\n"; std::exit(EXIT_FAILURE);}
  std::string str;
  std::vector<double> v;
  while(std::getline(i_file, str)) {v.push_back(std::stod(str));}
  double f = 0;
  std::for_each(v.begin(), v.end(), [&f](double i) {f += i;});
  std::for_each(v.begin(), v.end(), [f](double& i) {i /= f;});
  return v;
}

std::vector<double> v_L_Ka(const std::vector<double>& SUR_F, const std::vector<double>& FEC_F, const std::vector<int>& ac_f, double dispersal, const std::vector<double>& AC_FREQ_F) {
  int t; double exp_fec = 0, fec_sum = 0;
  for(t = 0; t < AC_FREQ_F.size(); ++t) {exp_fec += AC_FREQ_F[t] * FEC_F[t];}
  std::for_each(ac_f.begin(), ac_f.end(), [&fec_sum, FEC_F](int i) {fec_sum += FEC_F[i - 1];});
  double de = (1 - dispersal) * fec_sum + dispersal * ac_f.size() * exp_fec;
  if(de != 0) {
    std::vector<double> L;
    for(t = 0; t < ac_f.size(); ++t) {L.push_back((1 - dispersal) * FEC_F[ac_f[t] - 1] / de);}
    return L;
  } else {
    std::vector<double> L;
    for(t = 0; t < ac_f.size(); ++t) {L.push_back(0);}
    return L;
  }
}

std::vector<double> v_L_Kb(const std::vector<double>& SUR_M, const std::vector<double>& FEC_M, const std::vector<int>& ac_m, double mating, const std::vector<double>& AC_FREQ_M) {
  int t; double exp_fec = 0, fec_sum = 0;
  for(t = 0; t < AC_FREQ_M.size(); ++t) {exp_fec += AC_FREQ_M[t] * FEC_M[t];}
  std::for_each(ac_m.begin(), ac_m.end(), [&fec_sum, FEC_M](int i) {fec_sum += FEC_M[i - 1];});
  double de = mating * fec_sum + (1 - mating) * ac_m.size() * exp_fec;
  if(de != 0) {
    std::vector<double> L;
    for(t = 0; t < ac_m.size(); ++t) {L.push_back(mating * FEC_M[ac_m[t] - 1] / de);}
    return L;
  } else {
    std::vector<double> L;
    for(t = 0; t < ac_m.size(); ++t) {L.push_back(0);}
    return L;
  }
}

void AC_UPDATE(std::vector<int>& ac, const std::vector<int>& fate, int max_ac) {
  std::transform(ac.begin(), ac.end(), fate.begin(), ac.begin(), [max_ac](int a, int b) {return b == 0 ? 1 : (a + b <= max_ac ? a + b : 1);});
}

Eigen::MatrixXd R_UPDATE(const std::vector<int>& fate_f, const std::vector<int>& fate_m, const std::vector<double>& p_FM2f, const std::vector<double>& p_FM2m, const Eigen::MatrixXd& R_0) {
  Eigen::MatrixXd R_1 = Eigen::MatrixXd::Identity(R_0.rows(), R_0.cols());
  int i, j, r, c;
  for(i = 0; i < R_1.rows(); ++i) {
    for(j = i + 1; j < R_1.cols(); ++j) {
      if(i < fate_f.size() && j < fate_f.size()) { // F-F
        if(fate_f[i] == 0 && fate_f[j] == 0) { /* new F to new F */
          for(r = 0; r < R_0.rows(); ++r) {
            for(c = 0; c < R_0.cols(); ++c) {
              R_1(i, j) += p_FM2f[r] * p_FM2f[c] * R_0(r, c);
              R_1(j, i) += p_FM2f[r] * p_FM2f[c] * R_0(r, c);
            }
          }
        } else if(fate_f[i] == 0 && fate_f[j] == 1) { /* new F to old F */
          for(r = 0; r < R_0.rows(); ++r) {
            R_1(i, j) += p_FM2f[r] * R_0(r, j);
            R_1(j, i) += p_FM2f[r] * R_0(r, j);
          }
        } else if(fate_f[i] == 1 && fate_f[j] == 0) { /* old F to new F */
          for(c = 0; c < R_0.cols(); ++c) {
            R_1(i, j) += p_FM2f[c] * R_0(i, c);
            R_1(j, i) += p_FM2f[c] * R_0(i, c);
          }
        } else { /* old F to old F */
          R_1(i, j) = R_0(i, j);
          R_1(j, i) = R_0(i, j);
        }
      } else if(i < fate_f.size() && j >= fate_f.size()) { // F-M
        if(fate_f[i] == 0 && fate_m[j - fate_f.size()] == 0) { /* new F to new M */
          for(r = 0; r < R_0.rows(); ++r) {
            for(c = 0; c < R_0.cols(); ++c) {
              R_1(i, j) += p_FM2f[r] * p_FM2m[c] * R_0(r, c);
              R_1(j, i) += p_FM2f[r] * p_FM2m[c] * R_0(r, c);
            }
          }
        } else if(fate_f[i] == 0 && fate_m[j - fate_f.size()] == 1) { /* new F to old M */
          for(r = 0; r < R_0.rows(); ++r) {
            R_1(i, j) += p_FM2f[r] * R_0(r, j);
            R_1(j, i) += p_FM2f[r] * R_0(r, j);
          }
        } else if(fate_f[i] == 1 && fate_m[j - fate_f.size()] == 0) { /* old F to new M */
          for(c = 0; c < R_0.cols(); ++c) {
            R_1(i, j) += p_FM2m[c] * R_0(i, c);
            R_1(j, i) += p_FM2m[c] * R_0(i, c);
          }
        } else { /* old F to old M */
          R_1(i, j) = R_0(i, j);
          R_1(j, i) = R_0(i, j);
        }
      } else { // M-M
        if(fate_m[i - fate_f.size()] == 0 && fate_m[j - fate_f.size()] == 0) { /* new M to new M */
          for(r = 0; r < R_0.rows(); ++r) {
            for(c = 0; c < R_0.cols(); ++c) {
              R_1(i, j) += p_FM2m[r] * p_FM2m[c] * R_0(r, c);
              R_1(j, i) += p_FM2m[r] * p_FM2m[c] * R_0(r, c);
            }
          }
        } else if(fate_m[i - fate_f.size()] == 0 && fate_m[j - fate_f.size()] == 1) { /* new M to old M */
          for(r = 0; r < R_0.rows(); ++r) {
            R_1(i, j) += p_FM2m[r] * R_0(r, j);
            R_1(j, i) += p_FM2m[r] * R_0(r, j);
          }
        } else if(fate_m[i - fate_f.size()] == 1 && fate_m[j - fate_f.size()] == 0) { /* old M to new M */
          for(c = 0; c < R_0.rows(); ++c) {
            R_1(i, j) += p_FM2m[c] * R_0(i, c);
            R_1(j, i) += p_FM2m[c] * R_0(i, c);
          }
        } else { /* old M to old M */
          R_1(i, j) = R_0(i, j);
          R_1(j, i) = R_0(i, j);
        }
      }
    }
  }
  return R_1(Eigen::all, Eigen::all);
}

void OS_SEX_AC_AVG_R(int t, std::ostream& out, const Eigen::MatrixXd& R, const std::vector<int>& ac_f, const std::vector<int>& ac_m) {
  double r2F, r2M, r2FM;
  for(int r = 0; r < R.rows(); ++r) {
    if(r < ac_f.size()) { // a focal female
      r2F = (R(r, Eigen::seq(0, ac_f.size() - 1)).sum() - 1) / (ac_f.size() - 1); // f2F
      r2M = R(r, Eigen::seq(ac_f.size(), ac_f.size() + ac_m.size() - 1)).sum() / ac_m.size(); // f2M
      r2FM = (R.row(r).sum() - 1) / (R.cols() - 1);
      out << t << "\t" << "F" << "\t" << ac_f[r] << "\t" << r2F << "\t" << r2M << "\t" << r2FM << "\n";
    } else { // a focal male
      r2F = R(r, Eigen::seq(0, ac_f.size() - 1)).sum() / ac_f.size(); // m2F
      r2M = (R(r, Eigen::seq(ac_f.size(), ac_f.size() + ac_m.size() - 1)).sum() - 1) / (ac_m.size() - 1); // m2M
      r2FM = (R.row(r).sum() - 1) / (R.cols() - 1);
      out << t << "\t" << "M" << "\t" << ac_m[r - ac_f.size()] << "\t" << r2F << "\t" << r2M << "\t" << r2FM << "\n";
    }
  }
}

void SIMULATE(std::ostream& os, int T,  int cut, double df, double dm, double rho, int Nf, int Nm, 
  const std::vector<double>& SUR_F, const std::vector<double>& SUR_M,
  const std::vector<double>& FEC_F, const std::vector<double>& FEC_M,
  const std::vector<double>& lx_F, const std::vector<double>& lx_M,
  const std::vector<double>& AC_FREQ_F, const std::vector<double>& AC_FREQ_M,
  const std::vector<int>& AC_F, const std::vector<int>& AC_M, int max_ac_f, int max_ac_m) {
  /***** INITIALIZE AGE (OR AGE CLASSES) - DRAW FROM STABLE AGE CLASS DISTRIBUTIONS *****/
  std::vector<int> ac_f = SAMPLE(Nf, AC_F, lx_F);
  std::vector<int> ac_m = SAMPLE(Nm, AC_M, lx_M);
  std::vector<int> fate_f(Nf), fate_m(Nm);
  std::vector<double> v_L_Ka2a(Nf), v_L_Ka2b(Nf), v_L_Kb2(Nm);
  std::vector<double> v_H_Ka2a(Nf), v_H_Ka2b(Nf), v_H_Kb2a(Nm), v_H_Kb2b(Nm);
  /***** INITIALIZE RELATEDNESS MATRIX *****/
  Eigen::MatrixXd R = Eigen::MatrixXd::Identity(Nf + Nm, Nf + Nm);
  double p_native_f, p_native_m;
  std::vector<double> p_FM2f, p_FM2m;
  std::string o_dir = std::filesystem::current_path().string() + std::string("/out");
  if(!std::filesystem::is_directory(o_dir)) {std::filesystem::create_directory(o_dir);}
  for(int t = 0; t < T; ++t) {
    /***** CALCULATE THE PROBABILITIES OF INHERITANCE *****/
    v_L_Ka2a = v_L_Ka(SUR_F, FEC_F, ac_f, df, AC_FREQ_F);
    v_L_Ka2b = v_L_Ka(SUR_F, FEC_F, ac_f, dm, AC_FREQ_F);
    v_L_Kb2 = v_L_Kb(SUR_M, FEC_M, ac_m, rho, AC_FREQ_M);
    p_native_f = 0;
    p_native_m = 0;
    std::for_each(v_L_Ka2a.begin(), v_L_Ka2a.end(), [&p_native_f](double i) {p_native_f += i;});
    std::for_each(v_L_Ka2b.begin(), v_L_Ka2b.end(), [&p_native_m](double i) {p_native_m += i;});
    v_H_Ka2a = v_L_Ka2a;
    std::for_each(v_H_Ka2a.begin(), v_H_Ka2a.end(), [](double& i) {i *= .5;});
    v_H_Ka2b = v_L_Ka2b;
    std::for_each(v_H_Ka2b.begin(), v_H_Ka2b.end(), [](double& i) {i *= .5;});
    v_H_Kb2a = v_L_Kb2;
    std::for_each(v_H_Kb2a.begin(), v_H_Kb2a.end(), [p_native_f](double& i) {i *= (.5 * p_native_f);});
    v_H_Kb2b = v_L_Kb2;
    std::for_each(v_H_Kb2b.begin(), v_H_Kb2b.end(), [p_native_m](double& i) {i *= (.5 * p_native_m);});
    p_FM2f = v_H_Ka2a;
    p_FM2f.insert(p_FM2f.end(), v_H_Kb2a.begin(), v_H_Kb2a.end());
    p_FM2m = v_H_Ka2b;
    p_FM2m.insert(p_FM2m.end(), v_H_Kb2b.begin(), v_H_Kb2b.end());
    /***** SIMULATE THE FATES OF INDIVIDUALS *****/
    fate_f = FATE(ac_f, SUR_F);
    fate_m = FATE(ac_m, SUR_M);
    /***** CALCULATE INTER-INDIVIDUAL RELATEDNESS *****/
    R = R_UPDATE(fate_f, fate_m, p_FM2f, p_FM2m, R);
    if(t >= cut) {OS_SEX_AC_AVG_R(t, os, R, ac_f, ac_m);}
    AC_UPDATE(ac_f, fate_f, max_ac_f);
    AC_UPDATE(ac_m, fate_m, max_ac_m);
  }
}

int main(int argc, char* argv[]) {
  std::vector<double> SUR_F = READ_VECTOR_DOU("./par/SUR_F");
  std::vector<double> SUR_M = READ_VECTOR_DOU("./par/SUR_M");
  std::vector<double> FEC_F = READ_VECTOR_DOU("./par/FEC_F");
  std::vector<double> FEC_M = READ_VECTOR_DOU("./par/FEC_M");
  std::vector<double> lx_F = READ_VECTOR_DOU("./par/lx_F");
  std::vector<double> lx_M = READ_VECTOR_DOU("./par/lx_M");
  std::vector<double> AC_FREQ_F = STABLE_AC_FREQ("./par/lx_F"), AC_FREQ_M = STABLE_AC_FREQ("./par/lx_M");
  std::vector<int> AC_F(SUR_F.size()), AC_M(SUR_M.size());
  std::iota(std::begin(AC_F), std::end(AC_F), 1);
  std::iota(std::begin(AC_M), std::end(AC_M), 1);
  int max_ac_f = AC_F.size(), max_ac_m = AC_M.size(), Nf, Nm, T = 1e+06, cut = 9.9e+05;
  std::ofstream r;
  double df, dm, rho;

  int n_par_set = std::stoi(argv[1]);
  std::vector<double> NDM; /* NUMBERS OF FEMALES & MALES, FEMALE & MALE DISPERSAL RATES, RATE OF LOCAL MATING */

  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::string path_out = std::filesystem::current_path().root_path().string() + "out/";
  if(!std::filesystem::is_directory(path_out)) {
    path_out = std::filesystem::current_path().string() + "/out/";
    std::filesystem::create_directory(path_out);
  }

  for(int q = 0; q < n_par_set; ++q) {
    NDM = READ_VECTOR_DOU("./par/NDM_" + std::to_string(q));
    Nf = NDM[0], Nm = NDM[1], df = NDM[2], dm = NDM[3], rho = NDM[4];
    r.open("./out/r_____rep_" + std::to_string(rank) + "___" + std::to_string(Nf) + "+" + std::to_string(Nm) + "___" + std::to_string(df) + "+" + std::to_string(dm) + "___" + std::to_string(rho));
    r << "t\tsex\tage\tr2F\tr2M\tr2FM\n";
    SIMULATE(r, T, cut, df, dm, rho, Nf, Nm, SUR_F, SUR_M, FEC_F, FEC_M, lx_F, lx_M, AC_FREQ_F, AC_FREQ_M, AC_F, AC_M, max_ac_f, max_ac_m);
    r.close();
  }

  MPI_Finalize();

  return 0;
}