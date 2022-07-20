#include "GRF.hpp"
#include "FEM_system.hpp"
#include "FEM_solver.hpp"
#include "sampler.hpp"
#include "darcy_system_test.hpp"
#include <fstream>
#include "gtest/gtest.h"


using T = double;
const int K = 20;
double sigma = 1, jet = 1, l = 0.1;


//Tests for the FEM system and solver.
TEST(FEMTest, ZeroFieldToPerm) {
  //input a zero matrix into the field_to_perm routine and check that the resulting permeability is a matrix of ones
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  Darcy::System<T>::MT input;
  input = Darcy::System<T>::MT::Zero(K,K);
  Darcy::System<T>::MT expected;
  expected = Darcy::System<T>::MT::Ones(K,K);
  FEM::System<T> fsys(dsys,K);
  fsys.field_to_perm(input);
  //std::cout << fsys.get_PM() << std::endl << expected << std::endl;
  ASSERT_EQ(fsys.get_PM(),expected);
}

TEST(FEMTest, RandomFieldToPerm) {
  //input a randomly generated matrix (using eigen random) and check that field_to_perm provides the same result as by explicitly calling the elementwise exponential
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  Darcy::System<T>::MT input;
  input = Darcy::System<T>::MT::Random(K,K);
  Darcy::System<T>::MT expected;
  expected = exp(input.array()).matrix();
  FEM::System<T> fsys(dsys,K);
  fsys.field_to_perm(input);
  //std::cout << fsys.get_PM() << std::endl << expected << std::endl;
  ASSERT_EQ(fsys.get_PM(),expected);
}

TEST(FEMTest, SetupWithKnownField) {
  //Read in a known system matrix and rhs from a given matlab FEM routine
  //and assert equality to generated FEM system.
    const int n_int = (K - 1) * (K + 1);
    Darcy::System<T>::VTB boundary;
    boundary << 0, 1, 0, 1;
    Darcy_Test::System<T> dsys(boundary);
    typename Darcy::System<T>::MT input;
    input = Darcy::System<T>::MT::Ones(K + 1, K + 1);
    typename Darcy::System<T>::MT expected_A = Darcy::System<T>::MT::Zero(n_int, n_int);
    //Matrix from Matlab stored as a row Vector in A_K20.dat
    std::ifstream inA("../input_files/A_K20.dat");
    for (size_t i = 0; i < n_int; i++)
        for (size_t j = 0; j < n_int; j++)
            inA >> expected_A(i, j);
    inA.close();
    typename Darcy::System<T>::VT expected_b = Darcy::System<T>::VT::Zero(n_int);
    std::fstream inb("../input_files/b_K20.dat");
    for (size_t i = 0; i < n_int; i++)
        inb >> expected_b(i);
    inb.close();
    FEM::System<T> fsys(dsys,K+1);
    fsys.field_to_perm(input);
    fsys.setup();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Dense_A = fsys.get_A();
    //ASSERT_NEAR(Dense_A, expected_A, 0.00001);
    ASSERT_TRUE(expected_A.isApprox(Dense_A, 1e-5));
    //ASSERT_NEAR(fsys.get_b(), expected_b, 0.00001);
    ASSERT_TRUE(expected_b.isApprox(fsys.get_b(), 1e-4));
}

TEST(FEMTest, SolveWithKnownField) {
    //Same as SetupWithKnownField except also solve the system and compare results.
    //Test done on solution as Vector; for Matrix version uncomment lines.
    const int nvtx = (K + 1) * (K + 1);
    Darcy::System<T>::VTB boundary;
    boundary << 0, 1, 0, 1;
    Darcy_Test::System<T> dsys(boundary);
    typename Darcy::System<T>::MT input;
    input = Darcy::System<T>::MT::Ones(K + 1, K + 1);
    std::ifstream inSol("../input_files/sol_K20.dat");
    typename Darcy::System<T>::MT expected_sol;
    //expected_sol= Darcy::System<T>::MT::Zero(K + 1, K + 1);
    /*for (size_t i = 0; i < K + 1; i++)
        for (size_t j = 0; j < K + 1; j++)
            inSol >> expected_sol(i, j);*/
    expected_sol = Darcy::System<T>::MT::Zero(nvtx, 1);
    for (size_t i = 0; i < nvtx; i++)inSol >> expected_sol(i);
    inSol.close();
    FEM::System<T> fsys(dsys,K+1);
    fsys.field_to_perm(input);
    fsys.setup();
    FEM::Solver<T> fsol;
    fsol.solve(fsys);
    //typename Darcy::System<T>::MT SOL = fsys.get_p();
    //SOL.resize(K + 1, K + 1);
    //ASSERT_NEAR(SOL,expected_sol,0.00001)
    //ASSERT_TRUE(expected_sol.isApprox(SOL, 1e-5));
    //ASSERT_NEAR(fsys.get_p(), expected_sol, 0.00001);
    ASSERT_TRUE(expected_sol.isApprox(fsys.get_p(), 1e-5));
}

//Tests for the generator. Not sure how much can be done here, since it is not really a deterministic process unless the seeds are given.
TEST(GRFTest, GenerateFieldWithGivenSeed) {
  // Give the generator a seed and compare to an older realisation.
  std::fstream input("../input_files/field.data");
  int size; input >> size;
  GRF_matlib::MT<T> expected;
  expected = GRF_matlib::MT<T>::Zero(size,size);
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      input >> expected(i,j);
    }
  }
  input.close();
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  GRF::Generator<T> gen(size, dsys, true);
  gen.generate_field();
  GRF_matlib::MT<T> actual = gen.get_field();
  ASSERT_TRUE(expected.isApprox(actual, 1e-4));
}

TEST(GRFTest, FFT) {
  //fft the hilbert matrix and compare to matlabs fft2 routine. size of matrix is determined by the input file, which contains matlab result
  std::fstream input("../input_files/fft2.data");
  int size; input >> size;
  GRF_matlib::MTC<T> expected;
  expected = GRF_matlib::MTC<T>::Zero(size,size);
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      input >> expected(i,j);
    }
  }
  input.close();
  GRF_matlib::MTC<T> hilbert;
  GRF_matlib::MTC<T> actual;
  hilbert = GRF_matlib::MTC<T>::Zero(size,size);
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      hilbert(i,j) = static_cast<double>(1)/static_cast<double>(i+j+1);
    }
  }
  actual = GRF_matlib::fft(hilbert);
  //use Eigen's isApprox function here because there may be some small rounding
  //errors between our version and matlab, as well as sign differences on values
  //very close to 0
  ASSERT_TRUE(expected.isApprox(actual,1e-6));
}

TEST(GRFTest, IFFT) {
  //fft the hilbert matrix and compare to matlabs fft2 routine. size of matrix is determined by the input file, which contains matlab result
  std::fstream input("../input_files/ifft2.data");
  int size; input >> size;
  GRF_matlib::MT<std::complex<T>> expected;
  expected = GRF_matlib::MTC<T>::Zero(size,size);
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      input >> expected(i,j);
    }
  }
  input.close();
  GRF_matlib::MT<T> hilbert;
  GRF_matlib::MTC<T> actual;
  hilbert = GRF_matlib::MT<T>::Zero(size,size);
  for(int i = 0; i < size; i++) {
    for(int j = 0; j < size; j++) {
      hilbert(i,j) = static_cast<double>(1)/static_cast<double>(i+j+1);
    }
  }
  actual = GRF_matlib::ifft(hilbert);
  //use Eigen's isApprox function here because there may be some small rounding
  //errors between our version and matlab, as well as sign differences on values
  //very close to 0
  //lower threshold here for some reason that I don't understand
  ASSERT_TRUE(expected.isApprox(actual,1e-4));
}

TEST(GRFTest, FFTShiftEvenRowsEvenCols) {
  //call fft_shift on a matrix with even number of rows and columns and compare to expected shift
  //4x4 for simplicity here
  int size = 4;
  GRF_matlib::MT<T> matrix;
  matrix = GRF_matlib::MT<T>::Zero(size,size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      matrix(i,j) = i*size + j;
    }
  }
  //matrix is
  // 0  1  2  3
  // 4  5  6  7
  // 8  9 10 11
  //12 13 14 15
  GRF_matlib::MT<T> shifted_matrix = GRF_matlib::fft_shift(matrix);
  //shift by hand here, too lazy to figure out a good way to loop it (thats what fft_shift is for, anyway)
  matrix(0,0) = 10; matrix(0,1) = 11; matrix(0,2) = 8; matrix(0,3) = 9;
  matrix(1,0) = 14; matrix(1,1) = 15; matrix(1,2) = 12; matrix(1,3) = 13;
  matrix(2,0) = 2; matrix(2,1) = 3; matrix(2,2) = 0; matrix(2,3) = 1;
  matrix(3,0) = 6; matrix(3,1) = 7; matrix(3,2) = 4; matrix(3,3) = 5;
  ASSERT_EQ(matrix, shifted_matrix);
}

TEST(GRFTest, FFTShiftEvenRowsOddCols) {
  //call fft_shift on a matrix with even number of rows and odd number of columns and compare to expected shift
  int size = 4;
  GRF_matlib::MT<T> matrix;
  matrix = GRF_matlib::MT<T>::Zero(size,size+1);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size+1; j++) {
      matrix(i,j) = i*(size+1) + j;
    }
  }
  //matrix is
  // 0  1  2  3  4
  // 5  6  7  8  9
  //10 11 12 13 14
  //15 16 17 18 19
  GRF_matlib::MT<T> shifted_matrix = GRF_matlib::fft_shift(matrix);
  //shift by hand here, too lazy to figure out a good way to loop it (thats what fft_shift is for, anyway)
  matrix(0,0) = 13; matrix(0,1) = 14; matrix(0,2) = 10; matrix(0,3) = 11; matrix(0,4) = 12;
  matrix(1,0) = 18; matrix(1,1) = 19; matrix(1,2) = 15; matrix(1,3) = 16; matrix(1,4) = 17;
  matrix(2,0) = 3; matrix(2,1) = 4; matrix(2,2) = 0; matrix(2,3) = 1; matrix(2,4) = 2;
  matrix(3,0) = 8; matrix(3,1) = 9; matrix(3,2) = 5; matrix(3,3) = 6; matrix(3,4) = 7;
  ASSERT_EQ(matrix, shifted_matrix);
}

TEST(GRFTest, FFTShiftOddRowsEvenCols) {
  //call fft_shift on a matrix with odd number of rows and even number of columns and compare to expected shift
  int size = 4;
  GRF_matlib::MT<T> matrix;
  matrix = GRF_matlib::MT<T>::Zero(size+1,size);
  for (int i = 0; i < size+1; i++) {
    for (int j = 0; j < size; j++) {
      matrix(i,j) = i*(size) + j;
    }
  }
  //matrix is
  // 0  1  2  3
  // 4  5  6  7
  // 8  9 10 11
  //12 13 14 15
  //16 17 18 19

  GRF_matlib::MT<T> shifted_matrix = GRF_matlib::fft_shift(matrix);
  //shift by hand here, too lazy to figure out a good way to loop it (thats what fft_shift is for, anyway)
  matrix(0,0) = 14; matrix(0,1) = 15; matrix(0,2) = 12; matrix(0,3) = 13;
  matrix(1,0) = 18; matrix(1,1) = 19; matrix(1,2) = 16; matrix(1,3) = 17;
  matrix(2,0) = 2; matrix(2,1) = 3; matrix(2,2) = 0; matrix(2,3) = 1;
  matrix(3,0) = 6; matrix(3,1) = 7; matrix(3,2) = 4; matrix(3,3) = 5;
  matrix(4,0) = 10; matrix(4,1) = 11; matrix(4,2) = 8; matrix(4,3) = 9;
  ASSERT_EQ(matrix, shifted_matrix);
}

TEST(GRFTest, FFTShiftOddRowsOddCols) {
  //call fft_shift on a matrix with odd number of rows and columns and compare to expected shift
  int size = 4;
  GRF_matlib::MT<T> matrix;
  matrix = GRF_matlib::MT<T>::Zero(size+1,size+1);
  for (int i = 0; i < size+1; i++) {
    for (int j = 0; j < size+1; j++) {
      matrix(i,j) = i*(size+1) + j;
    }
  }
  //matrix is
  // 0  1  2  3  4
  // 5  6  7  8  9
  //10 11 12 13 14
  //15 16 17 18 19
  //20 21 22 23 24
  GRF_matlib::MT<T> shifted_matrix = GRF_matlib::fft_shift(matrix);
  //shift by hand here, too lazy to figure out a good way to loop it (thats what fft_shift is for, anyway)
  matrix(0,0) = 18; matrix(0,1) = 19; matrix(0,2) = 15; matrix(0,3) = 16; matrix(0,4) = 17;
  matrix(1,0) = 23; matrix(1,1) = 24; matrix(1,2) = 20; matrix(1,3) = 21; matrix(1,4) = 22;
  matrix(2,0) = 3; matrix(2,1) = 4; matrix(2,2) = 0; matrix(2,3) = 1; matrix(2,4) = 2;
  matrix(3,0) = 8; matrix(3,1) = 9; matrix(3,2) = 5; matrix(3,3) = 6; matrix(3,4) = 7;
  matrix(4,0) = 13; matrix(4,1) = 14; matrix(4,2) = 10; matrix(4,3) = 11; matrix(4,4) = 12;
  ASSERT_EQ(matrix, shifted_matrix);
}

TEST(GRFTest, ReducedCov) {
  //Pass in a darcy system with simple covariance function with small number of sample points.
  //Create the reduced covariance matrix and compare it to the expected matrix
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  int N = 3; T dx = 1.0/(N-1);
  //no padding used here, just for demo purposes, hardcoding the expected matrix
  GRF_matlib::MT<T> actual = GRF_matlib::reduced_cov<T>(N,N,dx,dx,&dsys);
  GRF_matlib::MT<T> expected(2*N-1,2*N-1);
  //fill expected from dsys, hardcoded to compare against what reduced_cov does
  Darcy::System<T>::VTX x;
  x(0)=-1;   x(1)=-1;   expected(0,0)=dsys.c(x);
  x(0)=-0.5;            expected(1,0)=dsys.c(x);
  x(0)=0;               expected(2,0)=dsys.c(x);
  x(0)=0.5;             expected(3,0)=dsys.c(x);
  x(0)=1;               expected(4,0)=dsys.c(x);
  x(0)=-1;   x(1)=-0.5; expected(0,1)=dsys.c(x);
  x(0)=-0.5;            expected(1,1)=dsys.c(x);
  x(0)=0;               expected(2,1)=dsys.c(x);
  x(0)=0.5;             expected(3,1)=dsys.c(x);
  x(0)=1;               expected(4,1)=dsys.c(x);
  x(0)=-1;   x(1)=0;    expected(0,2)=dsys.c(x);
  x(0)=-0.5;            expected(1,2)=dsys.c(x);
  x(0)=0;               expected(2,2)=dsys.c(x);
  x(0)=0.5;             expected(3,2)=dsys.c(x);
  x(0)=1;               expected(4,2)=dsys.c(x);
  x(0)=-1;   x(1)=0.5;  expected(0,3)=dsys.c(x);
  x(0)=-0.5;            expected(1,3)=dsys.c(x);
  x(0)=0;               expected(2,3)=dsys.c(x);
  x(0)=0.5;             expected(3,3)=dsys.c(x);
  x(0)=1;               expected(4,3)=dsys.c(x);
  x(0)=-1;   x(1)=1;    expected(0,4)=dsys.c(x);
  x(0)=-0.5;            expected(1,4)=dsys.c(x);
  x(0)=0;               expected(2,4)=dsys.c(x);
  x(0)=0.5;             expected(3,4)=dsys.c(x);
  x(0)=1;               expected(4,4)=dsys.c(x);
  ASSERT_TRUE(expected.isApprox(actual,1e-6));
}


//Tests for Sampler
TEST(SamplerTest, MeanAndVarianceOfOneSample) {
  //take one monte carlo sample and compare the mean and variance to the result of the solved system (should be equal resp. zero)
  //Use EXPECT_EQ instead of ASSERT_EQ so execution continues if theres a problem
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  FEM::System<T> fsys(dsys,K);
  FEM::Solver<T> fsol;
  GRF::Generator<T> gen(K+1,dsys);
  UQ::Sampler<T> sampler(1,gen,fsys,fsol,K);
  sampler.sample();
  Darcy::System<T>::MT result = fsys.get_p();
  result.resize(K+1,K+1);
  EXPECT_EQ(result,sampler.get_mean());
  Darcy::System<T>::MT expected_variance;
  expected_variance = Darcy::System<T>::MT::Zero(K+1,K+1);
  EXPECT_EQ(expected_variance,sampler.get_variance());
}

TEST(SamplerTest, MeanAndVarianceOnDirichletBoundary) {
  //take 10 samples (could be any number) and check the mean and variance on the dirichlet boundary. Should be equal to given pD function resp zero.
  //Use EXPECT_EQ
  Darcy::System<T>::VTB boundary;
  boundary << 0,1,0,1;
  Darcy_Test::System<T> dsys(boundary);
  FEM::System<T> fsys(dsys,K);
  FEM::Solver<T> fsol;
  GRF::Generator<T> gen(K+1,dsys);
  UQ::Sampler<T> sampler(10,gen,fsys,fsol,K);
  sampler.sample();
  Darcy::System<T>::VT actual_left_boundary, actual_right_boundary, expected_left_boundary, expected_right_boundary;
  actual_left_boundary = Darcy::System<T>::VT::Zero(K+1);
  actual_right_boundary = Darcy::System<T>::VT::Zero(K+1);
  expected_left_boundary = Darcy::System<T>::VT::Zero(K+1);
  expected_right_boundary = Darcy::System<T>::VT::Zero(K+1);
  Darcy::System<T>::VTX x;
  //fill actual variances, the expected ones are zero anyway
  for (int i = 0; i <= K; i++) {
    actual_left_boundary(i) = sampler.get_variance()(i,0);
    actual_right_boundary(i) = sampler.get_variance()(i,K);
  }
  EXPECT_EQ(actual_left_boundary, expected_left_boundary);
  EXPECT_EQ(actual_right_boundary, expected_right_boundary);
  //now we use the vectors for the mean
  //fill expected and actual mean values
  for (int i = 0; i <= K; i++) {
    x(0)=0;
    x(1)=static_cast<double>(i)/static_cast<double>(K);
    expected_left_boundary(i) = dsys.pD(x);
    x(0)=1;
    expected_right_boundary(i) = dsys.pD(x);
    actual_left_boundary(i) = sampler.get_mean()(i,0);
    actual_right_boundary(i) = sampler.get_mean()(i,K);
  }
  EXPECT_EQ(actual_left_boundary, expected_left_boundary);
  EXPECT_EQ(actual_right_boundary, expected_right_boundary);

}

int main(int argc, char* argv[]) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
