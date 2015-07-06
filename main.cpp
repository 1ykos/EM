#include <algorithm>
#include <array>
#include <cfloat>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <tgmath.h>
#include <unordered_map>
#include <utility>
#include <vector>
#include <functional>

#define N 2 // number of gaussians in the mixture model

using std::abs;
using std::array;
using std::bernoulli_distribution; 
using std::cin;
using std::copy;
using std::cout;
using std::endl;
using std::get;
using std::getline;
using std::ifstream;
using std::istream;
using std::knuth_b;
using std::list;
using std::lower_bound;
using std::lognormal_distribution;
using std::map;
using std::max;
using std::min;
using std::mt19937;
using std::multimap;
using std::normal_distribution;
using std::numeric_limits;
using std::ofstream;
using std::pair;
using std::random_device;
using std::random_shuffle;
using std::reverse_iterator;
using std::round;
using std::setw;
using std::shuffle;
using std::sort;
using std::sqrt;
using std::string;
using std::stringstream;
using std::swap;
using std::tuple;
using std::uniform_real_distribution;
using std::uniform_int_distribution;
using std::unique;
using std::unordered_map;
using std::vector;

const double pi = M_PI;
random_device urandom;
mt19937 urng(urandom());

template<typename T> inline const void kahan_sum(
  const T &x,
  T &c,
  T &sum
){
  const T y=x-c;
  const T t=sum+y;
  c = (t-sum)-y;
  sum=t;
}

template<typename T> inline const void mean_variance(
  const T &x,
  const T &w,
  T &sumw,
  T &mean,
  T &M2
){
  const T temp = w+sumw;
  const T delta = x-mean;
  const T R = delta*w/temp;
  mean += R;
  M2 += sumw*delta*R;
  sumw = temp;
}

template<typename T> inline const T variance(
  const T &M2,
  const T &sumw,
  const size_t &n
){
  return M2*n/(sumw*(n-1));
}

template<typename T> inline const T variance(
  const T &M2,
  const T &sumw
){
  return M2/sumw;
}

const vector<array<double,3>> ex_max(
  const vector<array<double,2>> &data,
  const size_t &n
){
  double mean=0;
  double sumw=0;
  double M2=0;
  for (auto it=data.begin(); it!=data.end(); ++it){
    mean_variance((*it)[0],(*it)[1],sumw,mean,M2);
  }
  const double s = sqrt(variance(M2,sumw,data.size()));
  vector<array<double,3>> model(n); // weight,mean,sigma
  // initialize the model with equal weights,
  // random points from the dataset as mean
  // and the square root of the over all variance
  for (size_t i=0; i!=n; ++i){
    model[i]={{1.0/n,data[i][0],s/n}};
  }
  vector<array<double,5>> temp(n); // sumw,mean,M2,T,TS
  double change;
  do{
    change=0.0;
    double TS = 0.0;
    for (auto t=temp.begin(); t!=temp.end(); ++t){
      (*t)[4]=0.0;
    }
    for (auto d=data.begin(); d!=data.end(); ++d){
      auto t=temp.begin();
      auto m=model.begin();
      for (;m!=model.end()&&t!=temp.end(); (++m,++t)){
        const double s2=(*m)[2]*(*m)[2]+(*d)[1]*(*d)[1];
        const double i =
          exp(-((*d)[0]-(*m)[1])*((*d)[0]-(*m)[1])/(2.0*s2))/sqrt(2*s2*pi);
        TS+=i;
        (*t)[4]+=i;
      }  
    }  
    for (auto t=temp.begin(); t!=temp.end(); ++t){
      (*t)[0]=(*t)[1]=(*t)[2]=(*t)[3]=0.0;
    }
    {
      auto t=temp.begin();
      auto m=model.begin();
      for (;m!=model.end()&&t!=temp.end(); (++m,++t)){
        change+=((*m)[0]-(*t)[4]/TS)*((*m)[0]-(*t)[4]/TS);
        (*m)[0]=(*t)[4]/TS;
      }
    }
    for (auto d=data.begin(); d!=data.end(); ++d){
      double T = 0.0;
      auto t=temp.begin();
      auto m=model.begin();
      for (;m!=model.end()&&t!=temp.end(); (++m,++t)){
        const double s2=(*m)[2]*(*m)[2]+(*d)[1]*(*d)[1];
        const double i =
          exp(-((*d)[0]-(*m)[1])*((*d)[0]-(*m)[1])/(2.0*s2))/sqrt(2*s2*pi);
        if (isnan(i)){
          cout << (*d)[0] << " "
               << (*m)[1] << " "
               << (*m)[2] << " "
               << (*d)[1] << endl;
          cout << "NAAAAN !!! !!! " << endl;
          exit(0);
        }
        T+=i;
        (*t)[3]=i;
      }
      t=temp.begin();
      m=model.begin();
      if (T<numeric_limits<double>::epsilon()) continue;
      for (;m!=model.end()&&t!=temp.end(); (++m,++t)){
        if ((*t)[3]<numeric_limits<double>::epsilon()) continue;
        const double x = (*d)[0];
        const double w = (*t)[3]*(*m)[0]/(T*(*d)[1]*(*d)[1]);
        mean_variance(x,w,(*t)[0],(*t)[1],(*t)[2]);
      }
    }
    {
      auto t=temp.begin();
      auto m=model.begin();
      for (;m!=model.end()&&t!=temp.end(); (++m,++t)){
        change+=((*m)[0]-(*t)[4]/TS)*((*m)[0]-(*t)[4]/TS);
        change+=((*m)[1]-(*t)[1])*((*m)[1]-(*t)[1]);
        change+=((*m)[2]-variance((*t)[2],(*t)[0]))
               *((*m)[2]-variance((*t)[2],(*t)[0]));
        (*m)[0]=(*t)[4]/TS;
        (*m)[1]=(*t)[1];
        (*m)[2]=variance((*t)[2],(*t)[0]);
      }
    }
    cout << change << endl;
  }while(change>1e-24);
  return model;
}
                                     

int main(int argc, char *argv[]){
  lognormal_distribution<double> lognormal;
  uniform_int_distribution<size_t> uniform(0,N-1);
  vector<array<double,2>> data;
  for (size_t i=0; i<256*N; ++i){
    const double s = lognormal(urng);
    data.push_back({normal_distribution<double>(uniform(urng),s)(urng),s});
//     data.push_back({normal_distribution<double>(uniform(urng),0.5)(urng),0.1});
  }
//   for (auto it=data.begin(); it!=data.end(); ++it){
//     cout << setw(16) << (*it)[0] << setw(16) << (*it)[1] << endl;
//   }
//   cout << " * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
  vector<array<double,3>> model = ex_max(data,N);
  cout << " * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
  cout << setw(16) << "weight" << setw(16) << "mean" << setw(16) << "sigma" << endl; 
  sort(model.begin(),model.end(),
       [](const array<double,3> &a, const array<double,3> &b){return a[1]<b[1];});
  for (auto it=model.begin(); it!=model.end(); ++it){
    cout << setw(16) << (*it)[0] << setw(16) << (*it)[1] << setw(16) << (*it)[2] << endl;    
  }
  return 0;  
}
