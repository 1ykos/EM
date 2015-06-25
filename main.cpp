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
using std::unique;
using std::unordered_map;
using std::vector;

const double pi = M_PI;

template<class T>
inline constexpr T pow(const T &base, unsigned const &exponent)
{
    // (parentheses not required in next line)
    return (exponent == 0)     ? 1 :
           (exponent % 2 == 0) ? pow(base, exponent/2)*pow(base, exponent/2) :
           base * pow(base, (exponent-1)/2) * pow(base, (exponent-1)/2);
}

static const double gaussian(const double& x, const double& m, const double s){
  return exp(-(x-m)*(x-m)/(2*s*s))/(s*sqrt(2*pi));
}
                        //weight,μ     ,σ
static const vector<tuple<double,double,double>> ex_max_gauss(
  const vector<tuple<double,double>> &data, n){
  vector<tuple<double,double,double>> distr0;
  distr0.reserve(min(n,data.size()));
  for(auto it=data.begin();it!=data.end()&&distr0.size()<n;++it){
    distr0.push_back(1.0,get<0>(*it),get<1>(*it));
  }
  for (auto it=distr0.begin();it!=distr0.end();++it){
    get<0>(*it)=1.0/distr0.size();
  }
  vector<tuple<double,double,double>> distr1=distr0;
  for (double change=numeric_limits<double>::infinity();change>1e-10;){
    //distr0 old, distr1 new
    distr0=distr1;
    for (auto it0=distr0.begin(),it1=distr1.begin();
           it0!=distr0.end()&&it1!=distr1.end(); (++it0,++it1)){
      douleb t_n=0;
      for (auto it0=distr0.begin(); it0!=distr0.begin(); ++it0){
        for (auto it=data.begin(); it!=data.end(); ++it){
          const double m1=get<0>(*it);
          const double s1=get<1>(*it);
          const double t =get<0>(*it0);
          const double m2=get<1>(*it0);
          const double s2=get<2>(*it0);
//           const double m = (m1/(s1*s1)+m2/(s2*s2))
//           /(1.0/(s1*s1)+1.0/(s2*s2));
//           const double s = s1*s1*s2*s2/(s1*s1+s2*s2);
//           integrate(gaussian(m1,m2,sqrt(s1**2+s2**2))*gaussian(x,m,s))
          t_n+=t*gaussian(m1,m2,sqrt(s1*s1+s2*s2));
        }
      }
      for (auto it0=distr0.begin(),it1=distr1.begin();
           it0!=distr0.end()&&it1!=distr1.end(); (++it0,++it1)){
        get<0>(*it1)=0;
        for (auto it=data.begin(); it!=data.end(); ++it){
          const double m1=get<0>(*it);
          const double s1=get<1>(*it);
          const double t =get<0>(*it0);
          const double m2=get<1>(*it0);
          const double s2=get<2>(*it0);
          get<0>(*it1)+=t*gaussian(m1,m2,sqrt(s1*s1+s2*s2));
        }
      }
      for (auto it0=distr0.begin(); it0!=distr0.begin(); ++it0){
        
        for (auto it=data.begin(); it!=data.end(); ++it){
          const double m1=get<0>(*it);
          const double s1=get<1>(*it);
          const double m2=get<1>(*it0);
          const double s2=get<2>(*it0);
//           const double m = (m1/(s1*s1)+m2/(s2*s2))
//           /(1.0/(s1*s1)+1.0/(s2*s2));
//           const double s = s1*s1*s2*s2/(s1*s1+s2*s2);
//           integrate(gaussian(m1,m2,sqrt(s1**2+s2**2))*gaussian(x,m,s))
          get<0>(*it0)=get<0>gaussian(m1,m2,sqrt(s1*s1+s2*s2));
        }
      }
      
    }
    for (auto it=data.begin(); it!=data.end(); ++it){
      for (auto it0=distr0.begin(),it1=distr1.begin();
           it0!=distr0.end()&&it1!=distr1.end(); (++it0,++it1)){
      }
      get<0>(*it1)=new wheight;
      get<1>(*it1)=new mean;
      for (auto it0=distr0.begin(),it1=distr1.begin();
           it0!=distr0.end()&&it1!=distr1.end(); (++it0,++it1)){
        get<2>(it1)=
        
      }
    }
  }
}

void generate_distributions(
  const unordered_map<size_t,vector<tuple<size_t,double,double>>> &uniques,
  unordered_map<size_t,PPD> &distr,
  const column_vector &starting_point
){
  cout << "generating partiality distributions" << endl;
  for (auto it0=uniques.begin(); it0!=uniques.end(); ++it0){
    double mean=0.0;
    for (auto it1=it0->second.begin();it1!=it0->second.end();++it1){
      mean += get<1>(*it1);
    }
    //initialization
    mean/=it0->second.size();
    distr[it0->first]={0.0,mean,4*mean,4*mean,0.75};
    //maximization
    for (size_t i=0; i<8; ++i){
      const double n = it0->second.size();
      double s_t1 =0;
      double s_t2 =0;
      double m1_n =0;
      double m2_n =0;
      for (auto it1=it0->second.begin();it1!=it0->second.end();++it1){
        const double t1=distr[it0->first].t1;
        const double t2=1-t1;
        const double s =starting_point(get<0>(*it1));
        const double x =get<1>(*it1)*s;
        //         cout << x << endl;
        const double m1=distr[it0->first].m1;
        const double m2=distr[it0->first].m2;
        const double s1=distr[it0->first].s1;
        const double s2=distr[it0->first].s2;
        const double e1=(s1*s1*128<(x-m1)*(x-m1))?t1*0.00311664:t1*exp(-(x-m1)*(x-m1)/(2*s1*s1))/(s1*sqrt(2*pi));
        const double e2=(s2*s2*128<(x-m2)*(x-m2))?t2*0.00311664:t2*exp(-(x-m2)*(x-m2)/(2*s2*s2))/(s2*sqrt(2*pi));
        //         const double t1i =(abs(e1)+abs(e2)<1e-10)?(e1>e2?1:0):e1/(e1+e2);
        //         const double t2i =(abs(e1)+abs(e2)<1e-10)?(e2>e1?1:0):e2/(e1+e2);
        const double t1i =e1/(e1+e2);
        const double t2i =e2/(e1+e2);
        //         cout << x << " " << m1 << " " << m2 << " " << s1 << " " << s2 << " " << e1 << " " << e2 << " " << t1i << " " << t2i << endl;
        s_t1 += t1i;
        s_t2 += t2i;
        m1_n += t1i*x;
        m2_n += t2i*x;
        //         cout << s_t1 << " " << s_t2 << " " << m1_n << " " << m2_n << endl;
      }
      //       cout << m1_n/s_t1 << " " << m2_n/s_t2 << endl;
      double s1_n =0;
      double s2_n =0;
      for (auto it1=it0->second.begin();it1!=it0->second.end();++it1){
        const double t1=distr[it0->first].t1;
        const double t2=1-t1;
        const double m1=distr[it0->first].m1;
        const double m2=distr[it0->first].m2;
        const double s =starting_point(get<0>(*it1));
        const double x =get<1>(*it1)*s;
        const double s1=distr[it0->first].s1;
        const double s2=distr[it0->first].s2;
        const double e1=(s1*s1*128<(x-m1)*(x-m1))?0.00311664:t1*exp(-(x-m1)*(x-m1)/(2*s1*s1))/(s1*sqrt(2*pi));
        const double e2=(s2*s2*128<(x-m2)*(x-m2))?0.00311664:t2*exp(-(x-m2)*(x-m2)/(2*s2*s2))/(s2*sqrt(2*pi));
        //         const double t1i =(abs(e1)+abs(e2)<1e-10)?(e1>e2?1:0):e1/(e1+e2);
        //         const double t2i =(abs(e1)+abs(e2)<1e-10)?(e2>e1?1:0):e2/(e1+e2);
        const double t1i =e1/(e1+e2);
        const double t2i =e2/(e1+e2);
        //         cout << t1 << " " << x << " " << m1 << " " << m2 << " " << s1 << " " << s2 << " " << e1 << " " << e2 << " " << t1i << " " << t2i << endl;
        s1_n += t1i*(x-m1_n/s_t1)*(x-m1_n/s_t1); // dot product
        s2_n += t2i*(x-m2_n/s_t2)*(x-m2_n/s_t2); // dot product
        //         cout << s1_n << " " << s2_n << endl;
      }
      //       cout << sqrt(abs(s1_n/s_t1)) << " " << sqrt(abs(s2_n/s_t2)) << endl;
      distr[it0->first].t1 = s_t1/n;
      distr[it0->first].m1 = m1_n/s_t1;
      distr[it0->first].m2 = m2_n/s_t2;
      distr[it0->first].s1 = sqrt(abs(s1_n/s_t1));
      distr[it0->first].s2 = sqrt(abs(s2_n/s_t2));
      if (distr[it0->first].m2<distr[it0->first].m1){
        swap(distr[it0->first].m2,distr[it0->first].m1);
        swap(distr[it0->first].s2,distr[it0->first].s1);
        distr[it0->first].t1=1-distr[it0->first].t1;
      }
      if (distr[it0->first].s1<max(distr[it0->first].m2,distr[it0->first].m1)/10){
        distr[it0->first].s1=max(distr[it0->first].m2,distr[it0->first].m1)/10;
      }
      if (distr[it0->first].s1<max(abs(distr[it0->first].m2),abs(distr[it0->first].m1))/10){
        distr[it0->first].s1=max(abs(distr[it0->first].m2),abs(distr[it0->first].m1))/10;
      }
      if (distr[it0->first].t1<0.015625) distr[it0->first].t1=0.015625;
      if (distr[it0->first].t1>0.984375) distr[it0->first].t1=0.984375;
      if (isnan(distr[it0->first].m1)||isnan(distr[it0->first].m2)||isnan(distr[it0->first].s1)||isnan(distr[it0->first].s2)){
        cout <<  "maan schon wieder ein nan :( " << endl;
        exit(-1);
      }
    }
    //     cout << mean << " " << distr[it0->first].m1 << " " << distr[it0->first].m2 << endl;
  }
}

int main(int argc, char *argv[]){
  
  return 0;  
}
