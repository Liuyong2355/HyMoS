/**
 * @file solution.h
 * @brief 二维问题解类，继承自二维vector
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-12
 */

# ifndef __SOLUTION__H__
# define __SOLUTION__H__

# include <vector>

template <typename DISTRIBUTION, typename MESH>
class Solution : public std::vector<std::vector<DISTRIBUTION>>{
	 public:
		  typedef std::vector<std::vector<DISTRIBUTION>> BASE;
		  typedef DISTRIBUTION dis_t;
		  typedef MESH mesh_t;

	 public:
		  Solution(){};
		  
		  Solution(size_t Nx, size_t Ny, int ORDER) : BASE(Nx,
					 std::vector<DISTRIBUTION>(Ny, DISTRIBUTION(ORDER))) {}

		  void Reinit(const Solution& sol) {
			  //取出行数
			  this->resize(sol.size());
			  #pragma omp parallel for			
			  for(size_t ny = 0; ny < sol.size(); ny++){
				  (*this)[ny].resize(sol[ny].size());
				  //列数sol[ny].size()
				  for(size_t nx = 0; nx < sol[ny].size(); nx++){
					  (*this)[ny][nx].Reinit(sol[ny][nx]);
				  }
			  }
		  }

		  unsigned int GetOrder(){
				return (*this)[0][0].GetOrder();
		  }

		  size_t n_ele(){
				size_t n_ele = (*this).size() * (*this)[0].size();
				return n_ele;
		  }

		  Solution& operator += (const Solution& sol){
			  #pragma omp parallel for
			  for(size_t ny = 0; ny < this->size(); ny++){
				  for(size_t nx = 0; nx < (*this)[ny].size(); nx++){
					  (*this)[ny][nx] += sol[ny][nx];
				  }
			  }
			  return *this;
		  }

		  Solution& operator -= (const Solution& sol){
			  #pragma omp parallel for
			  for(size_t ny = 0; ny < this->size(); ny++){
			  	 for(size_t nx = 0; nx < (*this)[ny].size(); nx++){
			  		  (*this)[ny][nx] -= sol[ny][nx];
			  	 }
			  }
			  return *this;
		  }

		  template <typename VALUE>
		  Solution& operator *= (const VALUE& s){
			  #pragma omp parallel for
			  for(size_t ny = 0; ny < this->size(); ny++){
				  for(size_t nx = 0; nx < (*this)[ny].size(); nx++){
					  (*this)[ny][nx] *= s;
				  }
			  }
			  return *this;
		  }

		  template <typename VALUE>
		  Solution& operator /= (const VALUE& s){
				(*this) *= 1.0 / s;

				return *this;
		  }
};

# endif
