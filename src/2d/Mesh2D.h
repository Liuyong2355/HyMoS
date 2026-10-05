/**
 * @file mesh2d.h
 * @brief 二维网格类
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-07
 */

#ifndef __MESH2D__H__
#define __MESH2D__H__

# include <vector>
# include <string>
# include <fstream>
# include <memory>

# include "ConfigFile.h"

template <typename VALUE>
class Mesh2D{
	public:
	    typedef VALUE point_t;
	    typedef std::vector<point_t> Base;
	private:
	    Base x_vec, y_vec;
	    size_t Nx, Ny;
	    point_t x_l, x_r, y_b, y_t;
	public:
		Mesh2D(){};
		Mesh2D(const size_t Nx, const size_t Ny){Reinit(Nx, Ny);}
		//Mesh2D(const std::string str){ReadMesh(str);}
		Mesh2D(const point_t& x0, const point_t& x1, 
			   const point_t& y0, const point_t& y1, 
			   const size_t Nx, const size_t Ny) : 
			x_l(x0), x_r(x1), y_b(y0), y_t(y1){Reinit(Nx, Ny);}

	    Mesh2D(std::istream& is){
			  is.read((char *)&Nx, sizeof(size_t));
			  x_vec.resize(Nx+1);
			  for(size_t i = 0; i < Nx+1; i++)
					is.read((char *)&(x_vec[i]), sizeof(double));
			  is.read((char *)&Ny, sizeof(size_t));
			  y_vec.resize(Ny+1);
			  for(size_t i = 0; i < Ny+1; i++)
					is.read((char *)&(y_vec[i]), sizeof(double));
			  x_l = x_vec[0];
			  x_r = x_vec[Nx];
			  y_b = y_vec[0];
			  y_t = y_vec[Ny];
		 }
		 
	    void Reinit(const size_t Nx, const size_t Ny){
			  x_vec.resize(Nx + 1);
			  point_t dx = (x_r - x_l) / Nx; 
			  for(size_t i = 0; i < Nx; i++){
					x_vec[i] = x_l + i * dx;
			  }
			  x_vec[Nx] = x_r;

		     y_vec.resize(Ny + 1);
		     point_t dy = (y_t - y_b) / Ny;
		     for(size_t i = 0; i < Ny; i++){
					y_vec[i] = y_b + i * dy;
			  } 
			  y_vec[Ny] = y_t;
		 }

		 point_t dx(const size_t i, const size_t j) const{
			  return (x_vec[i+1] - x_vec[i]);
		 }

		 point_t dy(const size_t i, const size_t j) const{
			  return (y_vec[j+1] - y_vec[j]);
		 }

	    point_t Area(const size_t i, const size_t j) const{ 
			  return (x_vec[i+1] - x_vec[i]) * (y_vec[j+1] - y_vec[j]);
		 }

		 void center(std::vector<double>& center, size_t nx, size_t ny) const {
			  center.resize(2);
			  center[0] = 0.5 * (x_vec[nx] + x_vec[nx+1]);
			  center[1] = 0.5 * (y_vec[ny] + y_vec[ny+1]);
		 }

		 size_t n_ele_x() const{
			  return Nx;
		 }

	    size_t n_ele_y() const{
			  return Ny;
		 }

		 size_t n_ele () const{
			  return Nx * Ny;
		 }

		 double DomainSize() const {
			  return (x_r - x_l) * (y_t - y_b);
		 }

		 void Coarser(){
			  Nx = Nx / 2;
			  Ny = Ny / 2;
			  this->Reinit(Nx, Ny);
		 }

		 void Finer(){
			  Nx = Nx * 2;
			  Ny = Ny * 2;
			  this->Reinit(Nx, Ny);
		 }

		 void read_config(const ConfigFile& cf){
			  int m_t = (int)cf.Value("Mesh", "Type");
		     if( m_t == 0 ){
				   Nx = (int)cf.Value("Mesh","Nx");
				   x_l = (double)cf.Value("Mesh","Lx");
				   x_r = x_l + (double)cf.Value("Mesh","Lenx");

					Ny = (int)cf.Value("Mesh","Ny");
				   y_b = (double)cf.Value("Mesh","Ly");
				   y_t = y_b + (double)cf.Value("Mesh","Leny");

					Reinit(Nx, Ny);
			  }
		 }

		 void dump_data(std::ostream& os) const {
			  os.write((const char *)&(Nx), sizeof(size_t));
			  for(size_t i = 0; i < Nx+1; i++)
					os.write((const char *)&(x_vec)[i], sizeof(double));
			  os.write((const char *)&(Ny), sizeof(size_t));
			  for(size_t i = 0; i < Ny+1; i++)
					os.write((const char *)&(y_vec)[i], sizeof(double));
		 }
};
#endif
