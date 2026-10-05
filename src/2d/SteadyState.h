/**
 * @file SteadyState.h
 * @brief 稳态问题类
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-12
 */

# ifndef __STEADYSTATE__H__
# define __STEADYSTATE__H__

# include <memory>
# include <vector>
# include <sstream>
# include <iostream>
# include <sys/time.h>
# include <omp.h>
# include <assert.h>
# include "Controller.h"
# include "ConfigFile.h"
# include "config.h"

template <typename SOLVER >
class SteadyState{
	public:
	    typedef SOLVER solver_t; 
	    typedef typename solver_t::Residual Residual;
	    typedef typename Residual::Problem Problem;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	    enum {dim = dis_t::dim};
	public:
	    bool is_stop_run;
	    bool is_exit;
	    bool is_read;
	private:
	    int ORDER;
	    value_t tol;
	    int n_thread;		  
	    std::unique_ptr<solver_t> p_solver;
	    ptr_mesh_t p_mesh;
	    std::unique_ptr<sol_t> sol;
	    std::unique_ptr<sol_t> rhs;
	    double res_norm;
	private:
	    bool is_converged() const {return res_norm < tol;}
	    void allocMemory();
	    void freeMemory();
	public:
		SteadyState() : is_stop_run(true),
						is_exit(false), is_read(false),
		  			p_solver(new solver_t()),
		  			p_mesh(new mesh_t()) {};

		void initialize(){
		  	if(sol.get() == NULL) alloc_memory();
		  	p_solver->initialize(*sol, *p_mesh);
		  	rhs->Reinit(*sol);
		}

		void setMaxwellian(){
		  	if(sol.get()==NULL) alloc_memory();
		  	std::cerr << "Set the solution as Maxwellian ... " << std::endl;
		  	Problem::initialize(*sol, *p_mesh);
		  	rhs->Reinit(*sol);
		  	std::cerr << "Set initial value OK! " << std::endl;
		}

		void stop_run() { is_stop_run = true; }

		void run(){
		  	if( !is_read ) {
		  		 std::cerr << "Please read the config file (Default: input.txt) \
		  			  or load data (Default: DATA.dat) from previous rtest first!" 
		  			  << std::endl;
		  	}

		  	std::cerr << "Entering run ... " << std::endl;
		  	is_stop_run = false;

		  	timeval start_time;
		  	gettimeofday( &start_time, NULL );

		  	static unsigned int n_run = 0;    /**< 存储 run 函数调用次数 */
		  	++n_run;
		  	if( n_run == 1 ){
		  		 std::ofstream clear_file( "res_step", std::ios::out );
		  		 clear_file.close();
		  	}
		  	std::ofstream os_res( "res_step", std::ios::app );
		  	os_res.precision( 12 );
		  	os_res << "% the n-th run: n = " << n_run << std::endl;

		  	if(sol.get() == NULL) initialize();

		  	unsigned int steps = 0;
		  	//OPENMP parallel
		  	omp_set_num_threads(n_thread);

		  	Residual::getResidualPL2Norm(res_norm, *sol, *rhs, *p_mesh);
		  	os_res << steps << " " << res_norm << " " << 0 << std::endl;

		  	std::cout.precision(12);
		  	std::cout << "steps = " << steps 
		  			    << ":  residual = " << res_norm
		  		       << std::endl;

		  	while((! (is_stop_run || is_converged() || is_exit))){
		  		 steps++;
		  		 p_solver->solve(*sol, *rhs, *p_mesh);
		  		 p_solver->postProcessing(*sol, *p_mesh);
		  		 Residual::getResidualPL2Norm(res_norm, *sol, *rhs, *p_mesh);
		  		 timeval end_time;
		  		 gettimeofday( &end_time, NULL );
		  		 std::cout << "steps = " << steps 
					 	   << ": residual = " << res_norm
		  			       << std::endl;
		  		 os_res << steps << " "
		  			    << res_norm << " "
		  				<< getElapsedTime(start_time, end_time) << std::endl ;
		  		 getControl();
		  	}

		  	if( is_exit ){
		  	std::cerr << "The Pro is stopped!" << std::endl;
		  	} else if (is_stop_run){
		  	std::cerr << "Run is stopped." << std::endl;
		  	}

		  	timeval end_time;
		  	gettimeofday( &end_time, NULL );

		  	std::cerr << "The elapsed steps is " << steps << std::endl;
		  	std::cerr << "The elapsed time is " << getElapsedTime(start_time,
		  			  end_time) << std::endl ;

		  	os_res << "The elapsed steps is " << steps << std::endl;
		  	os_res << "The elapsed time is " << getElapsedTime(start_time,
		  			  end_time) << std::endl;

		  	if(is_converged()){
		  		 std::stringstream file_data;
		  		 file_data << "DATA" << n_run << "th.dat";
		  		 dump_data( file_data.str() ) ;
		  		 Problem::outputSolution(*sol, *p_mesh);
		  	}

		  	os_res << std::endl;
		  	os_res.close();
		}

		void read_config(const std::string& filename){
			std::cerr << "Reading config file " << filename << " ... " << std::flush;
			ConfigFile cf( filename );
			try{
				 ORDER = (int)cf.Value("System","ORDER");
				 tol = (value_t)cf.Value("Error", "Tol" );
				 n_thread = (int)cf.Value("System", "n_thread");

	 		    p_mesh->read_config( cf );

				 Problem::read_config( p_mesh, cf );

				 p_solver->read_config( cf );

				 is_read = true;
				 std::cerr << "OK!" << std::endl;
			 }
			 catch (std::string s){
				  std::cerr << s << std::endl;
				  exit(1);
			 }
		}
		  
	   void alloc_memory(){
			size_t Nx = p_mesh->n_ele_x();
			size_t Ny = p_mesh->n_ele_y();
			sol.reset(new sol_t(Nx, Ny, ORDER));
			rhs.reset(new sol_t(Nx, Ny, ORDER));
		}

		double getElapsedTime(const timeval& _start, const timeval& _end){
			 return _end.tv_sec - _start.tv_sec + 
					 (_end.tv_usec - _start.tv_usec) / 1000000.;
		}

		void reset_tol(const double _tol){
		    std::cout << " Reset tol as :" << _tol << std::endl ;
		    tol = _tol;
		}

		void reset_order(int order){
			std::cerr << "Reset order as " << order << std::endl;
			ORDER = order;
		    size_t Nx = p_mesh->n_ele_x();
			size_t Ny = p_mesh->n_ele_y();
			for(size_t ny = 0; ny < Ny; ny++){
			    for(size_t nx = 0; nx < Nx; nx++){
					(*sol)[nx][ny].ResetOrder(ORDER);
			   	   	(*rhs)[nx][ny].ResetOrder(ORDER);
			    }
			}
		}

		void print_config(){
			Problem::print_config();
			std::cerr << " Order = " << ORDER << std::endl;
			std::cerr << " Nx = " << p_mesh->n_ele_x() << std::endl;
			std::cerr << " Ny = " << p_mesh->n_ele_y() << std::endl;
			std::cerr << " tol = " << tol << std::endl;
		    p_solver->print_config();
			std::cerr << std::flush;
		}

		void dump_data(const std::string& filename){
			std::cerr << "Dumping data to file " << filename << " ... " << std::
			     flush;
			std::ofstream os(filename.c_str(), std::ios::binary);
			p_mesh->dump_data(os);
			dump_config(os);
			dump_solution(os);
			Problem::dump_data(os);
		    p_solver->dump_data(os);
		}

		void dump_config(std::ostream& os) const {
			os.write((const char *)&ORDER, sizeof(int));
			os.write((const char *)&tol, sizeof(double));
			os.write((const char *)&n_thread, sizeof(int));
		}

		void dump_solution(std::ostream& os) const {
			size_t Nx = p_mesh->n_ele_x();
			size_t Ny = p_mesh->n_ele_y();
			size_t n_ele = Nx * Ny;
			std::cout << " dump solution " << std::endl;
			os.write((const char *)&(n_ele), sizeof(size_t));
			for(size_t ny = 0; ny < Ny; ny++){
			    for(size_t nx = 0; nx < Nx; nx++){
					const dis_t& dis = (*sol)[nx][ny];
					for(size_t j = 0; j < dim; j++){
					  os.write((const char *)&(dis.Center()[j]),
								sizeof(double));
					}
					os.write((const char *)&(dis.Scaling()),
						  sizeof(double));
					typename dis_t::ConstIterator the_entry = dis.begin(),
						end_entry = dis.end();
			   	 	for (; the_entry != end_entry; ++ the_entry){
			   	     os.write((const char *)&(*the_entry), sizeof(double));
					}
				}
			}
		}

		void load_data(const std::string& filename){
			std::cerr << "Loading data from file " << filename
			   		  << " ... " << std::flush;
			//free_memory();

			std::ifstream is(filename.c_str(), std::ios::binary);

				p_mesh.reset(new mesh_t(is));

			    load_config(is);
				alloc_memory();
				load_solution(is, *sol);
			(*rhs).Reinit(*sol);
			Problem::load_data(is);
			p_solver.reset(new solver_t(p_mesh, is));
			is_read = true;
			std::cerr << "OK!" << std::endl;
			print_config();
		}

		void load_config(std::istream& is){
			is.read((char *)&ORDER, sizeof(int));
			is.read((char *)&tol, sizeof(double));
			is.read((char *)&n_thread, sizeof(int));
		}

		void load_solution(std::istream& is, sol_t& _sol) const {
			size_t n_ele;
			is.read((char *)&n_ele, sizeof(size_t));
			size_t Nx = p_mesh->n_ele_x();
			size_t Ny = p_mesh->n_ele_y();
			std::cout << " Nx = " << Nx
			   		  << " Ny = " << Ny
			   		  << " n = " << n_ele 
			   		  << " n_ele = " << _sol.n_ele() << std::endl;
			assert((n_ele == _sol.n_ele()));
			for(size_t ny = 0; ny < Ny; ny++){
			     for(size_t nx = 0; nx < Nx; nx++){
			   		dis_t& dis = _sol[nx][ny];
			   		for(size_t j = 0; j < dim; j++)
			   			 is.read((char *)&(dis.Center()[j]), sizeof(double));
			   		is.read((char *)&(dis.Scaling()), sizeof(double));
			   		dis.ResetOrder(ORDER);
			   		typename dis_t::Iterator the_entry = dis.begin(),
			   					end_entry = dis.end();
			   		for(; the_entry != end_entry; the_entry++){
			   			 is.read((char *)&(*the_entry), sizeof(double));
			   		}
			     }
			}
		}
};

# endif
