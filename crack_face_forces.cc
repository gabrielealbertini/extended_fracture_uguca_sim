#include "static_communicator_mpi.hh"
#include "uca_parameter_reader.hh"
#include "material.hh"
#include "uca_simple_mesh.hh"

#include "unimat_normal_interface.hh"
//#include "bimat_interface.hh"
#include "fracture_law.hh"
//#include "barras_law.hh"

#include <sys/time.h>
#include <stdio.h>
#include <iostream>
#include <cmath>

using namespace uguca;

int main(int argc, char *argv[]) {

  // communicator for parallel simulation
  StaticCommunicatorMPI * comm = StaticCommunicatorMPI::getInstance();
  int world_rank = comm->whoAmI();
  int world_size = comm->getNbProc();

  struct timeval t0, t1;
  gettimeofday(&t0, NULL);

  // get input file name
  std::string fname;
  if(argc<2) {
    if (world_rank==0) {
      std::cerr << "Not enough arguments:"
		<< " ./crack_face_forces <input_file>" << std::endl;
    }
    return 1;
  }
  else {
    fname = argv[1];
  }
  
  if (world_rank==0)
    std::cout << "simulation_code = weak-interface" << std::endl;
  

  // read input file
  ParameterReader data;
  data.readInputFile(fname);
  
  // mesh
  double x_length   = data.get<double>("x_length");
  int nb_x_elements = data.get<int>("nb_x_elements");
  
  int dim = 2;
  double z_length = 0.0;
  double nb_z_elements = 1;
  if (data.has("z_length")) {
    dim = 3;
    z_length   = data.get<double>("z_length");
    nb_z_elements = data.get<int>("nb_z_elements");  
  }
  
  SimpleMesh mesh = dim==2 ?
    SimpleMesh(x_length, nb_x_elements) :
    SimpleMesh(x_length, nb_x_elements,
	       z_length, nb_z_elements);
  
  
  // constitutive interface law
  FractureLaw law(mesh,
		  data.get<double>("tauc"),
		  data.get<double>("dc"));

  // materials
  Material top_mat = dim==3 ?
    Material(data.get<double>("E_top"),
	     data.get<double>("nu_top"),
	     data.get<double>("rho_top")) :
    Material(data.get<double>("E_top"),
	     data.get<double>("nu_top"),
	     data.get<double>("rho_top"),
	     data.get<bool>("pstress_top"));
  top_mat.readPrecomputedKernels();
   
  // interface
  UnimatNormalInterface interface(mesh, top_mat, law);

  // external loading
  interface.getLoad().component(0).setAllValuesTo(data.get<double>("shear_load"));
  interface.getLoad().component(1).setAllValuesTo(data.get<double>("normal_load"));

  // time step
  double duration = data.get<double>("duration");
  double time_step = data.get<double>("tsf") * interface.getStableTimeStep();
  int nb_time_steps = duration/time_step;
  interface.setTimeStep(time_step);

  // initialization
  interface.init();

  // notch
  double * X = mesh.getLocalCoords()[0];
  NodalFieldComponent & tau_max = law.getTauMax();
  double notch_length = data.get<double>("notch_length");
  
  for (int i=0;i<mesh.getNbLocalNodes(); ++i) 
    if (std::abs(X[i] - x_length/2.) < (notch_length))
      tau_max(i) = 0.0;
  
  // pressure on crack face
  double pressure = data.get<double>("pressure");
  double pressure_length = data.get<double>("pressure_length");
  
  NodalFieldComponent & norm_load = interface.getLoad().component(1);
  for (int i=0;i<mesh.getNbLocalNodes(); ++i)
  if (std::abs(X[i] - x_length/2.) < pressure_length)
    norm_load(i) = pressure;
    
  // interface heterogeneity
  if (data.has("het_tauc")) {
    NodalFieldComponent & dc = law.getDc();
    double xstart = data.get<double>("het_xstart");
    double xend   = data.get<double>("het_xend");
    double het_tauc   = data.get<double>("het_tauc");
    double het_dc   = data.get<double>("het_dc");
      
    for (int i=0;i<mesh.getNbLocalNodes(); ++i) {
      if ((std::abs(X[i] - x_length/2.) > xstart) &&
	  (std::abs(X[i] - x_length/2.) < xend)) {
	tau_max(i) = het_tauc;
	dc(i) = het_dc;
      }
    }
  }

  // prevent crack growth on left
  bool unilateral_crack_growth = false;
  if (data.has("unilateral_crack_growth")) 
    unilateral_crack_growth = data.get<bool>("unilateral_crack_growth");


  // dumping
  std::string simulation_name = data.get<std::string>("simulation_name");
  std::string dump_folder = data.get<std::string>("dump_folder");
  if (world_rank==0)
    std::cout << "output_folder = " << dump_folder << std::endl;

  // copy input file
  if (world_rank==0) {
    std::ifstream src_ifile(fname, std::ios::binary);
    std::ofstream dst_ifile(dump_folder + simulation_name + ".in", std::ios::binary);
    dst_ifile << src_ifile.rdbuf();
    src_ifile.close();
    dst_ifile.close();
  }

  std::string dumper_bname = simulation_name;
  if (world_rank==0) 
    std::cout << "dumper_bname = " << dumper_bname << std::endl;
  
  interface.initDump(dumper_bname,
		     dump_folder,
		     Dumper::Format::Binary);
    
  interface.registerDumpFields(data.get<std::string>("dump_fields"));
  unsigned int dump_int = std::max(1, int(data.get<double>("dump_interval")/time_step));

  // setup of restart to dump
  Restart restart(data.get<std::string>("simulation_name"),
		  data.get<std::string>("restart_dir"),
		  Restart::Format::ASCII);
  
  restart.registerIO(interface.getTop().getDisp());
  //restart.registerIO(interface.getBot().getDisp());
  
  double restart_dump_time = 0.0;
  if (data.has("restart_dump_time")) {
    restart_dump_time = data.get<double>("restart_dump_time");
    if (restart_dump_time > nb_time_steps){
      std::cout<<"restart_dump_time > nb_time_steps !!!"<<std::endl;
      return 1;
    }
  }
  // quasi dynamic integration
  double quasi_dynamic_start_time = 0.0;
  double quasi_dynamic_end_time = 0.0;
  if (data.has("quasi_dynamic_start_time")){
    quasi_dynamic_start_time = data.get<double>("quasi_dynamic_start_time");
    quasi_dynamic_end_time = data.get<double>("quasi_dynamic_end_time");
  }

  // time stepping
  int s=1;

  // check if restart from file is required
  if (data.has("restart_load_time")) {
    s = data.get<double>("restart_load_time")/time_step;
     
    // different from dumper not to overwrite it
    Restart restart_load = restart;
    restart_load.initIO(data.get<std::string>("restart_name"),
			data.get<std::string>("restart_dir"),
			Restart::Format::ASCII);
    restart_load.load(s);
    if (world_rank==0)
      std::cout<<"load restart -- s=" << s << "/" << nb_time_steps << std::endl<<std::endl;
    bool dynamic = false;
    interface.advanceTimeStep(dynamic);
  }
  else { // only dump when not restarted
    interface.dump(0,0);
  }

  //nucleation
  double *tau_max_0nuc = new double[mesh.getNbLocalNodes()];
 
  double nuc_r        = data.get<double>("nuc_r");
  double nuc_tstart   = data.get<double>("nuc_tstart");
  double nuc_tend     = data.get<double>("nuc_tend");
  double nuc_x_center = data.get<double>("nuc_x_center");
  double nuc_z_center = dim==3 ? data.get<double>("nuc_z_center") : 0.0;
  double nuc_tauc = data.get<double>("nuc_tauc");
  bool   nuc_sym  = data.get<bool>("nuc_sym");
  double nuc_xc0  = data.get<double>("nuc_xc0");
  
  int hf_dump_tstart = nuc_tstart/time_step;
  
  // main time loop
  for (; s<=nb_time_steps; ++s) {
    
    if (world_rank==0) {
      double seconds = t1.tv_sec + t1.tv_usec*1e-6;
      gettimeofday(&t1, NULL);
      seconds = t1.tv_sec + t1.tv_usec*1e-6-seconds;
      long minutes = seconds/60;
      long hours = minutes/60;
      seconds -= minutes*60;
      minutes -=hours*60;
      printf("s=%d/%d wall time %dh %02dmin %2.3fsec\n",
	     s,nb_time_steps,int(hours),int(minutes),seconds); 
      std::cout.flush();
    }
    
    // nucleation
    if (data.has("nuc_r")) {
      if (s==int(nuc_tstart/time_step)){
	// copy tau_max
	for (int i=0;i<mesh.getNbLocalNodes(); ++i)
	  tau_max_0nuc[i] = tau_max(i);
	// unilateral prop
	if (unilateral_crack_growth){
	  for (int i=0;i<mesh.getNbLocalNodes(); ++i) {
	    if (X[i] < x_length/2. - notch_length) {
	      tau_max(i) = 100.0*tau_max(i);
	    }
	  }
	}
      }
      double * Z = dim==3 ? mesh.getLocalCoords()[2] : NULL;
      
      for (int i=0;i<mesh.getNbLocalNodes(); ++i) {
	if (!nuc_sym)
	  if (X[i]<x_length/2.)
	    continue;

	double px = std::abs(X[i] - x_length/2.) - nuc_x_center;
	double pz = dim==3 ? std::abs(Z[i] - z_length/2.) - nuc_z_center : 0.0;
	double r_outher = std::min(std::max(0.0,
					    (s*time_step - nuc_tstart)/(nuc_tend-nuc_tstart)),
				   1.0)*nuc_r;
	double r_inter = std::max(0.0,r_outher - nuc_xc0);
	double radius = std::sqrt(px*px+pz*pz);

	if (radius < r_outher) {
	  double nuc_tauc_here = nuc_tauc;
	  if (radius > r_inter) {
	    double d = (r_outher-radius)/nuc_xc0;
	    nuc_tauc_here = std::max(nuc_tauc, tau_max_0nuc[i] - d*(tau_max_0nuc[i]-nuc_tauc));
	  }
	  
	  tau_max(i) = std::min(tau_max(i),nuc_tauc_here);
	}
      }
    }
    
    // time integration
    bool dynamic = true;
    //  

    if ((s <= int(quasi_dynamic_end_time/time_step)) &&
	(s > int(quasi_dynamic_start_time/time_step))) {
      dynamic = false;
    }

  interface.advanceTimeStep(dynamic);

    // dump
    if (s>hf_dump_tstart) {
      if (s % dump_int == 0)
	interface.dump(s,s*time_step);
    }
    else {
      if (s % (10*dump_int) == 0)
	interface.dump(s,s*time_step);
    }
    
    // dump restart
    if (s  == int(restart_dump_time/time_step)){
      if (world_rank==0) {  
	std::cout<<"dump restart -- s=" << s << "/" << nb_time_steps << std::endl;
	std::cout.flush();
      }
      restart.dump(s);
    }
  }

  if (world_rank==0) {
    gettimeofday(&t1, NULL);
    double seconds = t1.tv_sec - t0.tv_sec +(t1.tv_usec-t0.tv_usec) *1e-6;
    long minutes = seconds/60;
    long hours = minutes/60;
    seconds -= minutes*60;
    minutes -=hours*60;
    std::cout<<"went to the end"<<std::endl;
    std::cout<<"mpi "<< world_size<<std::endl;
    printf("wall time %dh %02dmin %2.3fsec\n",int(hours),int(minutes),seconds);
  }
  delete[] tau_max_0nuc;

  comm->finalize();
  return 0;
}
