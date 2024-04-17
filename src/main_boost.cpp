#include <fsm/fill_spill_merge.hpp>

#include <richdem/common/Array2D.hpp>
#include <richdem/common/gdal.hpp>
#include <richdem/misc/misc_methods.hpp>
#include <richdem/ui/cli_options.hpp>

#include <string>
#include <fstream>

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

namespace rd = richdem;
namespace dh = richdem::dephier;


template<class elev_t>
void save_dh(richdem::dephier::DepressionHierarchy<elev_t> &deps, const char * filename){
    // make an archive
    std::ofstream ofs(filename);
    boost::archive::binary_oarchive oa(ofs);
    oa << deps;
}

//template<class elev_t>
void restore_dh(std::vector<richdem::dephier::Depression<double>> &deps, const char * filename){
    // open the archive
    std::ifstream ifs(filename);
    boost::archive::binary_iarchive ia(ifs);
    ia >> deps;
}


int main(int argc, char **argv){
  CLI::App app("Fill-Spill-Merge Example Program");

  std::string topography_filename;
  std::string output_prefix;
  double      surface_water_level = std::numeric_limits<double>::quiet_NaN();
  std::string surface_water_filename;
  double      ocean_level;
  std::string save_dh_filename;
  std::string load_dh_filename;

  size_t num;
  size_t len;
  app.add_option("topography", topography_filename, "Topography to run FSM on")->required();
  app.add_option("output",     output_prefix,       "Path of GeoTiff output file")->required();
  const auto swl_ptr = app.add_option("--swl", surface_water_level, "Surface water level as a numeric constant");
  app.add_option("--swf", surface_water_filename, "File containing surface water levels")->excludes(swl_ptr);
  app.add_option("ocean_level",         ocean_level, "Elevation of the ocean")->required();
  const auto save_dh_ptr = app.add_option("--save_dh",    save_dh_filename, "Filename where you would like the depression hierarchy to be saved for reuse (optional, requires Boost)");
  app.add_option("--load_dh",    load_dh_filename, "Filename where you would like the depression hierarchy to be loaded from, if it has been computed before (optional, requires Boost)")->excludes(save_dh_ptr);


  CLI11_PARSE(app, argc, argv);

  std::cout<<"m Input DEM           = "<<topography_filename   <<std::endl;
  std::cout<<"m Output prefix       = "<<output_prefix         <<std::endl;
  std::cout<<"m Surface water level = "<<surface_water_level   <<std::endl;
  std::cout<<"m Surface water file  = "<<surface_water_filename<<std::endl;
  std::cout<<"m Ocean level         = "<<ocean_level           <<std::endl;
  std::cout<<"m Filename to save DH for reuse         = "<<save_dh_filename <<std::endl;
  std::cout<<"m Filename to open previously saved DH  = "<<load_dh_filename <<std::endl;


  rd::Timer timer_io;
  timer_io.start();
  rd::Array2D<double> topo(topography_filename);
  timer_io.stop();

  std::cout<<"m Data width  = "<<topo.width ()<<std::endl;
  std::cout<<"m Data height = "<<topo.height()<<std::endl;
  std::cout<<"m Data cells  = "<<topo.numDataCells()<<std::endl;

  rd::Array2D<double> wtd(topography_filename);
  if(surface_water_filename.empty()){
    wtd = rd::Array2D<double>(topo,surface_water_level);
    //All cells have the same amount of water
  } else {
    wtd = rd::Array2D<double>(surface_water_filename);
  }

  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );      //No cells are part of a depression
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);      //No cells flow anywhere

  wtd.setNoData(topo.noData());

  //Label the ocean cells. This is a precondition for using
  //`GetDepressionHierarchy()`.
  rd::BucketFillFromEdges<rd::Topology::D8>(topo, label, ocean_level, dh::OCEAN);

  //Make NoData cells also ocean cells. Ocean has no water on it to begin with.
  #pragma omp parallel for
  for(unsigned int i=0;i<label.size();i++){
    if(topo.isNoData(i) || label(i)==dh::OCEAN){
      label(i) = dh::OCEAN;
      wtd  (i) = 0;
    }
  }

  rd::Timer timer_calc;
  timer_calc.start();

  //Generate flow directions, label all the depressions, and get the hierarchy
  //connecting them

  std::vector<richdem::dephier::Depression<double>> deps;

  //We are not loading in a DH, so we need to calculate one
  if(load_dh_filename.empty()){
    deps = dh::GetDepressionHierarchy<double,rd::Topology::D8>(topo, label, flowdirs);
  }
  //We are loading in a DH, so load in the data here and do not calculate it.
  //We need to load in 3 pieces of data: the DH, the labels array, and the flowdirs array.
  else{


    try{
      restore_dh(deps, (load_dh_filename+"_DH.txt").c_str());           //load the DH

      label = rd::Array2D<dh::dh_label_t>(load_dh_filename+"_label.tif");     //load the labels array

      std::ifstream infile(load_dh_filename+"_flowdirs.txt");                 //load the flowdirs

      // Check if the file is opened successfully
      if (infile.is_open()) {
        // Read the data from the file into the array
        int num;
        int i = 0;
        while (infile >> num) {
            flowdirs(i++) = num;
        }
        // Close the file
        infile.close();
      }
    }
    catch(std::exception& e){
      throw std::runtime_error("No such file! Check your filename in --load_dh!");
    }
  }

  //run FSM on the DH:
  dh::FillSpillMerge(topo, label, flowdirs, deps, wtd);

  timer_calc.stop();

  timer_io.start();

  //If we are saving the data from this depression hierarchy to reuse with FSM later, do so here:
  if(!save_dh_filename.empty()){

    //save the actual DH:
    save_dh(deps, (save_dh_filename+"_DH.txt").c_str());

    //save the labels array:
    label.projection = topo.projection;
    label.saveGDAL(save_dh_filename+"_label.tif");


    //save the flow directions array:
    std::ofstream outfile_flowdirs(save_dh_filename+"_flowdirs.txt");

    // Check if the file is opened successfully
    if (outfile_flowdirs.is_open()) {
        // Write the contents of the array to the file
        for (int i = 0; i < sizeof(flowdirs) / sizeof(flowdirs(0)); ++i) {
            outfile_flowdirs << flowdirs(i) << " ";
        }

        // Close the file
        outfile_flowdirs.close();
        std::cout << "Array saved to file successfully." << std::endl;
    } else {
        std::cerr << "Unable to open file." << std::endl;
    }
  }

  //Output the water table depth
  wtd.saveGDAL(output_prefix+"-wtd.tif");

  for(unsigned int i=0;i<topo.size();i++)
    if(!topo.isNoData(i))
      wtd(i) += topo(i);

  //Output the new height of the hydraulic surface
  wtd.saveGDAL(output_prefix+"-hydrologic-surface-height.tif");

  timer_io.stop();

  std::cout<<"Finished."<<std::endl;
  std::cout<<"IO time   = "<<timer_io.accumulated()  <<" s"<<std::endl;
  std::cout<<"Calc time = "<<timer_calc.accumulated()<<" s"<<std::endl;


  return 0;
}
