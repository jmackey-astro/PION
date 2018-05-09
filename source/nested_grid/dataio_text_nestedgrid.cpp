/// \file dataio_text_nestedgrid.cpp
/// \author Jonathan Mackey
///
/// Class definitions for ASCII Text data I/O for a nested grid
/// 

// ##################################################################
// ##################################################################


int dataio_text::OutputData(
      const string outfile,
      vector<class GridBaseClass *> &cg,  ///< address of vector of grid pointers.
      class SimParams &SimPM,  ///< pointer to simulation parameters
      const long int counter   ///< number to stamp file with (e.g. timestep)
      )
{
  if (!cg[0])
    rep.error("dataio_text::output_ascii_data() null pointer to grid!",cg[0]);

  for (int l=0; l<SimPM.grid_nlevels; l++) {

    // for now write a different file for each level in the nested grid.
    ostringstream temp; temp << outfile << "_level";
    temp.width(2); temp.fill('0');
    temp << l;
    dataio_text::gp = cg[l];

    cout <<"dataio_text::OutputData() writing data.\n";
    string fname = set_filename(temp.str(), counter);
    int err = output_ascii_data(fname, SimPM);
    cout <<"dataio_text::OutputData() written data.\n";

  }
  return err;
}



