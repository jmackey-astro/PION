//
// this contains an old class i wrote for holding the radiation source parameters.
// it turned out to be too complicated for its own good, so i got rid of it, but
// it may return someday, so i am keeping it here.
//

// *******************************************************************
// *************** RADIATION SOURCES FOR RAY-TRACING *****************
// *******************************************************************

#define RT_SRC_MONO 1 ///< monochromatic radiation source.
#define RT_SRC_DIFF 2 ///< not really a source, a ray to calculate column density along.

///
/// If we are doing Ray Tracing, we want to store a list of sources
/// into a globally accessible place, so put them here
///
class RadiationSources {
 public:
  RadiationSources();
  ~RadiationSources();
  ///
  /// For I/O this is the number of sources we expect to add:
  ///
  int Nsources;
  ///
  /// Return number of sources 
  ///
  int nsources();
  ///
  /// Add a source, with position and strength
  ///
  void add_src(const double,  ///< Strength
	       const double *, ///< position vector (MAX_DIM elements please!)
				 const int,      ///< type of source
				 const bool       ///< true if at infinity.
	       );
  ///
  /// Get source strength
  /// Assigns negative value on error.
  ///
  void src_strength(const int, ///< src_id
		    double *   ///< strength (output).
		    );
  ///
  /// Get source position, physical units. (returns MAX_DIM array).
  /// Assigns huge negative value on error
  ///
  void src_position(const int, ///< src_id
		    double *   ///< position (MAX_DIM array) (output).
		    );
  ///
  /// Reset source position (requires MAX_DIM array).
  /// This is needed since the RT algorithm moves sources to cell vertices.
  ///
  int reset_position(const int, ///< src_id
		     const double * ///< position (MAX_DIM array).
		     );
	///
	/// Get source type
	///
	int src_type(const int, ///< src_id
							 int *      ///< src type (output).
							 );
	///
	/// returns true if source is at infinity.
	///
	bool src_at_infinity(const int ///< src id.
	                     );
	///
	/// returns true if there is an Nth source, where N is input.
	///
	bool src_exists(const int ///< src id
									);
 private:
  int nsrc;
  std::vector<struct rad_src_info *> slist;
};
extern class RadiationSources RSP;

// *******************************************************************
// *************** RADIATION SOURCES FOR RAY-TRACING *****************
// *******************************************************************
class RadiationSources RSP;



//------------------------------------------------
//------------ Radiation Sources -----------------
//------------------------------------------------
RadiationSources::RadiationSources()
{
  nsrc=0;
  Nsources=0;
  slist.clear();
}

RadiationSources::~RadiationSources()
{
  for (unsigned int i=0;i<slist.size();i++) {
    if (slist[i]!=0)
      slist[i] = mem.myfree(slist[i]);
  }
  slist.clear();
}

int RadiationSources::nsources()
{
  return nsrc;
}

void RadiationSources::add_src(const double  s, ///< Strength
			       const double *p,  ///< position
						 const int type,      ///< type of source
						 const bool at_inf      ///< true if at infinity.
			       )
{
  struct rad_src_info *rs=0;
  rs = mem.myalloc(rs,1);
  rs->strength = s; 
  for (int v=0;v<MAX_DIM;v++) rs->position[v] = p[v];
  rs->id = slist.size();

	rs->type = type;
	rs->at_infinity = at_inf;

  slist.push_back(rs);
  nsrc++;
  return;
}

void RadiationSources::src_strength(const int id, ///< src_id
				    double *s     ///< strength
				    )
{
  if      (static_cast<unsigned int>(id) >= slist.size())
    *s = -VERY_LARGE_VALUE;
  else if (!(slist[id]))
    *s = -VERY_LARGE_VALUE;
  else
    *s = slist[id]->strength;

  return;
}

void RadiationSources::src_position(const int id, ///< src_id
				    double *p     ///< strength
				    )
{
  if      (static_cast<unsigned int>(id) >= slist.size())
    for (int v=0;v<MAX_DIM;v++)
      p[v] = -VERY_LARGE_VALUE;
  else if (!(slist[id]))
    for (int v=0;v<MAX_DIM;v++)
      p[v] = -VERY_LARGE_VALUE;
  else
    for (int v=0;v<MAX_DIM;v++)
      p[v] = slist[id]->position[v];

  return;
}

int RadiationSources::reset_position(const int id, ///< src_id
				     const double *p ///< position (MAX_DIM array).
				     )
{
  if      (static_cast<unsigned int>(id) >= slist.size())
    return 1;
  else if (!(slist[id]))
    return 2;
  else {
    for (int v=0;v<MAX_DIM;v++) {
      slist[id]->position[v] = p[v];
			//
			// Also reset position in global data array (for data I/O)
			//
			SimPM.RS.sources[id].position[v] = p[v];
		}
	}
  
  return 0;
}

int RadiationSources::src_type(const int id, ///< src_id
				    	int *type ///< src type (output).
				     	)
{
  if (static_cast<unsigned int>(id) >= slist.size()) {
		rep.error("Requested source type for source which doesn't exist",id);
	}
	*type = slist[id]->type;
	return 0;
}

bool RadiationSources::src_at_infinity(const int id)
{
  if (static_cast<unsigned int>(id) >= slist.size()) {
		rep.error("Requested source type for source which doesn't exist",id);
	}
	return slist[id]->at_infinity;
}

bool RadiationSources::src_exists(const int id)
{
	int t = -99999;
	try {
		t = slist.at(id)->id;
	}
	catch (out_of_range) {
		return false;
	}
	//
	// If there was no out-of-range then the source exists.
	//
	cout <<"source id="<<t<<" exists?\n";
	return true;
}

//------------------------------------------------
// from dataio.cc:
    //
    // Now make sure RSP is properly populated.
    //
    //for (int i=0; i<SimPM.RS.Nsources; i++) {
    //  if (RSP.src_exists(i)) {
    //    cout <<"Source "<<i<<" already in RSP.  Ignoring it.\n";
    //  }
    //  else {
    //    RSP.add_src(SimPM.RS.sources.at(i).strength,
    //                SimPM.RS.sources.at(i).position,
    //                SimPM.RS.sources.at(i).type,
    //                SimPM.RS.sources.at(i).at_infinity
    //                );
    //  }
    //} // loop over sources
    
    //
    // Check we got all the sources:
    //
    //if (RSP.Nsources != RSP.nsources()) {
    //  rep.error("Got the wrong number of radiation sources!",
    //            RSP.Nsources-RSP.nsources());
    //}


