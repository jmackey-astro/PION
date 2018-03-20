/// \file comm_files.cc
///
/// \brief Contains comms class for multi-process communication using 
/// files (no MPI!).
///
/// \author Jonathan Mackey
/// \date 2009-01-27.
/// 
/// Modifications:\n
///  - 2010.11.15 JM: replaced endl with c-style newline chars.
/// - 2012.05.15 JM: Added function for global-reduce (max/min/sum) of arrays,
///    but it is not yet implemented.  If I ever need it, I will write it...
/// - 2015.01.26 JM: Removed mpiPM (no longer global), added COMM
///    setup, and added get_rank_nproc() function

#ifdef PARALLEL
#ifdef USE_FILE_COMMS



#ifdef INTEL
#include <mathimf.h
#else
#include <cmath>
#endif // INTEL
#include <sstream>
using namespace std;


class comms_base *COMM = new comm_files ();


// ##################################################################
// ##################################################################



comm_files::comm_files()
{
  cout <<"*** comm_files constructor. ***\n";
  //
  // standard name for all the files:
  //
  comm_files::dir = "filecomms/";
  comm_files::ini = "INITIALISED";
  comm_files::bar = "_BARRIER_";
  go_order        = "master_says_go";
  go_response     = "_going";
  comm_global     = "GLOBAL_";
  celldata        = "cell_data_step_";
  doubledata      = "dubl_data_step_";
  finished_with_file = "_FINISHED";
  myrank = -1;
  nproc = -1;
  barrier_count=0;
}



// ##################################################################
// ##################################################################



comm_files::~comm_files()
{
  cout <<"*** comm_files  destructor. ***\n";
}



// ##################################################################
// ##################################################################



int comm_files::init(int *argc,   ///< number of program arguments.
		     char ***argv ///< character list of arguments.
		     )
{
  //
  // Get rank and nproc from command-line args.
  //
  string args[(*argc)];
  for (int i=0;i<(*argc);i++) {
    args[i] = (*argv)[i];
  }

  myrank = nproc = -1;
  for (int i=0;i<(*argc); i++) {
    if      (args[i].find("myrank=") != string::npos) {
      cout <<"found myrank as "<<i<<"th arg.\n";
      string rank= args[i].substr(7);
      myrank = atoi(rank.c_str());
    }
    else if (args[i].find("nproc=") != string::npos) {
      cout <<"found nproc as "<<i<<"th arg.\n";
      string np= args[i].substr(6);
      nproc = atoi(np.c_str());
    }
  }
  if (myrank<0 || nproc<0)
    rep.error("failed to find rank and nproc in command-line args",myrank);

  //
  // reduce number of args by 2, since myrank and nproc must be last two args.
  //
  (*argc) -= 2;

  //
  // If I'm the root processor, create the directory, and then a file in it
  //
  if (myrank==0) {
    //
    // create a directory for all the files we will create and destroy
    //
    string tmp="mkdir "+dir;
    system(tmp.c_str());
    //
    // create a file to let everyone know that I am present
    //
    string s(dir); s += ini;
    ofstream outf(s.c_str());
    outf.close();
  }
    
  else {
    //
    // wait for initialised file...
    //
    string s(dir); s += ini;
    wait_for_file(s);
  }

  ostringstream oss;
  oss<<dir<<ini<<"_rank_"<<myrank<<"_of_"<<nproc<<"_present_and_correct";
  ofstream outf(oss.str().c_str());
  outf.close();
  //
  // Now make sure everyone exists and is ready to go:
  //
  barrier("comm_files_init_end");

  cout << "comm_files::init():  rank: "<<myrank<<" nproc: "<<nproc<<"\n";
  return 0;
}


// ##################################################################
// ##################################################################



int comm_mpi::get_rank_nproc(
        int *r, ///< rank.
	int *n  ///< nproc
	)
{
  *r = comm_files::myrank;
  *n = comm_files::nproc;
  return err;
}



// ##################################################################
// ##################################################################



int comm_files::finalise()
{
  //
  // In case stuff didn't get initialised right...
  //
  if (myrank <0 || nproc<0) return 1;

  //
  // Otherwise set up a barrier and wait.
  //
  barrier("comm_files_finalise");
  if (myrank==0) {
    cout <<"comm_files::finalise() removing _initialised_ file\n";
    string s(dir); s += ini;
    remove(s.c_str());
  }
  cout << "comm_files::finalise():  rank: "<<myrank<<" nproc: "<<nproc<<"\n";
  return 0;
}



// ##################################################################
// ##################################################################



int comm_files::abort()
{
  ostringstream s; s<<dir<<"ERROR_ABORTABORTABORT";
  ofstream outf(s.str().c_str());
  outf.close();
  return 0;
}



// ##################################################################
// ##################################################################



int comm_files::barrier(const std::string msg)
{
#ifdef TESTING
  cout <<"comm_files::barrier(): "<<msg<<"\n";
#endif
  //
  // Need a counter on the barrier name, so that processes don't jump two barriers
  // while the master is waiting for the others.
  //
  comm_files::barrier_count ++;
  ostringstream tmp; tmp <<dir<<"rank_"<< myrank << bar << barrier_count << msg;
  ofstream outf(tmp.str().c_str());
  outf.close();

  //
  // Need to identify barriers uniquely -- one way is to have a counter for
  // what barrier we are on.
  //
  ostringstream bname; bname << bar << barrier_count << msg;

  if (myrank==0) {
    master_wait_on_slave_files("barrier",bname.str());
    master_send_go_signal(bname.str());
  }
  else {
    slave_recv_go_signal(bname.str());
  }

  remove(tmp.str().c_str());
#ifdef TESTING
  cout <<"comm_files::barrier(): barrier crossed: "<<msg<<"\n";
#endif
  return 0;
}



// ##################################################################
// ##################################################################



double comm_files::global_operation_double(const string s, ///< Either Max or Min
					   const double d  ///< this process's max/min value.
					   )
{
  double local=d, global=0.0;

  //
  // Every process writes their value to file.
  //
  ostringstream f; f<<dir<<"rank_"<<myrank<<"_"<<s;
  fs.acquire_lock(f.str());
  ofstream outf(f.str().c_str(),ios_base::binary);
  outf.write(reinterpret_cast<char *>(&local),sizeof(double));
  outf.close();
  fs.release_lock(f.str());

  //
  // Root process reads all these in and does the calculation.
  //
  if (myrank==0) {
    //
    // Read in all values:
    //
    double vals[nproc];
    for (int r=0;r<nproc;r++) {
      f.str(""); f<<dir<<"rank_"<<r<<"_"<<s;
      wait_for_file(f.str());
      fs.acquire_lock(f.str());
      ifstream infile(f.str().c_str(), ios_base::binary);
      //infile >> vals[r];
      infile.read(reinterpret_cast<char *>(&(vals[r])),sizeof(double));
#ifdef TESTING
      cout <<"global_operation_double() "<<s<<" got value "<<vals[r]<<" from file "<<f.str()<<"\n";
#endif
      infile.close();
      fs.release_lock(f.str());
      remove(f.str().c_str());
    }
    
    //
    // Do the operation on all the values:
    //
    if      (s=="MAX") {
      global = -HUGEVALUE;
      for (int r=0;r<nproc;r++) global = max(global,vals[r]);
    }
    else if (s=="MIN") {
      global = HUGEVALUE;
      for (int r=0;r<nproc;r++) global = min(global,vals[r]);
    }
    else if (s=="SUM") {
      global=0.0;
      for (int r=0;r<nproc;r++) global += vals[r];
    }
    else 
      rep.error("comm_files:global_operation_double: Bad identifier",s);
  
    //
    // Write the global value to file.
    //
    f.str(""); f<<dir<<comm_global<<barrier_count<<s;
    fs.acquire_lock(f.str());
    ofstream outfile(f.str().c_str(),ios_base::binary);
    //outfile << global;
    outfile.write(reinterpret_cast<char *>(&global),sizeof(double));
    outfile.close();
    fs.release_lock(f.str());
#ifdef TESTING
    cout <<"global_operation_double() "<<s<<" : global value = "<<global<<"\n";
#endif
  } // if myrank==0

  else {
    //
    // wait for global value to be written by root proc:
    //
    f.str(""); f<<dir<<comm_global<<barrier_count<<s;
    wait_for_file(f.str());
    if (fs.file_is_locked(f.str())) {
      usleep(FDELAY_USECS);
      while (fs.file_is_locked(f.str())) {
	usleep(FDELAY_USECS);
      }
    }
    ifstream infile(f.str().c_str(), ios_base::binary);
    if (!infile.is_open()) rep.error("failed to open infile for global max,min",f.str());
    //infile >>global;
    infile.read(reinterpret_cast<char *>(&global),sizeof(double));
    infile.close();
#ifdef TESTING
    cout <<"global_operation_double() "<<s<<" : read global value = "<<global<<"\n";
#endif
  } // not root proc.


  barrier("global_operation_double__end");
  if (myrank==0) {
    remove(f.str().c_str());
  }
  
  return global;
}



// ##################################################################
// ##################################################################


// ##################################################################
// ##################################################################



void comm_files::global_op_double_array(
                      const std::string s, ///< MAX,MIN,SUM
                      const size_t  Nel,   ///< Number of elements in array.
		      double *data         ///< pointer to this process's data array.
		      )
{
  rep.error("comm_files:global_op_double_array: Not implemented in file-comms",);
  return;
}




// ##################################################################
// ##################################################################





int comm_files::broadcast_data(const int sender,       ///< rank of sender.
			       const std::string type, ///< Type of data INT,DOUBLE,etc.
			       const int n_el,         ///< number of elements
			       void *data              ///< pointer to data.
			       )
{
  int err=0;
  string fn="broadcast_";
  ostringstream f;
  f<<dir<<comm_global<<fn<<type;

  if (myrank==0) {
    //
    // Open global file, cast pointer to correct type, and write data to file.
    //
    fs.acquire_lock(f.str());
    ofstream outfile(f.str().c_str(), ios_base::binary);
    if      (type=="DOUBLE") {
      double *d = static_cast<double *>(data);
      for (int i=0; i<n_el; i++) 
	outfile.write(reinterpret_cast<char *>(&(d[i])),sizeof(double));
      //outfile<<d[i];
    }
    else if (type=="INT") {
      int *d = static_cast<int *>(data);
      for (int i=0; i<n_el; i++)
	outfile.write(reinterpret_cast<char *>(&(d[i])),sizeof(int));
      //outfile<<d[i];
    }
    else if (type=="FLOAT") {
      float *d = static_cast<float *>(data);
      for (int i=0; i<n_el; i++)
	outfile.write(reinterpret_cast<char *>(&(d[i])),sizeof(float));
      //outfile<<d[i];
    }
    else rep.error("Bad type of data to send",type);
    outfile.close();
    fs.release_lock(f.str());

  }
  else {
    //
    // wait for file.
    //
    wait_for_file(f.str());
    fs.acquire_lock(f.str());
    ifstream infile(f.str().c_str(), ios_base::binary);
#ifdef TESTING
    cout <<"broadcast_data() "<<type<<": received data: [";
#endif
    if      (type=="DOUBLE") {
      double *d = static_cast<double *>(data);
      for (int i=0; i<n_el; i++) 
	infile.read( reinterpret_cast<char *>(&(d[i])),sizeof(double));
      //infile>>d[i];
#ifdef TESTING
      for (int i=0; i<n_el; i++) cout << d[i]<<", ";
#endif
    }
    else if (type=="INT") {
      int *d = static_cast<int *>(data);
      for (int i=0; i<n_el; i++)
	infile.read( reinterpret_cast<char *>(&(d[i])),sizeof(int));
      //infile>>d[i];
#ifdef TESTING
      for (int i=0; i<n_el; i++) cout << d[i]<<", ";
#endif
    }
    else if (type=="FLOAT") {
      float *d = static_cast<float *>(data);
      for (int i=0; i<n_el; i++)
	infile.read( reinterpret_cast<char *>(&(d[i])),sizeof(float));
      //infile>>d[i];
#ifdef TESTING
      for (int i=0; i<n_el; i++) cout << d[i]<<", ";
#endif
    }
    else rep.error("Bad type of data to send",type);
    infile.close(); cout <<"] ... "<<n_el<<" elements\n";
    fs.release_lock(f.str());
  }
  
  barrier("comm_files__broadcast_data__end");
  if (myrank==0) {
    remove(f.str().c_str());
  }

  return err;
}


int comm_files::send_cell_data(
      const int to_rank,    ///< rank to send to.
      std::list<cell *> *l, ///< list of cells to get data from.
      long int nc,          ///< number of cells in list (extra checking!)
      const int ndim, ///< ndim
      const int nvar, ///< nvar
      string &id,           ///< identifier for send, for tracking delivery later.
      const int comm_tag    ///< comm_tag, to say what kind of send this is.
      )
{
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::send_cell_data() starting. \n";
#endif //TESTING
  //
  // First initialise everything and check we have cells to get data from.
  //
  if(!id.empty()) id.erase();
  if (nc==0 || (l->empty()) ) {
    cout <<myrank<<"\t Nothing to send to rank: "<<to_rank<<" !!!\n";
    return 1;
  }
  if (to_rank<0 || to_rank>nproc)
    rep.error("to_rank is out of bounds",to_rank);
  list<cell *>::iterator c=l->begin();
  int err=0;

  //
  // Determine size of send buffer needed
  //
  int unitsize = 
    ndim*sizeof(CI.get_ipos(*c,0))
      + nvar*sizeof((*c)->P[0])
	+ sizeof((*c)->id);
  long int totalsize = 0;
  totalsize = sizeof(int) + nc*unitsize;

  //
  // Allocate memory for the record of the send.
  //
  struct sent_info *si = 0;
#ifdef TESTING
  si = mem.myalloc(si,1, "comm_files:send_cell_data: si");
#else
  si = mem.myalloc(si,1);
#endif //TESTING


  //
  // ALL SAME AS MPI VERSION UP TO HERE, NOW WE PACK+SEND DATA DIFFERENTLY IN WHAT FOLLOWS
  //
  // write data to file.  Format = [comm_tag,n_cells, totalsize(bytes), ALL_CELLS[c->id,c->x,c->Ph]]
  //
  int ct=0;
  ostringstream f;
  f<<dir<<celldata<<SimPM.timestep<<"_rank_"<<myrank<<"_to_rank_"<<to_rank;
  //
  // Create record of the send.
  //
  si->comm_tag = comm_tag;
  si->from_rank = myrank;
  si->to_rank = to_rank;
  //si->data = 0;
  si->type = COMM_CELLDATA;
  id = f.str();
  si->id = id;
  comm_files::sent_list.push_back(si);  

  fs.acquire_lock(f.str());
  ofstream outfile(f.str().c_str(), ios_base::binary);
  outfile.write(reinterpret_cast<char *>(&(si->comm_tag)),sizeof(int));
  outfile.write(reinterpret_cast<char *>(&(nc))          ,sizeof(long int));
  outfile.write(reinterpret_cast<char *>(&(totalsize))   ,sizeof(long int));
  int ipos[MAX_DIM];
  do {
    CI.get_ipos(*c,ipos);
    outfile.write(reinterpret_cast<char *>(&((*c)->id)),sizeof(int));
    for (int i=0;i<ndim;i++) 
      outfile.write(reinterpret_cast<char *>(&(ipos[i])),sizeof(int));
    for (int v=0;v<nvar;v++)
      outfile.write(reinterpret_cast<char *>(&((*c)->Ph[v])),sizeof(double));
    ct++;
    ++c;  // next cell in list.
  } while ( c != l->end() );

  //
  // Check that we packed the right amount of data:
  //
  //  cout <<myrank<<"\tcomm_pack_send_data: bufsiz: ";
  //  cout <<totalsize<<"  nc="<<nc<<" ct="<<ct<<"\n";
  if (ct != nc) rep.error("Length of list doesn't match nc",ct-nc);
  if (err) rep.error("MPI_Pack returned abnormally",err);


  //
  // Close file and release lock
  //
  outfile.close();
  fs.release_lock(f.str());
  f.str("");
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::send_cell_data() comm_tag="<<comm_tag<<" nc="<<nc<<" \n";
  cout <<"rank: "<<myrank<<"  comm_files::send_cell_data() returning.\n";
#endif //TESTING
  return 0;
}


int comm_files::wait_for_send_to_finish(string &id ///< identifier for the send we are waiting on.
				      )
{
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::wait_for_send_to_finish() starting\n";
#endif //TESTING

  //
  // Find the send in the list of active sends, based on the identifier string.
  //
  int el=0;
  list<sent_info *>::iterator i;
  struct sent_info *si=0;
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::wait_for_send_to_finish() more than one send, so finding in list.\n";
  cout <<"\t\tsend id="<<id<<"\n";
#endif //TESTING

  for (i=sent_list.begin(); i!=sent_list.end(); ++i) {
    si = (*i); el++;
    if (si->id==id) break;
  }

  if (i==sent_list.end())
    rep.error("Failed to find send with id:",id);
#ifdef TESTING
  cout <<"found send id="<<si->id<<" and looking for id="<<id<<"\n";
  cout <<"rank: "<<myrank<<"  comm_files::wait_for_send_to_finish() found this send.\n";
#endif //TESTING

  //
  // Now we have the record of the send, so we wait for receiver to finish with the file.
  //
  wait_for_peer_to_read_file(si->id);
  remove((si->id).c_str());
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::wait_for_send_to_finish() peer has read file, and I deleted it.\n";
#endif

#ifdef TESTING
  si = mem.myfree(si,"comm_files::wait_for_send_to_finish(): si");
#else
  si = mem.myfree(si);
#endif
  sent_list.erase(i);

#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::wait_for_send_to_finish() returning\n";
#endif //TESTING
  return 0;
}


int comm_files::look_for_data_to_receive(int *from_rank, ///< rank of sender
					 string &id,     ///< identifier for receive.
					 int *comm_tag,  ///< comm_tag associated with data.
					 const int type  ///< type of data we are looking for.
					 )
{
  int err=0;
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::look_for_data_to_receive() starting\n";
#endif //TESTING
  //
  // Create a new received info record.
  //
  //
  struct recv_info *ri=0;
#ifdef TESTING
  ri = mem.myalloc(ri,1, "comm_files:look_for_data_to_receive: ri");
#else
  ri = mem.myalloc(ri,1);
#endif
  ri->to_rank   = myrank;
  if (type!=COMM_CELLDATA && type!=COMM_DOUBLEDATA)
    rep.error("only know two types of data to look for!",type);
  ri->type = type;

  //
  // Now look for data being sent to us:
  // This is a very inefficient method, but speed it not of the essence here...
  //
#ifdef TESTING
  cout <<"rank: "<<myrank<<"  comm_files::look_for_data_to_receive() looking for source\n";
#endif //TESTING
  ostringstream f; bool found=false;
  do {
    usleep(FDELAY_USECS);
    usleep(FDELAY_USECS);
    for (int rank=0; rank<nproc; rank++) {
      if (!found) {
	if (type==COMM_CELLDATA) {
	  f.str(""); f<<dir<<celldata  <<SimPM.timestep<<"_rank_"<<rank<<"_to_rank_"<<myrank;
	}
	else if (type==COMM_DOUBLEDATA) {
	  f.str(""); f<<dir<<doubledata<<SimPM.timestep<<"_rank_"<<rank<<"_to_rank_"<<myrank;
	}
	if (fs.file_exists(f.str()) && !fs.file_exists(f.str()+finished_with_file)) {
	  //
	  // Need to check for 'finished' file (in if statement above) in case there 
	  // is a queue of data and we could re-find the data we just finished reading.
	  //
#ifdef TESTING
	  cout <<"comm_files::look_for_data_to_receive: found data file : "<< f.str()<<"\n";
#endif
	  found=true;
	  *from_rank = rank;
	  id = f.str();
	}
      } // if (!found) look for file.
    } // loop over all ranks
  } while (!found);

  //
  // Wait for file to be ready, and read comm_tag from start of file.
  //
  wait_for_file(f.str());
  fs.acquire_lock(f.str());
  ifstream infile(f.str().c_str(), ios_base::binary);
  if (!infile.is_open()) rep.error("failed to open file for reading",f.str());
  infile.read( reinterpret_cast<char *>(comm_tag),sizeof(int));

#ifdef TESTING
  cout <<"comm_files::look_for_data_to_receive: got comm_tag: "<<*comm_tag<<" from file.\n";
#endif
  infile.close();
  fs.release_lock(f.str());

  ri->id     = id;
  ri->comm_tag  = *comm_tag;
  ri->from_rank = *from_rank;
  comm_files::recv_list.push_back(ri);

#ifdef TESTING
  cout <<"comm_files::look_for_data_to_receive: returning.\n";
#endif //TESTING
  return err;
}

int comm_files::receive_cell_data(
      const int from_rank,  ///< rank of process we are receiving from.
      std::list<cell *> *l, ///< list of cells to get data for. 
      const long int nc,    ///< number of cells in list (extra checking!)
      const int ndim, ///< ndim
      const int nvar, ///< nvar
      const int comm_tag,   ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
      const string &id      ///< identifier for receive, for any book-keeping that might be needed.
      )
{
  //  int err=0;
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: starting.\n";
#endif //TESTING

  //
  // Find the recv in the list of active receives.  I use a string identifier for
  // sends and receives, and the look_for_data() function returns a string which
  // should be passed to this function.
  //
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: recv_list size="<<recv_list.size()<<"\n";
#endif //TESTING
  if (recv_list.empty()) rep.error("Call look4data before receive_data",recv_list.size());

  struct recv_info *info=0;
  list<recv_info *>::iterator i;
  for (i=recv_list.begin(); i!=recv_list.end(); ++i) {
    info = (*i);
    if (info->id==id) break;
  }
  if (i==recv_list.end())
    rep.error("Failed to find recv with id:",id);
#ifdef TESTING
  cout <<"found recv id="<<info->id<<" and looking for id="<<id<<"\n";
  cout <<"comm_files::receive_cell_data: found recv id\n";
#endif //TESTING

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank) rep.error("from_ranks don't match",from_rank-info->from_rank);
  if (comm_tag  != info->comm_tag ) rep.error("Bad comm_tag",comm_tag);
  if (id        != info->id       ) rep.error("Bad id",id);

  //
  // DIFFERENT FROM MPI COMMS BELOW THIS POINT
  //
  // Open file:
  //
  fs.acquire_lock(info->id);
  ifstream infile((info->id).c_str(), ios_base::binary);
  if (!infile.is_open()) rep.error("failed to open file for reading",info->id);

  //
  // Skip comm_tag, read in n_cells, totalsize, and then start reading data into cells.
  // Format = [comm_tag, n_cells, totalsize(bytes), ALL_CELLS[c->id,c->x,c->Ph]]
  //
  int tmp=0; int tmp2=0;
  long int n_cells=0, totalsize=0;
  infile.read( reinterpret_cast<char *>(&(tmp))      ,sizeof(int));
  infile.read( reinterpret_cast<char *>(&(n_cells))  ,sizeof(long int));
  infile.read( reinterpret_cast<char *>(&(totalsize)),sizeof(long int));
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: got comm:"<<tmp<<" n_cells="<<n_cells<<" size="<<totalsize<<"\n";
#endif //TESTING
  if (n_cells != nc)
    rep.error("comm_files::receive_cell_data: n_cells has unexpected value",n_cells-nc);

  //
  // For data, i send id,position,state_vec, but id and position aren't used, so I just
  // read them into temp arrays here.  i used them for debugging to make sure the cells
  // were ordered the same way by sender and receiver.  I should really get rid of them
  // now but haven't bothered.
  //
  list<cell*>::iterator c=l->begin();
  int ct=0;
  do {
    if (c==l->end()) rep.error("Got too many cells!",ct);
    infile.read( reinterpret_cast<char *>(&(tmp))      ,sizeof(int));
    for (int i=0;i<ndim;i++)
      infile.read( reinterpret_cast<char *>(&(tmp2))      ,sizeof(int));
    for (int v=0;v<nvar;v++)
      infile.read( reinterpret_cast<char *>(&((*c)->Ph[v]))      ,sizeof(double));
    ++c;
    ct++;
  } while (c!=l->end());
  if (ct != nc) rep.error("got wrong number of cells", ct-nc);


  //
  // Close file, release lock, and create 'finished' file
  //
  infile.close();
  fs.release_lock(info->id);

  ostringstream f; f<<info->id<<finished_with_file;
  ofstream outf(f.str().c_str());
  outf.close(); f.str("");
  //
  // DIFFERENT FROM MPI COMMS ABOVE HERE
  //

  //
  // Free memory and delete entry from recv_list
  //
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: freeing memory\n";
#endif //TESTING
  if (recv_list.size()==1) {
    recv_list.pop_front();
  }
  else
    rep.error("recv list is big!",recv_list.size());
#ifdef TESTING
  info= mem.myfree(info,"comm_files::receive_cell_data() info");
#else
  info= mem.myfree(info);
#endif //TESTING

#ifdef TESTING
  cout <<"comm_files::receive_cell_data: returning\n";
#endif //TESTING
  return 0;
}


int comm_files::send_double_data(
        const int to_rank,   ///< rank to send to.
	const long int n_el, ///< size of buffer, in number of doubles.
	const double *data,  ///< pointer to double array.
	string &id,          ///< identifier for send, for tracking delivery later.
	const int comm_tag   ///< comm_tag, to say what kind of send this is.
				 )
{
  if (!data) rep.error("comm_files::send_double_data() null pointer!",data);

  //
  // Allocate memory for a record of the send
  //
  struct sent_info *si=0;
#ifdef TESTING
  si = mem.myalloc(si,1, "comm_files:send_double_data: si");
#else
  si = mem.myalloc(si,1);
#endif //TESTING

  //
  // filename for send: I tag with timestep just to be sure.  Set string identifier
  // to be the filename.  This id is used by the caller to track the send.
  //
  ostringstream f;
  f<<dir<<doubledata<<SimPM.timestep<<"_rank_"<<myrank<<"_to_rank_"<<to_rank;
  id = f.str();

  //
  // Add info to send record.
  //
  si->comm_tag  = comm_tag;
  si->from_rank = myrank;
  si->to_rank   = to_rank;
  si->type      = COMM_DOUBLEDATA;
  si->id        = id;
  comm_files::sent_list.push_back(si);

  //
  // open file and write data to file:
  // I cast some of the data to (const char *) because it is const data.
  // Doesn't seem to make any difference.
  //
  fs.acquire_lock(id);
  ofstream outfile(id.c_str(), ios_base::binary);
  outfile.write(reinterpret_cast<char *>(&(si->comm_tag)),sizeof(int));
  outfile.write(reinterpret_cast<const char *>(&(n_el)),sizeof(long int));
  //outfile << comm_tag << n_el;
  for (int i=0;i<n_el;i++) 
    outfile.write(reinterpret_cast<const char *>(&(data[i])),sizeof(double));
  //outfile << data[i];
  outfile.close();
  fs.release_lock(id);

  //
  // All done, so return.
  //
  return 0;
}

int comm_files::receive_double_data(const int from_rank, ///< rank of process we are receiving from.
				    const int comm_tag,  ///< comm_tag: what sort of comm we are looking for (PER,MPI,etc.)
				    const string &id,    ///< identifier for receive, for any book-keeping that might be needed.
				    const long int nel,  ///< number of doubles to receive
				    double *data         ///< Pointer to array to write to (must be already initialised).
				    )
{
#ifdef TESTING
  cout <<"comm_files::receive_double_data: starting.\n";
#endif //TESTING

  //
  // Find the recv in the list of active receives, based on the id string passed to
  // the function.  This should have been obtained by the look_for_data() function.
  //
#ifdef TESTING
  cout <<"comm_files::receive_double_data: recv_list size="<<recv_list.size()<<"\n";
#endif //TESTING
  if (recv_list.empty()) rep.error("Call look4data before receive_data",recv_list.size());

  struct recv_info *info=0;
  list<recv_info *>::iterator i;
  for (i=recv_list.begin(); i!=recv_list.end(); ++i) {
    info = (*i);
    if (info->id==id) break;
  }
  if (i==recv_list.end())
    rep.error("Failed to find recv with id:",id);
#ifdef TESTING
  cout <<"found recv id="<<info->id<<" and looking for id="<<id<<"\n";
  cout <<"comm_files::receive_double_data: found recv id\n";
#endif //TESTING

  //
  // Check that everthing matches.
  //
  if (from_rank != info->from_rank) rep.error("from_ranks don't match",from_rank-info->from_rank);
  if (comm_tag  != info->comm_tag ) rep.error("Bad comm_tag",comm_tag);
  if (id        != info->id       ) rep.error("Bad id",id);
  if (info->type!= COMM_DOUBLEDATA) rep.error("data is not double array!",info->type);

  //
  // DIFFERENT FROM MPI COMMS BELOW HERE
  //
  // Open file:
  //
  fs.acquire_lock(info->id);
  ifstream infile((info->id).c_str(), ios_base::binary);

  //
  // Skip comm_tag, read in n_el, and then start reading data into array.
  // Format = [comm_tag,n_el, data[]]
  //
  int tmp=0; long int ct=0;
  infile.read( reinterpret_cast<char *>(&(tmp)),sizeof(int));
  infile.read( reinterpret_cast<char *>(&(ct)) ,sizeof(long int));
  //  infile >> tmp >> ct;
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: got comm:"<<tmp<<" n_el="<<ct<<"\n";
#endif //TESTING

  if (nel != ct)
    rep.error("comm_files::receive_cell_data: n_el has unexpected value",nel-ct);
  for (int v=0;v<nel;v++)
    infile.read( reinterpret_cast<char *>(&(data[v])),sizeof(double));
    //infile >> data[v];

  //
  // Close file, release lock, and create 'finished' file
  //
  infile.close();
  fs.release_lock(info->id);
  ostringstream f; f<<info->id<<finished_with_file;
  ofstream outf(f.str().c_str());
  outf.close(); f.str("");

  //
  // DIFFERENT FROM MPI COMMS ABOVE HERE
  //

  //
  // Free memory and delete entry from recv_list
  //
#ifdef TESTING
  cout <<"comm_files::receive_cell_data: freeing memory\n";
#endif //TESTING
  if (recv_list.size()==1) {
    recv_list.pop_front();
  }
  else
    rep.error("recv list is big!",recv_list.size());
#ifdef TESTING
  info= mem.myfree(info,"comm_files::receive_cell_data() info");
#else
  info= mem.myfree(info);
#endif //TESTING

  return 0;
}


#ifdef SILO
int comm_files::silo_pllel_init(const int n_files,    ///< number of files to write
				const std::string iotype, ///< READ or WRITE
				const std::string session_id, ///< identifier for this read/write.
				int *group_rank,    ///< rank of group (nth file).
				int *rank_in_group  ///< rank in group (nth proc in file i).
				)
{
  comm_files::silo_id   = session_id;
  comm_files::num_files = n_files;
  comm_files::silo_filetype = DB_PDB;
  comm_files::silo_iotype = iotype;

  //
  // Silo's parallel interface writes nf files by np processes.  So I need 
  // to create nf groups of processes, and assign each process a rank in its
  // group.  In principle this code is very easy to modify to work for any file
  // I/O method, if I wrote file open/close/create functions that could be
  // redefined.  The wait_for_file() and finish_with_file() functions
  // could take void ** pointers instead of DBfile **, so that it really 
  // is a general interface.
  //
  int ngrp=nproc/num_files; // integer division
  if (!GS.equalD(static_cast<double>(ngrp), static_cast<double>(nproc)/num_files))
    rep.error("comm_files::silo_pllel_init() choose number of files to divide into nproc evenly!",
	      static_cast<double>(nproc)/num_files);
  
  *group_rank    = myrank/ngrp;
  *rank_in_group = myrank % ngrp;
  comm_files::grp_rank    = *group_rank;
  comm_files::rank_in_grp = *rank_in_group;
#ifdef TESTING
  cout <<"comm_files::silo_pllel_init() grp_rank="<<grp_rank<<" and rank_in_grp="<<rank_in_grp<<"\n";
#endif
  barrier("silo_pllel_init__end");
  return 0;
}

int comm_files::silo_pllel_wait_for_file(const std::string s_id, ///< identifier for this read/write.
					 const std::string s_filename, ///< File Name
					 const std::string s_dir, ///< Directory to open in file.
					 DBfile **dbfile ///< pointer that file gets returned in.
					 )
{
#ifdef TESTING
  cout <<"comm_files::silo_pllel_wait_for_file() opening file: "<<s_filename<<" into directory: "<<s_dir<<"\n";
#endif
  if (s_id!=silo_id) rep.error(s_id, silo_id);
  if (*dbfile) rep.error("please pass in null file pointer!",*dbfile);


  //
  // Silo's parallel interface has one process reading from/writing to each file at a time, so
  // the first process can start immediately, and the others have to wait for a signal that
  // the previous one has finished with the file.  All the read/write operations are then done
  // with no communication by other code.  The wait_for_file() and finish_with_file() functions
  // could take void ** pointers instead of DBfile **, so that it really is a general interface.
  //
  if (rank_in_grp==0) {
#ifdef TESTING
    cout <<"comm_files::silo_pllel_wait_for_file() i'm first in group, creating/opening file and returning.\n";
#endif

    if      (silo_iotype=="READ") {
      //
      // Open file for reading.
      //
      *dbfile = DBOpen(s_filename.c_str(), silo_filetype, DB_READ);
      if (!(*dbfile)) rep.error("open silo file failed.",*dbfile);
      DBSetDir(*dbfile,"/");
      DBSetDir(*dbfile,s_dir.c_str());
    }

    else if (silo_iotype=="WRITE") {
      //
      // Create file for writing.
      //
      *dbfile = DBCreate(s_filename.c_str(), DB_CLOBBER, DB_LOCAL, "comment about data", silo_filetype);
      if (!(*dbfile)) rep.error("create silo file failed.",*dbfile);
      
      //
      // Create directory in file, and move to that directory.
      //
      DBSetDir(*dbfile,"/");
      DBMkDir(*dbfile,s_dir.c_str());
      DBSetDir(*dbfile,s_dir.c_str());
    }
    else rep.error("need either READ or WRITE as iotype",silo_iotype);
  }

  else {
    //
    // Wait for my turn to open file -- previous process will write a 'finished' file.
    //
    ostringstream f;
    f<<dir<<silo_id<<"_grp_"<<grp_rank<<"_rank_"<<rank_in_grp-1<<finished_with_file;
    wait_for_file(f.str());
    fs.acquire_lock(f.str());
#ifdef TESTING
    cout <<"comm_files::silo_pllel_wait_for_file(): found file: "<<f.str()<<", so my turn to write file.\n";
#endif
    remove(f.str().c_str());
    fs.release_lock(f.str());

    if      (silo_iotype=="READ") {
      //
      // Open file for reading.
      //
      *dbfile = DBOpen(s_filename.c_str(), silo_filetype, DB_READ);
      if (!(*dbfile)) rep.error("open silo file failed.",*dbfile);
      DBSetDir(*dbfile,"/");
      DBSetDir(*dbfile,s_dir.c_str());
    }
    else if (silo_iotype=="WRITE") {
      //
      // Open file for writing.
      //
      *dbfile = DBOpen(s_filename.c_str(), silo_filetype, DB_APPEND);
      if (!(*dbfile)) rep.error("open silo file failed.",*dbfile);
      DBSetDir(*dbfile,"/");
      DBMkDir(*dbfile,s_dir.c_str());
      DBSetDir(*dbfile,s_dir.c_str());
    }
    else rep.error("need either READ or WRITE as iotype",silo_iotype);
  }

  return 0;
}

int comm_files::silo_pllel_finish_with_file(const std::string s_id, ///< identifier for this read/write.
					    DBfile **dbfile ///< pointer to file we have been working on.
					    )
{
#ifdef TESTING
  cout <<"comm_files::silo_pllel_finish_with_file() passing file on to next proc.\n";
#endif
  if (s_id!=silo_id) rep.error(s_id, silo_id);
  if (!(*dbfile)) rep.error("file pointer i s null!",*dbfile);

  //
  // Last process in each group needs to close the file and return.
  // All others need to write a 'finished' message for the next process
  // in the group to read and take its turn with the file.
  //

  int ngrp=nproc/num_files; // integer division
  if (rank_in_grp==(ngrp-1) || myrank==(nproc-1)) {
    //
    // I am last file in group, so just close file and return
    //
    DBClose(*dbfile); *dbfile=0;
#ifdef TESTING
    cout <<"comm_files::silo_pllel_finish_with_file() I'm last proc in group, returning.\n";
#endif
  }
  else {
    //
    // Have someone to pass file onto, so close silo file and create handoff file.
    //
    DBClose(*dbfile); *dbfile=0;
    ostringstream f;
    f<<dir<<silo_id<<"_grp_"<<grp_rank<<"_rank_"<<rank_in_grp<<finished_with_file;
    fs.acquire_lock(f.str());
    ofstream outfile(f.str().c_str());
    outfile.close();
    fs.release_lock(f.str());
#ifdef TESTING
    cout <<"comm_files::silo_pllel_finish_with_file() closed file and created finished file, so returning.\n";
#endif
  }

  silo_id.erase();
  return 0;
}

void comm_files::silo_pllel_get_ranks(const std::string id, ///< identifier for this read/write.
				      const int proc, ///< processor rank
				      int *group, ///< rank of group processor is in.
				      int *rank   ///< rank of processor within group.
				      )
{
  if (id!=silo_id) rep.error(id, silo_id);

  int ngrp=nproc/num_files; // integer division
  *group = proc/ngrp;
  *rank  = proc % ngrp;
  return;
}
#endif // SILO


/*************************************************************************************/
/********************* PRIVATE FUNCTIONS -- PRIVATE FUNCTIONS ************************/
/*************************************************************************************/
void comm_files::wait_for_file(const string filename)
{
  //
  // Waits for filename to appear.
  //
  int  stop=0;
  do {
    usleep(FDELAY_USECS);
    if (!fs.file_exists(filename)) {
#ifdef TESTING
      cout <<"   rank "<<myrank<<": ...still waiting for file "<<filename<<"\n";
#endif
    }
    else {
#ifdef TESTING
      cout <<"   rank "<<myrank<<": recieved signal ("<<filename<<") - proceeding\n";
#endif
      stop=1;
    }
  } while (!stop);
  return;
}

void comm_files::wait_for_file_to_disappear(const string filename)
{
  //
  // Does exactly what it says on the tin.
  //
  int  stop=0;
  do {
    usleep(FDELAY_USECS);
    if (fs.file_exists(filename)) {
#ifdef TESTING
      cout <<"   rank "<<myrank<<": ...still waiting for file "<<filename<<" to be deleted\n";
#endif
    }
    else {
#ifdef TESTING
      cout <<"   rank "<<myrank<<": file "<<filename<<" has been deleted -- proceeding.\n";
#endif
      stop=1;
    }
  } while (!stop);
  return;
}


void comm_files::master_send_go_signal(const std::string identifier)
{
  //
  // Informs all other nodes that we are ready to proceed to the next stage.
  //
  //
  // Send the go signal
  //
  ostringstream tmp; tmp <<dir<<go_order<<identifier;
  fs.acquire_lock(tmp.str());
  ofstream outf(tmp.str().c_str());
  outf.close();
  fs.release_lock(tmp.str());
  //
  // Wait for replies (don't need to wait for rank 0 - that's us!)
  //
  string response = go_response+identifier;
  master_wait_on_slave_files("master_send_go_signal()",response);

  //
  // All slaves should have replied - clean up the files.
  //
  remove(tmp.str().c_str());
  for (int rank=1; rank<nproc; rank++) {
    tmp.str(""); tmp <<dir<<"rank_"<<rank<<response;
    //sprintf(filename,"filecomms/slave_%03d_going",rank);
    remove(tmp.str().c_str());
  }
  tmp.str("");
  return ;
}

void comm_files::slave_recv_go_signal(const std::string identifier)
{
  //
  // Recieves the go signal from the master node. (BLOODY WELL BETTER HOPE
  //  THAT IT DOESN'T PROCEED UNTIL IT GETS THE SIGNAL!!!!)
  //
  //DbgMsg("node %02i,    slave %03d: waiting for go signal\n",gd->M_rank,gd->M_rank);
  ostringstream tmp; tmp <<dir<<go_order<<identifier;
  wait_for_file(tmp.str());

  //DbgMsg("node %02i,    slave %03d: acknowledging signal\n",gd->M_rank,gd->M_rank);
  tmp.str(""); tmp<<dir<<"rank_"<<myrank<<go_response<<identifier;
  ofstream outf(tmp.str().c_str());
  outf.close();

  //DbgMsg("node %02i,    slave %03d: commencing graft\n",gd->M_rank,gd->M_rank);
  return;
}

void comm_files::master_wait_on_slave_files(const std::string msg,       ///< message.
					    const std::string identifier ///< what slave files are called.
					    )
{
#ifdef TESTING
  cout <<"master_wait_on_slave_files: "<<msg<<": starting.\n";
#endif
  ostringstream tmp;

  //
  // Don't need to wait for my own signal, so start at 1!
  //
  for (int rank=1;rank<nproc;rank++) {
    tmp.str(""); tmp<<dir<<"rank_"<<rank<<identifier;
    wait_for_file(tmp.str());
#ifdef TESTING
    cout <<"master_wait_on_slave_files:"<<msg<<" got file:"<<tmp.str()<<"\n";
#endif
  }

  tmp.str("");
  return;
}

void comm_files::wait_for_peer_to_read_file(const std::string f ///< name of file
					    )
{
  ostringstream tmp; tmp<<f<<finished_with_file;
  wait_for_file(tmp.str());
#ifdef TESTING
  cout <<"wait_for_peer_to_read_file: got file: "<<tmp.str()<<" ...deleting it\n";
#endif
  remove(tmp.str().c_str());
  return;
}


#endif // USE_FILE_COMMS
#endif //PARALLEL

