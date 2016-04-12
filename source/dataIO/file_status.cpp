/// \file file_status.cpp
///
/// - 2013-02-19 JM: hived off from dataio.cc to its own file.
/// - 2015.01.15 JM: Added new include statements for new PION version.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"
#include "tools/reporting.h"
#include "tools/mem_manage.h"
#ifdef TESTING
#include "tools/command_line_interface.h"
#endif // TESTING

#include "dataIO/file_status.h"

using namespace std;

//------------------------------------------------------
//----------- FILE_STATUS CLASS (LOCK/UNLOCK)-----------
//------------------------------------------------------

bool file_status::file_exists(string fname) 
{
  ifstream inf(fname.c_str());
  if (inf) {
//    cout <<":\t file exists.\n";
    inf.close();
    return true;
  }
  else return false;
}

bool file_status::file_is_locked(string fname) 
{
  fname += "_lock";
  ifstream inf(fname.c_str());
  if (inf) {
//    cout <<"file_status:\t file is locked.\n";
    inf.close();
    return true;
  }
  else return false;
}

int file_status::acquire_lock(string fname)
{
  struct timespec waittime; waittime.tv_sec = 1; waittime.tv_nsec = 50000000;
  while (file_is_locked(fname)) {
    nanosleep(&waittime,0);
//    cout <<"file_status:\t File locked... waiting\n";
  }
  if(file_lock(fname)) rep.error("Failed to lock file",1);
  return(0);
}
int file_status::release_lock(string fname)
{
  file_unlock(fname);
  return 0;
}

int file_status::file_lock(string fname)
{
  fname += "_lock";
  ofstream outf(fname.c_str());
  if(!outf) {
    return 1;
  }
  outf <<"locking file"<<"\n";
  outf.close();
//  cout <<"file_status:\t locking file.\n";
  return 0;  
}

void file_status::file_unlock(string fname) 
{
  fname += "_lock";
  ifstream inf(fname.c_str());
  if (!inf) {cout <<"file_status:\t File is not locked!\n"; return;}
  inf.close();
  if( remove(fname.c_str()) != 0 )
    rep.error("Error deleting lock file",remove(fname.c_str()));
  else
//    cout<<"file_status:\t Unlocked file.\n";
  return; 
}

//
// Return list of files in a given directory.
//
int file_status::get_dir_listing (
      const string dir,    ///< directory to list.
      list<string> *files  ///< list to put filenames in.
      )
{
  cout <<"get_dir_listing() reading directory: "<<dir<<"\n";
  DIR *dp=0;
  struct dirent *dirp=0;
  if((dp  = opendir(dir.c_str())) == 0) {
    cout << "Error(" << errno << ") opening " << dir << "\n";
    return errno;
  }
  
  while ((dirp = readdir(dp)) != 0) {
    //string temp=dirp->d_name;files->push_back(temp);
    files->push_back(string(dirp->d_name));
    //cout <<"\tget_dir_listing() file: "<<string(dirp->d_name)<<"\n";
  }
  closedir(dp);
  cout <<"get_dir_listing() done."<<"\n";
  return 0;
}

//
// Given a directory and a string to match files to (may be blank),
// return sorted list of files.
//
int file_status::get_files_in_dir(
      const string dir,    ///< directory to list.
      const string str,    ///< string that files start with
      list<string> *files  ///< list to put filenames in.
      )
{
  cout <<"get_files_in_dir(): starting.\n";
  int err=0;
  if (!files->empty())
    cout <<"WARNING: list of files is not empty, adding to end of list.\n";

  //
  // get directory listing
  //
  err += get_dir_listing(dir,files);
  cout <<"get_files_in_dir(): got dir listing with "<<files->size()<<" elements.\n";

  //
  // remove elements that don't begin with a given substring
  //
  if (!str.empty()) {
    cout <<"get_files_in_dir(): looking for substring in filenames: "<<str<<"\n";
    list<string>::iterator i=files->begin();
    if (i!=files->end()) { // check that dir listing is not empty...
      do {
        //cout <<*i<<"\n";
	//if ((*i).find(str) == string::npos) {
	//
	// If filename doesn't start with str, then delete it.
	//
	if ((*i).find(str) != 0) {
	  //cout <<"removing file "<<*i<<" from list.\n";
	  files->erase(i);
	  //files->remove(i);
	  i=files->begin();
	}
	else {
          i++;
          //cout <<"Keeping file "<<*i<<" in list.\n";
        }
      } while (i!=files->end());
    }
  }
  else
    cout <<"get_files_in_dir(): No substring, so not removing any elements. returning...\n";

  //
  // sort remaining elements.
  //
  files->sort();
  cout <<"get_files_in_dir(): done.\n";
  return err;
}


