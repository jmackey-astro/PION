/// file_status.h
///
///  - 2010-02-03 JM: Added virtual destructor for file_status so
///     Intel compiler will stop complaining.
///
#ifndef FILE_STATUS_H
#define FILE_STATUS_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <dirent.h>
#include <errno.h>
#include <list>
#include <string>

/** \brief Some functions for querying and locking files.
 *
 * Note the acquire_lock() function has synchronization problems, and
 * needs to be written in a better way to get an instantaneous lock.  As
 * it is, two processes can (and do) lock the same file and both think
 * they are the only processes accessing it.
 * */
class file_status {
  public:
    virtual ~file_status() {}

    /** \brief returns true if a file exists, and false if not. */
    bool file_exists(std::string  ///< filename.
    );
    /** \brief Gain control of the file; if file is already locked then
     * wait until it is unlocked, then lock it and return.
     *
     * WARNING This function should be regarded as broken.
     * */
    int acquire_lock(std::string  ///< filename.
    );
    /** \brief If this process locked the file, then this function will
     * release the lock.
     * */
    int release_lock(std::string  ///< filename.
    );
    /** \brief Check if the file is locked.*/
    bool file_is_locked(std::string  ///< filename.
    );

    ///
    /// Given a directory and a string to match files to (may be blank),
    /// return sorted list of files.
    ///
    virtual int get_files_in_dir(
        const std::string,       ///< directory to list.
        const std::string,       ///< string that files start with
        std::list<std::string>*  ///< list to put filenames in.
    );

  protected:
    ///
    /// Return list of files in a given directory.
    ///
    virtual int get_dir_listing(
        const std::string,       ///< directory to list.
        std::list<std::string>*  ///< list to put filenames in.
    );

  private:
    /** \brief Lock file -- this should happen synchronously with having
     * determined the file is unlocked, but currently it doesn't.
     * */
    int file_lock(std::string  ///< filename.
    );
    /** \brief If I have locked the file, then unlock it. */
    void file_unlock(std::string  ///< filename.
    );
};

#endif  // FILE_STATUS_H
