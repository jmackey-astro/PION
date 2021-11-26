///
/// \file mem_manage.h
/// \author Jonathan Mackey
///
/// Declares a class for managing memory, but I'm not sure how useful
/// it really is.
///
/// Modifications:
/// - 2015.01.08 JM: created file, moved class from global.h

#ifndef MEM_MANAGE_H
#define MEM_MANAGE_H

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include <new>


#include <spdlog/spdlog.h>

// ------------------ MEMORY MANAGEMENT -------------------
/** \brief This class contains generic functions for dynamically allocating
 * and freeing memory.  Hopefully everywhere in the code will use this and
 * not have to write it all out every time.
 */
class memory_management {
public:
  memory_management() {}
  ~memory_management() {}
  /** \brief Allocates memory for an array of any valid type.
   *
   * If pointer passed in is not null (=0), a warning is displayed and
   * the pointer is returned unchanged.  If it is null, then the
   * function attempts to initialise it, and if it fails, it exits the
   * code, since memory allocation failures are generally non-recoverable
   * (in the sense that if I ask for memory, I need it for something...).
   */
  ///
  /// Turns out this isn't particularly useful since I've read in various
  /// places that allocation requests NEVER fail unless they eat up
  /// all the RAM and SWAP, and maybe even more.  OS's are designed to
  /// at least try to deal with any memory requests.
  ///
  template<class T>
  T *myalloc(
      T *ptr,              ///< uninitialised null pointer.
      const long int n_el  ///< number elements to initialise.
  )
  {
    try {
      ptr = new T[n_el];
    }
    catch (const std::bad_alloc &) {
      spdlog::error("mem_alloc() pointer initialisation failed");
      spdlog::error(
          "{}: {}", "mem_alloc() pointer initialisation failed.",
          fmt::ptr(ptr));
    }
    return ptr;
  }

  /** \brief Delete a pointer to an array.
   *
   * If the pointer is already null (=0), print a warning and don't
   * try to delete it again!  Otherwise delete it and set it to zero.
   */
  template<class T>
  T *myfree(T *ptr)
  {
    //
    // if ptr is already null, can't free anything, so just return.
    //
    if (!ptr) {
    }
    else {
      delete[] ptr;
      ptr = 0;
    }
    return 0;
  }

  /** \brief Delete a pointer to single object (no [] in delete).
   *
   * If the pointer is already null (=0), print a warning and don't
   * try to delete it again!  Otherwise delete it and set it to zero.
   */
  template<class T>
  T *myfree_single(T *ptr)
  {
    //
    // if ptr is already null, can't free anything, so just return.
    //
    if (!ptr) {
    }
    else {
      delete ptr;
      ptr = 0;
    }
    return ptr;
  }
};

extern class memory_management mem;
/************************* MEMORY MANAGEMENT ***********************/
/*******************************************************************/

#endif  // MEM_MANAGE_H
