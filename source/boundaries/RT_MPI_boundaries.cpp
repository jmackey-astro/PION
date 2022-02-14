/// \file RT_MPI_boundaries.cpp
/// \brief Class definitions for raytracing across a parallel
///   decomposed domain, by communicating with neighbouring MPI
///   processes
/// \author Jonathan Mackey
///
/// Modifications :\n
/// - 2018.08.09 JM: moved code from parallel uniform grid.

#include "defines/functionality_flags.h"
#include "defines/testing_flags.h"

#include "boundaries/RT_MPI_boundaries.h"
#include "tools/mem_manage.h"

#ifdef SPDLOG_FWD
#include <spdlog/fwd.h>
#endif
#include <spdlog/spdlog.h>
/* prevent clang-format reordering */
#include <fmt/ranges.h>

#include <sstream>
#include <vector>
using namespace std;

//#define RT_TESTING

// ##################################################################
// ##################################################################



RT_MPI_bc::~RT_MPI_bc()
{
#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc destructor.");
#endif
  //
  // Delete cells created for RT receive boundaries.
  // Have to delete lists of cells, and boundary_data, and maybe cells.
  //
  for (unsigned int isrc = 0; isrc < RT_source_list.size(); isrc++) {
    struct RT_source_comms_info *rc = &(RT_source_list[isrc]);
    if (!rc->RT_recv_list.empty()) {
      struct boundary_data *b;
      for (unsigned int ii = 0; ii < rc->RT_recv_list.size(); ii++) {
        int n = rc->RT_recv_list[ii].dir;
        b     = rc->RT_recv_list[ii].RT_bd;
        //
        // If a pre-existing boundary, just delete the list.
        //
        if (n == dir_XN || n == dir_XP || n == dir_YN || n == dir_YP
            || n == dir_ZN || n == dir_ZP) {
          // cout <<"\tNot deleting standard boundary data for dir:
          // "<<n<<"\n";
          if (!b) {
            // cout <<"\t    Zero boundary data pointer for dir:
            // "<<n<<"\n";
          }
          else if (!b->data.empty()) {
            // cout <<"\t    Emptying list.\n";
            b->data.clear();  // No dynamically allocated memory in
                              // these boundaries.
          }
          else {
            // cout <<"\t    Empty boundary list (!)(?) Probably
            // shouldn't happen.\n";
          }
        }  // if pre-existing.

        // cout <<"\t\tNow deleting boundary data...\n";
        if (b) {
          b = mem.myfree(b);
          // RT_recv_list[ii].RT_bd =
          // mem.myfree(RT_recv_list[ii].RT_bd);
        }
        // cout <<"\t\tFinished with this Receive boundary\n";
      }
      rc->RT_recv_list.clear();
    }

    //
    // Delete send lists.  No cells to delete, just the lists, and
    // boundary_data
    //
    if (!rc->RT_send_list.empty()) {
      // cout <<"\tDeleting RT Send List.\n";
      struct boundary_data *b;
      for (unsigned int ii = 0; ii < rc->RT_send_list.size(); ii++) {
        b = rc->RT_send_list[ii].RT_bd;
        if (!b) {
          // cout <<"\t\tZero send boundary data pointer for dir:
          // "<<RT_send_list[i].dir<<"\n";
        }
        else if (!b->data.empty()) {
          // cout <<"\t\tEmptying send list.\n";
          b->data.clear();  // No dynamically allocated memory in
                            // these boundaries.
        }
        else {
          // cout <<"\t\tEmpty send boundary list (!)(?) Probably
          // shouldn't happen.\n";
        }
        b = mem.myfree(b);
      }
      rc->RT_send_list.clear();
      // cout <<"\tDeleting RT Send List... done!\n";
    }
  }

  return;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::Setup_RT_Boundaries(
    class SimParams &par,       ///< pointer to simulation parameters
    class Sub_domain &ppar,     ///< domain decomposition info
    class GridBaseClass *grid,  ///< pointer to grid.
    const int src_id,
    struct rad_src_info &RS)
{
#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::Setup_RT_Boundaries() starting!");
#endif

  // ----------------------- SANITY CHECKS --------------------------
  //
  // Check that send/recv lists for this source have not already been set up!
  // Element i of the list should NOT correspond to source i, unless by
  // chance. Look through the ids in RT_source_list, just to make sure it is
  // not already there.
  //
  int sle = -1;  // source list element.
  for (int v = 0; v < static_cast<int>(RT_source_list.size()); v++) {
    if (RT_source_list[v].source_id == src_id) sle = v;
  }
  if (sle >= 0)
    spdlog::error(
        "{}: {}", "source is already present in RT_source_list!", sle);
  // ----------------------- SANITY CHECKS --------------------------

  // ------------------------ SETUP BOUNDARIES ----------------------
  //
  // Now we need to add a new send list and a new recv list for this source.
  //
  struct RT_source_comms_info this_src_comms;
  this_src_comms.source_id = src_id;
  this_src_comms.RT_recv_list.clear();
  this_src_comms.RT_send_list.clear();

  //
  // Now we call the appropriate function, depending on what type of
  // source it is.  The only two relevant types of source, as far as
  // boundary communication is concerned, are sources at infinity
  // and those at a finite position, because the boundary
  // communication differs between them.
  //
  int err = 0;
  if (!RS.at_infinity) {
    err += setup_RT_finite_ptsrc_BD(
        par, ppar, grid, this_src_comms.source_id, RS,
        this_src_comms.RT_recv_list, this_src_comms.RT_send_list);
  }
  else {
    err += setup_RT_infinite_src_BD(
        par, ppar, grid, this_src_comms.source_id, RS,
        this_src_comms.RT_recv_list, this_src_comms.RT_send_list);
  }

  //
  // Now we add this_src_comms to the list.
  //
  RT_source_list.push_back(this_src_comms);

  // ------------------------ SETUP BOUNDARIES ----------------------

#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::Setup_RT_Boundaries() finished!");
#endif
  return err;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::Receive_RT_Boundaries(
    class SimParams &par,       ///< pointer to simulation parameters
    class Sub_domain &ppar,     ///< domain decomposition info
    class GridBaseClass *grid,  ///< pointer to grid.
    const int src_id,           ///< source id
    struct rad_src_info &RS)
{
#ifdef RT_TESTING
  spdlog::info("\tReceive_RT_Boundaries() src={}: Starting\n", src_id);
#endif
  int err = 0;

  //
  // We need the recv-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle = -1;  // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
      spdlog::debug("\tRecv_RT_Boundaries() src={}", src_id << ": no procs to recv from.\n";
#endif
    return 0;
  }
  else {
    for (int v = 0; v < static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle < 0)
    spdlog::error("{}: {}", "source not found in RT_source_list", sle);
#ifdef RT_TESTING
  else
    spdlog::debug("\tRecv_RT_Boundaries() src={}: receiving data.\n", src_id);
#endif

  //  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
  //    spdlog::error("{}: {}", "Can't receive boundary for src_id b/c boundary
  //    was never set up",src_id);
  //  }
  //  if (RT_source_list[src_id].source_id != src_id) {
  //    spdlog::error("{}: {}", "source id for boundary indexed by src_id is
  //    wrong",src_id);
  //  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);

  //
  // See how many boundaries to receive, and create a list of them
  //
  int n_recv = i_src->RT_recv_list.size();
  if (n_recv == 0) {
#ifdef RT_TESTING
    spdlog::info(
        "RT_MPI_bc::Receive_RT_Boundaries() nothing to receive... returning.");
#endif
  }
  else {
    std::vector<int> r_dirs(n_recv);  // List of directions.
    std::vector<bool> r_done(
        n_recv);  // Will set these to true as I receive them.
#ifdef RT_TESTING
    spdlog::debug("\tReceive_RT_Boundaries() src={}: to recv from: ", src_id);
#endif
    stringstream s_array;
    for (int i = 0; i < n_recv; i++) {
      r_dirs[i] = i_src->RT_recv_list[i].dir;
      r_done[i] = false;
#ifdef RT_TESTING
      s_array << " " << r_dirs[i] << "[rank=" << i_src->RT_recv_list[i].rank
              << "]";
#endif
    }
#ifdef RT_TESTING
    spdlog::debug("{}", s_array.get());
#endif

    //
    // Now receive the n_recv boundaries, in any order:
    //
    double *buf = 0;
    struct boundary_data *b;
    int from_rank, comm_tag, ct, i_recv;
    string id;
    for (int i = 0; i < n_recv; i++) {
      id.erase();
      //
      // Find a message that has been sent to us.
      //
      from_rank = -1;
      err       = par.levels[0].sub_domain.look_for_data_to_receive(
          &from_rank, id, &comm_tag, BC_RTtag, COMM_DOUBLEDATA);
      if (err)
        spdlog::error("{}: {}", "Receive_RT_Boundaries() look4data", err);
      if (comm_tag != BC_RTtag) {
        spdlog::error("{}: {}", "Getting data, but not RT", comm_tag);
      }

      // Set pointer to boundary we are receiving, and note which one
      // it is in i_recv.
      b         = 0;
      bool test = false;
      i_recv    = -1;
      for (int v = 0; v < n_recv; v++) {
        if (i_src->RT_recv_list[v].rank == from_rank) {
          b      = i_src->RT_recv_list[v].RT_bd;
          i_recv = v;
          if (!test)
            test = true;
          else
            spdlog::error(
                "{}: {}", "2 boundaries from same rank!!!", from_rank);
        }
      }
      if (!b)
        spdlog::error(
            "{}: {}", "No Receive Boundary for this proc!", from_rank);
      if (r_dirs[i_recv] != i_src->RT_recv_list[i_recv].dir) {
        spdlog::error(
            "{}: {}", "direction mismatch",
            i_recv);  // this is paranoid checking!
      }
      if (r_done[i_recv])
        spdlog::error("{}: {}", "This boundary already received!!!", i_recv);

      //
      // See how many doubles we are expecting, and allocate memory.
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      ct = b->data.size() * 2 * RS.NTau;
      if (ct < 1) {
        spdlog::error("data size = {}, NTau = {}", b->data.size(), RS.NTau);
        spdlog::error("{}: {}", "Empty boundary!", ct);
      }
      vector<pion_flt> buf(ct);
#ifdef RT_TESTING
      spdlog::debug(
          "\tReceive_RT_Boundaries() rad_src={}: getting {} cells from proc. {}",
          src_id, ct, from_rank);
#endif

      //
      // Receive data into buffer.
      //
      err = par.levels[0].sub_domain.receive_double_data(
          from_rank,  ///< rank of process we are receiving from.
          comm_tag,   ///< comm_tag: what sort of comm we are looking for
                      ///< (PER,MPI,etc.)
          id,         ///< identifier for receive, for any book-keeping that
                      ///< might be needed.
          ct,         ///< number of doubles to receive
          buf         ///< Pointer to array to write to (must be already
                      ///< initialised).
      );
      if (err)
        spdlog::error("{}: {}", "Receive_RT_Boundaries() getdata failed", err);

      //
      // Loop over cells in boundary and copy received col values into
      // c->col. Hopefully they are ordered the same way in packing as
      // unpacking!!!
      //
      list<cell *>::iterator c = b->data.begin();
      int count                = 0;
      double tau[MAX_TAU];
      for (short unsigned int v = 0; v < MAX_TAU; v++)
        tau[v] = 0.0;

      for (c = b->data.begin(); c != b->data.end(); ++c) {
        if (count >= ct)
          spdlog::error("{}: {}", "too many cells!!!", count - ct);
#ifdef RT_TESTING
        if (count < 100) {
          spdlog::debug(
              "col = {} for cell {}: cell-id={}, pos=[{},{}", buf[count], count,
              (*c)->id, (*c)->pos[XX], (*c)->pos[YY]);
          if (par.ndim > 2)
            spdlog::debug(",{}]", (*c)->pos[ZZ]);
          else
            spdlog::debug("]");
        }
#endif
        // get column to and through cell.
        for (short unsigned int v = 0; v < RS.NTau; v++) {
          tau[v] = std::max(buf[count], 0.0);
          if (buf[count] < -SMALLVALUE) {
            spdlog::debug(
                "RECV neg tau cell : {}",
                std::vector<int>((*c)->pos, (*c)->pos + par.ndim));
            spdlog::debug("tau={} ... correcting to zero, ", buf[count]);
            CI.print_cell(*c);
          }
          count++;
        }
        CI.set_col(*c, src_id, tau);

        // get column through cell.
        for (short unsigned int v = 0; v < RS.NTau; v++) {
          tau[v] = std::max(buf[count], 0.0);
          if (buf[count] < -SMALLVALUE) {
            spdlog::debug(
                "RECV neg dtau cell : {}",
                std::vector<int>((*c)->pos, (*c)->pos + par.ndim));
            spdlog::debug("tau={} ... correcting to zero, ", buf[count]);
            CI.print_cell(*c);
          }
          count++;
        }
        CI.set_cell_col(*c, src_id, tau);

        (*c)->rt = false;
        // HACK
        /*
        cell *p=*c;
        if (grid->idx()==2 && (p->pos[XX]==127 || p->pos[XX]==129 ||
        p->pos[XX]==131 || p->pos[XX]==133)) { cout <<"RT MPI BC RECV:
        dir="<<b->dir<<": "; CI.print_cell(p);
        }
        */
        // HACK

      }  // loop over boundary cells.
      if (count != ct) spdlog::error("{}: {}", "BIG ERROR!", count - ct);

      //
      // Record boundary as received, so we know if we get the
      // same data twice.
      //
      r_done[i_recv] = true;
      if (err)
        spdlog::error(
            "{}: {}", "Error receiving boundary data from rank: ", from_rank);
    }  // Loop over boundaries.
  }    // else have boundaries.

#ifdef RT_TESTING
  spdlog::debug("\tReceive_RT_Boundaries() src={}: Finished\n", src_id);
#endif
  return err;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::Send_RT_Boundaries(
    class SimParams &par,       ///< pointer to simulation parameters
    class Sub_domain &ppar,     ///< domain decomposition info
    class GridBaseClass *grid,  ///< pointer to grid.
    const int src_id,           ///< source id
    struct rad_src_info &RS)
{
  int err = 0;
#ifdef RT_TESTING
  spdlog::debug("\tSend_RT_Boundaries() src={}: Starting\n", src_id);
#endif
  //
  // We need the send-boudary-list for src_id, which is found by looking
  // through the ids in RT_source_list
  //
  int sle = -1;  // source list element.
  if (RT_source_list.empty()) {
#ifdef RT_TESTING
    spdlog::debug(
        "\tSend_RT_Boundaries() src={}: no procs to send to.\n", src_id);
#endif
    return 0;
  }
  else {
    for (int v = 0; v < static_cast<int>(RT_source_list.size()); v++) {
      if (RT_source_list[v].source_id == src_id) sle = v;
    }
  }
  if (sle < 0)
    spdlog::error("{}: {}", "source not found in RT_source_list", sle);
#ifdef RT_TESTING
  else
    spdlog::debug("\tSend_RT_Boundaries() src={}: sending data.\n", src_id);
#endif

  //  if (RT_source_list.size()<=static_cast<unsigned int>(src_id)) {
  //    spdlog::error("{}: {}", "Can't send boundary for src_id b/c boundary was
  //    never set up",src_id);
  //  }
  //  if (RT_source_list[src_id].source_id != src_id) {
  //    spdlog::error("{}: {}", "source id for boundary indexed by src_id is
  //    wrong",src_id);
  //  }
  struct RT_source_comms_info *i_src = &(RT_source_list[sle]);

  //
  // See how many boundaries to send:
  // If none, just skip on to the par.levels[0].sub_domain.barrier() call at the
  // end. If we have some, do the sends.
  //
  int n_send = i_src->RT_send_list.size();
  if (n_send == 0) {
#ifdef RT_TESTING
    spdlog::info(
        "RT_MPI_bc::Send_RT_Boundaries() nothing to send...  returning.");
#endif
  }
  else {
    //
    // Allocate memory for send buffers.
    //
#ifdef RT_TESTING
    spdlog::debug(
        "\tSend_RT_Boundaries() src={}: Allocating {} data buffers.", src_id,
        n_send);
#endif
    std::string *id = 0;
    id              = mem.myalloc(id, n_send);
    //
    // Loop over Boundaries, put c->col into data, send it.
    // I'm going to use non-blocking sends, so they all send at once.
    //
    for (int i = 0; i < n_send; i++) {
      struct boundary_data *b = i_src->RT_send_list[i].RT_bd;

#ifdef RT_TESTING
      spdlog::debug(
          "\tSend_RT_Boundaries() src={}: Sending boundary in dir {}", src_id,
          i_src->RT_send_list[i].dir);
#endif
      //
      // How many cell column densities to send:
      // For sources with more than one optical depth, we need to
      // account for this.
      //
      const int nc = b->data.size() * 2 * RS.NTau;
      vector<double> data(nc);

      //
      // Put columns into data
      //
      list<cell *>::iterator c = b->data.begin();
      int count                = 0;
      double tau[MAX_TAU];
      for (short unsigned int v = 0; v < MAX_TAU; v++)
        tau[v] = 0.0;

      for (c = b->data.begin(); c != b->data.end(); ++c) {
        if (count >= nc)
          spdlog::error("{}: {}", "too many cells!!!", count - nc);

          // HACK
          /*
          cell *p=*c;
          if (grid->idx()==2 && (p->pos[XX]==127 || p->pos[XX]==129))
          { cout <<"RT MPI BC SEND: dir="<<b->dir<<": ";
            CI.print_cell(p);
          }
          */
          // HACK

#ifdef RT_TESTING
        if (count < 100) {
          spdlog::debug("send data [{}]: col[0] = ", i);
          CI.get_col(*c, src_id, tau);
          spdlog::debug(
              "{} for cell {}: pos=[{},{}", tau[0], count, (*c)->pos[XX],
              (*c)->pos[YY]);
          if (par.ndim > 2)
            spdlog::debug(",{}]\n", (*c)->pos[ZZ]);
          else
            spdlog::debug("]");
        }
#endif
        // column to and through cell.
        CI.get_col(*c, src_id, tau);
        for (short unsigned int v = 0; v < RS.NTau; v++) {
          data[count] = std::max(tau[v], 0.0);
          if (tau[v] < -SMALLVALUE) {
            spdlog::debug(
                "SEND neg tau cell : {}",
                std::vector<int>((*c)->pos, (*c)->pos + par.ndim));
            spdlog::debug("tau={} ", tau[v]);
            CI.print_cell(*c);
            // tau[v]=0.0;
          }
          count++;
        }

        // column through cell.
        CI.get_cell_col(*c, src_id, tau);
        for (short unsigned int v = 0; v < RS.NTau; v++) {
          data[count] = std::max(tau[v], 0.0);
          if (tau[v] < -SMALLVALUE) {
            // Occurs if ionization fraction is slightly >1,
            // because of discretisation or roundoff errors.
#ifdef RT_TESTING
            spdlog::debug(
                "SEND neg dtau cell : {}",
                std::vector<int>((*c)->pos, (*c)->pos + par.ndim));
            spdlog::debug("tau={} ", tau[v]);
            CI.print_cell(*c);
#endif
            // tau[v]=0.0;
          }
          count++;
        }
      }

      //
      // Send data
      //
#ifdef RT_TESTING
      spdlog::debug(
          "Send_BC[{}]: Sending {} doubles to proc: {} in direction {}", i, nc,
          i_src->RT_send_list[i].rank, i_src->RT_send_list[i].dir);
#endif
      err += par.levels[0].sub_domain.send_double_data(
          i_src->RT_send_list[i].rank,  ///< rank to send to.
          nc,       ///< size of buffer, in number of doubles.
          data,     ///< vector of data
          id[i],    ///< identifier for send, for tracking delivery later.
          BC_RTtag  ///< comm_tag, to say what kind of send this is.
      );
      if (err) spdlog::error("{}: {}", "Send_BC[] send_data failed.", err);
#ifdef RT_TESTING
      spdlog::debug("Send_BC[{}]: send returned with id={}", i, id[i]);
#endif
    }  // loop over all sends to send data.

    //
    // Loop over all sends, and wait until they complete!
    //
#ifdef RT_TESTING
    spdlog::debug(
        "\tSend_RT_Boundaries() src={}: Waiting for sends to complete.\n",
        src_id);
#endif
    for (int i = 0; i < n_send; i++) {
      err += par.levels[0].sub_domain.wait_for_send_to_finish(id[i]);
    }

    //
    // delete array of ids
    //
    id = mem.myfree(id);
  }  // else (have boundaries to send).

  //
  // Wait for all processors to finish RT for this source, so no messages get
  // mixed up! Not sure if this is strictly necessary... so should try to do
  // without it, once I know the algorithm is sound.
  //
#ifdef RT_TESTING
  spdlog::debug(
      "\tSend_RT_Boundaries() src={}: Waiting for all procs to finish.\n",
      src_id);
#endif
  par.levels[0].sub_domain.barrier();
  return err;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::setup_RT_infinite_src_BD(
    class SimParams &par,       ///< pointer to simulation parameters
    class Sub_domain &ppar,     ///< domain decomposition info
    class GridBaseClass *grid,  ///< pointer to grid.
    const int src_id,
    struct rad_src_info &RS,
    std::vector<struct RT_boundary_list_element> &RT_recv_list,
    std::vector<struct RT_boundary_list_element> &RT_send_list)
{
#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::setup_RT_diffuse_radn_src_BD() starting.");
#endif
  int err = 0;
  //
  // Make sure src is at infinity, and has the correct type.
  //
  if (!RS.at_infinity) {
    spdlog::error(
        "{}: {}", "Source for diffuse radiation is not at infinity", src_id);
  }
  //
  // If the previous function returned at all, then the source
  // exists, so no need to check that.  For a source at infinity
  // there can be only one send and one receive domain (or none, if
  // current domain an edge domain).
  //

  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
    spdlog::error(
        "{}: {}", "diffuse radn src recv-list not empty!", RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
    spdlog::error(
        "{}: {}", "diffuse radn src send-list not empty!", RT_send_list.size());
  }

  //
  // Get source direction:
  //
  enum direction srcdir = RT_src_at_infty_direction(par, src_id, RS);
  if (srcdir == NO) {
    spdlog::error("{}: {}", "src position is not at infinity", srcdir);
  }
  enum direction recv_dir = srcdir;
  enum direction send_dir = grid->OppDir(srcdir);
  //
  // Recv boundary in the source direction, send boundary in opp.dir.
  // Can have only one send and one recv source, so the lists should
  // get exactly one element.
  //
  int recv_proc = ppar.get_neighbour_rank(recv_dir);
  int send_proc = ppar.get_neighbour_rank(send_dir);

  //
  // Set up the lone recv boundary list element.  Note we leave an
  // empty vector if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempR;
  tempR.RT_bd = 0;
  switch (recv_dir) {
    case XN:
      tempR.dir = static_cast<int>(dir_XN);
      break;
    case XP:
      tempR.dir = static_cast<int>(dir_XP);
      break;
    case YN:
      tempR.dir = static_cast<int>(dir_YN);
      break;
    case YP:
      tempR.dir = static_cast<int>(dir_YP);
      break;
    case ZN:
      tempR.dir = static_cast<int>(dir_ZN);
      break;
    case ZP:
      tempR.dir = static_cast<int>(dir_ZP);
      break;
    case NO:
      spdlog::debug(
          "\t No processor in receive direction d={}: proc={}", recv_dir,
          recv_proc);
      tempR.dir = -1;
      break;
    default:
      spdlog::debug(
          "\t No processor in receive direction d={}: proc={}", recv_dir,
          recv_proc);
      tempR.dir = -1;
      break;
  }
  tempR.rank = recv_proc;
  if (tempR.rank >= 0) {
#ifdef RT_TESTING
    spdlog::debug(
        "\t\tFound diffuse-recv-proc in dir={} (R.dir={}), rank={}, setting up boundary data.\n",
        recv_dir, tempR.dir, tempR.rank);
#endif
    RT_recv_list.push_back(tempR);
    err = setup_RT_recv_boundary(grid, RT_recv_list.front());
    if (err)
      spdlog::error(
          "{}: {}", "failed to set up diffuse radn recv boundary", err);
  }
  else {
#ifdef RT_TESTING
    spdlog::debug(
        "\t\tDiffuse-recv-proc doesn't exist in dir={}; I must be at edge of domain.\n",
        recv_dir);
#endif
  }

  //
  // Set up the lone send boundary list element.  Note we leave an
  // empty vector if there is no neighbour domain to add.
  //
  struct RT_boundary_list_element tempS;
  tempS.RT_bd = 0;
  switch (send_dir) {
    case XN:
      tempS.dir = static_cast<int>(dir_XN);
      break;
    case XP:
      tempS.dir = static_cast<int>(dir_XP);
      break;
    case YN:
      tempS.dir = static_cast<int>(dir_YN);
      break;
    case YP:
      tempS.dir = static_cast<int>(dir_YP);
      break;
    case ZN:
      tempS.dir = static_cast<int>(dir_ZN);
      break;
    case ZP:
      tempS.dir = static_cast<int>(dir_ZP);
      break;
    case NO:
#ifdef RT_TESTING
      spdlog::debug(
          "\t No processor in send direction d={}: proc={}", send_dir,
          send_proc);
#endif
      tempS.dir = -1;
      break;
    default:
#ifdef RT_TESTING
      spdlog::debug(
          "\t No processor in send direction d={}: proc={}", send_dir,
          send_proc);
#endif
      tempS.dir = -1;
      break;
  }
  tempS.rank = send_proc;
  if (tempS.rank >= 0) {
#ifdef RT_TESTING
    spdlog::debug(
        "\t\tFound diffuse-send-proc in dir={}, rank={}, setting up boundary data.",
        send_dir, tempS.rank);
#endif
    RT_send_list.push_back(tempS);
    err = setup_RT_send_boundary(grid, RT_send_list.front());
    if (err)
      spdlog::error(
          "{}: {}", "failed to set up diffuse radn send boundary", err);
  }
  else {
#ifdef RT_TESTING
    spdlog::debug(
        "\t\tDiffuse-send-proc doesn't exist in dir={}; I must be at edge of domain.",
        send_dir);
#endif
  }

#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::setup_RT_diffuse_radn_src_BD() returning.");
#endif
  return err;
}



// ##################################################################
// ##################################################################



enum direction RT_MPI_bc::RT_src_at_infty_direction(
    class SimParams &par,  ///< pointer to simulation parameters
    const int src_id,
    struct rad_src_info &RS)
{
  //
  // Get source position vector; compare values to find the
  // very large value which identifies the direction of the
  // infinite source.
  //
  enum direction srcdir = NO;
  for (int v = 0; v < par.ndim; v++) {
    if (RS.pos[v] < -1.e99) srcdir = static_cast<direction>(2 * v);
    if (RS.pos[v] > 1.e99) srcdir = static_cast<direction>(2 * v + 1);
  }
  if (srcdir == NO) {
    spdlog::warn(
        "WARNING: RT source is not at infinity! RT_src_at_infty_direction()");
    spdlog::debug("Srcpos : {}", RS.pos);
  }
  return srcdir;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::setup_RT_finite_ptsrc_BD(
    class SimParams &par,       ///< pointer to simulation parameters
    class Sub_domain &ppar,     ///< domain decomposition info
    class GridBaseClass *grid,  ///< pointer to grid.
    const int src_id,
    struct rad_src_info &RS,
    std::vector<struct RT_boundary_list_element> &RT_recv_list,
    std::vector<struct RT_boundary_list_element> &RT_send_list)
{
#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::setup_RT_finite_ptsrc_BD() starting.");
#endif
  //
  // Check that send and recv lists are empty!
  //
  if (!RT_recv_list.empty()) {
    spdlog::error(
        "{}: {}", "Monochromatic point src recv-list not empty!",
        RT_recv_list.size());
  }
  if (!RT_send_list.empty()) {
    spdlog::error(
        "{}: {}", "Monochromatic point src send-list not empty!",
        RT_send_list.size());
  }

  //
  // First we find the source:
  //
  // Source should be at a cell corner, so we need to be careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processors!), so
  // this should work fine.
  //
  std::vector<enum direction> srcdir(par.ndim);
  std::array<double, MAX_DIM> srcpos;
  for (int v = 0; v < par.ndim; v++) {
    srcpos[v] = RS.pos[v];
  }

  std::array<int, MAX_DIM> i_srcpos;
  CI.get_ipos_vec(srcpos, i_srcpos);

  for (int i = 0; i < par.ndim; i++) {
    if (i_srcpos[i] < grid->iXmin(static_cast<axes>(i))) {
#ifdef RT_TESTING
      spdlog::debug(
          "*** axis={} srcpos={} xmin={} src is off grid in neg.dir.", i,
          srcpos[i], grid->Xmin(static_cast<axes>(i)));
#endif
      srcdir[i] = static_cast<direction>(2 * i);
    }
    else if (i_srcpos[i] > grid->iXmax(static_cast<axes>(i))) {
#ifdef RT_TESTING
      spdlog::debug(
          "*** axis={} srcpos={} xmax={} src is off grid in pos.dir.\n", i,
          srcpos[i], grid->Xmax(static_cast<axes>(i)));
#endif
      srcdir[i] = static_cast<direction>(2 * i + 1);
    }
    else {
#ifdef RT_TESTING
      spdlog::debug(
          "*** axis={} srcpos={} xmin={} xmax={} src is on local domain.", i,
          srcpos[i], grid->Xmin(static_cast<axes>(i)),
          grid->Xmax(static_cast<axes>(i)));
#endif
      srcdir[i] = static_cast<direction>(-1);
    }
  }

  //
  // Choose processors I need to receive data from and send data to.
  //
  // First see if direction of source off grid has a neigbour
  // processor andif the opposite direction has a neighbour
  // processor.
  //
  // recv procs:
  bool recv_proc_exists[MAX_DIM];
  for (int i = 0; i < MAX_DIM; i++)
    recv_proc_exists[i] = false;
  for (int i = 0; i < par.ndim; i++) {
    enum direction posdir = static_cast<direction>(2 * i + 1);
    enum direction negdir = static_cast<direction>(2 * i);
    if ((srcdir[i] == negdir)
        && (!pconst.equalD(
               grid->Xmin(static_cast<axes>(i)),
               grid->level_Xmin(static_cast<axes>(i)))))
      recv_proc_exists[i] = true;
    else if (
        (srcdir[i] == posdir)
        && (!pconst.equalD(
               grid->Xmax(static_cast<axes>(i)),
               grid->level_Xmax(static_cast<axes>(i)))))
      recv_proc_exists[i] = true;
    else
      recv_proc_exists[i] = false;
  }

  //
  // send procs:
  //
  // Source should be at a cell corner, so we need to be more careful
  // with the less than/greater than questions.  We will define the
  // source to be on the domain if it is within the domain or on the
  // domain boundary.  The raytracer moves the source to the nearest
  // cell vertex (which should be consistent across processes), so
  // this should work fine.
  //
  // So we only send data to a neighbour if no part of it is in a
  // source plane.
  //
  bool send_proc_exists[2 * MAX_DIM];
  for (int i = 0; i < 2 * MAX_DIM; i++)
    send_proc_exists[i] = false;
  for (int i = 0; i < par.ndim; i++) {
    enum direction posdir = static_cast<direction>(2 * i + 1);
    enum direction negdir = static_cast<direction>(2 * i);
    if ((i_srcpos[i] > grid->iXmin(static_cast<axes>(i)))
        && (!pconst.equalD(
               grid->Xmin(static_cast<axes>(i)),
               grid->level_Xmin(static_cast<axes>(i)))))
      send_proc_exists[negdir] = true;
    else
      send_proc_exists[negdir] = false;
    // either doesn't exist, or we don't need it.

    if ((i_srcpos[i] < grid->iXmax(static_cast<axes>(i)))
        && (!pconst.equalD(
               grid->Xmax(static_cast<axes>(i)),
               grid->level_Xmax(static_cast<axes>(i)))))
      send_proc_exists[posdir] = true;
    else
      send_proc_exists[posdir] = false;
    // either doesn't exist, or we don't need it.
  }
#ifdef RT_TESTING
  spdlog::debug("send_proc_exists : {}", send_proc_exists);
#endif

  //
  // Get ranks and directions of MPI processes to receive from.
  //
  struct RT_boundary_list_element temp;
  temp.RT_bd = 0;
  for (int d = 0; d < par.ndim; d++) {
    if (recv_proc_exists[d] == true) {
      temp.rank = ppar.get_neighbour_rank(srcdir[d]);
      temp.dir  = srcdir[d];
      RT_recv_list.push_back(temp);
    }
  }

  //
  // Choose processes I need to send data to.  Could be all directions if
  // the source is within my subdomain
  //
  temp.RT_bd = 0;
  for (int dd = 0; dd < 2 * par.ndim; dd++) {
    if (send_proc_exists[dd]) {
      temp.rank = ppar.get_neighbour_rank(dd);
      temp.dir  = static_cast<enum direction>(dd);
      RT_send_list.push_back(temp);
    }
  }

#ifdef RT_TESTING
  spdlog::debug(
      "send list: {} boundaries.\nrecv list: {} boundaries.\n",
      RT_send_list.size(), RT_recv_list.size());
  stringstream s_array;
  s_array << "SEND LIST DIRECTIONS: ";
  for (unsigned int i = 0; i < RT_send_list.size(); i++)
    s_array << RT_send_list[i].dir << " ";
  s_array << "\n";
  s_array << "RECV LIST DIRECTIONS: ";
  for (unsigned int i = 0; i < RT_recv_list.size(); i++)
    s_array << RT_recv_list[i].dir << " ";
  spdlog::debug("{}", s_array);
#endif

  //
  // Set up receive boundaries.  First initialise all boundary_data pointers
  // to zero, then setup faces.
  //
  int err = 0;
  for (unsigned int i = 0; i < RT_recv_list.size(); i++)
    RT_recv_list[i].RT_bd = 0;
  //
  // First set up face boundaries.
  //
  for (unsigned int i = 0; i < RT_recv_list.size(); i++) {
    int d = RT_recv_list[i].dir;
    if (d == dir_XN || d == dir_XP || d == dir_YN || d == dir_YP || d == dir_ZN
        || d == dir_ZP) {
      err += setup_RT_recv_boundary(grid, RT_recv_list[i]);
    }
  }
  //
  // Make sure we initialised all the boundaries correctly!
  //
  for (unsigned int i = 0; i < RT_recv_list.size(); i++) {
    if (!(RT_recv_list[i].RT_bd)) {
      spdlog::error(
          "{}: {}", "Missed out on receive boundary!", RT_recv_list[i].dir);
    }
  }
  if (err) spdlog::error("{}: {}", "failed to setup RT recv boundaries", err);

  for (unsigned int i = 0; i < RT_send_list.size(); i++) {
    RT_send_list[i].RT_bd = 0;
    err += setup_RT_send_boundary(grid, RT_send_list[i]);
  }
  for (unsigned int i = 0; i < RT_send_list.size(); i++) {
    if (!(RT_send_list[i].RT_bd)) {
      spdlog::error("{}: {}", "Missed out on send boundary!", i);
    }
  }
  if (err) spdlog::error("{}: {}", "failed to setup RT send boundaries", err);

#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::setup_RT_monochromatic_ptsrc_BD() finished!");
#endif
  return 0;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::setup_RT_recv_boundary(
    class GridBaseClass *grid,          ///< pointer to grid.
    struct RT_boundary_list_element &b  ///< boundary info.
)
{
#ifdef RT_TESTING
  spdlog::info("RT_MPI_bc::setup_RT_recv_boundary() starting (dir={}).", b.dir);
#endif
  int err = 0;

  //
  // First set up boundary data, since it starts as an uninitialised pointer.
  //
  if (b.RT_bd)
    spdlog::error(
        "{}: {}", "Already set up RT boudary data!", fmt::ptr(b.RT_bd));
  b.RT_bd = mem.myalloc(b.RT_bd, 1);
  enum direction d1;

  if (grid->BC_bd.size() == 0)
    spdlog::error(
        "{}: {}", "setup regular boundaries before RT ones!!!", b.dir);

  switch (b.dir) {
      //
      // Each face at a time, already have boundaries set up, so need
      // to connect the boundary cells to each other laterally.
      //
    case dir_XN:
      d1 = XN;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;
    case dir_XP:
      d1 = XP;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;
    case dir_YN:
      d1 = YN;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;
    case dir_YP:
      d1 = YP;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;
    case dir_ZN:
      d1 = ZN;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;
    case dir_ZP:
      d1 = ZP;
      err += RT_populate_recv_boundary(b.RT_bd, grid->BC_bd[d1], d1);
      break;

    default:
      spdlog::error("{}: {}", "bad dir!", b.dir);
  }

#ifdef RT_TESTING
  spdlog::info(
      "RT_MPI_bc::setup_RT_recv_boundary() returning (dir={}).", b.dir);
#endif
  return err;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::RT_populate_recv_boundary(
    struct boundary_data *b,         ///< pointer to RT boundary data.
    const struct boundary_data *b2,  ///< pointer to BC boundary data.
    const enum direction offdir      ///< face dir
)
{
  // if (par.ndim==1) return 0; // no connections to make!
  if (!b || !b2)
    spdlog::error(
        "{}: {}", "RT_MPI_bc::RT_populate_recv_boundary() Null b/b2 pointer",
        fmt::ptr(b));

#ifdef RT_TESTING
  spdlog::debug(
      "RT_populate_recv_boundary() offdir={},  ondir={}", offdir, b2->ondir);
#endif  // RT_TESTING

  //
  // Add all boundary cells to the RT boundary list.
  //
  list<cell *>::const_iterator bpt = b2->data.begin();
  do {
    b->data.push_back(*bpt);
#ifdef RT_TESTING
    spdlog::debug("RT_populate_recv_boundary() cpos= ");
    spdlog::debug(" : {}", (*bpt)->pos);
#endif  // RT_TESTING
    ++bpt;
  } while (bpt != b2->data.end());
  return 0;
}



// ##################################################################
// ##################################################################



int RT_MPI_bc::setup_RT_send_boundary(
    class GridBaseClass *grid,               ///< pointer to grid.
    struct RT_boundary_list_element &send_b  ///< boundary info.
)
{
#ifdef RT_TESTING
  spdlog::debug(
      "RT_MPI_bc::setup_RT_send_boundary() starting (dir={}).\n", send_b.dir);
#endif
  int err = 0;
  //
  // get a pointer to the existing grid boundary.
  //
  boundary_data *grid_b = grid->BC_bd[send_b.dir];
  if (!grid_b) {
    spdlog::error(
        "{}: {}", "RT boundary in dir with no real boundary!", send_b.dir);
  }

  //
  // Set up boundary data.
  //
  if (send_b.RT_bd) {
    spdlog::error(
        "{}: {}", "Already set up RT boudary data!", fmt::ptr(send_b.RT_bd));
  }
  send_b.RT_bd = mem.myalloc(send_b.RT_bd, 1);

  //
  // Now get every cell in the grid boundary, and add the 'Nbc'th
  // cell in the on-grid-direction to the send-boundary list.
  //
  list<cell *>::const_iterator bpt = grid_b->data.begin();
  cell *target                     = 0;
  do {
    target = (*bpt);
    for (int v = 0; v < grid_b->depth; v++)
      target = grid->NextPt(target, grid_b->ondir);

    send_b.RT_bd->data.push_back(target);
#ifdef RT_TESTING
    spdlog::debug("setup_RT_send_boundary() cpos= ");
    spdlog::debug("  : {}", target->pos);
#endif  // RT_TESTING
    ++bpt;
  } while (bpt != grid_b->data.end());

#ifdef RT_TESTING
  spdlog::debug(
      "RT_MPI_bc::setup_RT_send_boundary() returning (dir={}).", send_b.dir);
#endif
  return err;
}



// ##################################################################
// ##################################################################
