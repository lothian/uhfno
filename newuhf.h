#ifndef _newUHF_h_
#define _newUHF_h_

#include "libscf_solver/uhf.h"
#include "libmints/mints.h"
#include "boost/shared_ptr.hpp"

namespace psi { namespace scf {

class NewUHF : public UHF {
public:
  NewUHF(Options &options, boost::shared_ptr<PSIO> psio);
  virtual ~NewUHF();

protected:
  virtual void compute_nos();
  virtual void finalize();
  virtual void common_init();
};

}} // psi::scf

#endif
