#include "libmints/mints.h"
#include "libscf_solver/uhf.h"
#include "newuhf.h"

namespace psi { namespace scf {

NewUHF::NewUHF(Options &options, boost::shared_ptr<PSIO> psio): UHF(options, psio) { }
NewUHF::~NewUHF() {}

void NewUHF::compute_nos()
{
  // Build S^1/2
  SharedMatrix eigvec = factory_->create_shared_matrix("L");
  SharedVector eigval(factory_->create_vector());

  S_->diagonalize(eigvec, eigval);
  int *dimpi = eigval->dimpi();
  for(int h=0; h < nirrep_; h++)
    for(int i=0; i < dimpi[h]; i++) {
      double shalf = sqrt(eigval->get(h,i));
      eigval->set(h, i, shalf);
    }

  SharedMatrix shalf = factory_->create_shared_matrix("sqrt S eigenvalues");
  shalf->set_diagonal(eigval);

  SharedMatrix TMP = factory_->create_shared_matrix("shalf * L+");
  TMP->gemm(false, true, 1.0, shalf, eigvec, 0.0);
  SharedMatrix SHalf = factory_->create_shared_matrix("L * shalf * L+");
  SHalf->gemm(false, false, 1.0, eigvec, TMP, 0.0);

  // Diagonalize S^1/2 Dt S^1/2
  TMP->gemm(false, false, 1.0, SHalf, Dt_, 0.0);
  SharedMatrix SDS = factory_->create_shared_matrix("S^1/2 Dt S^1/2");
  SDS->gemm(false, false, 1.0, TMP, SHalf, 0.0);

  SharedMatrix UHF_NOs = factory_->create_shared_matrix("UHF NOs");
  SharedVector UHF_NOONs(factory_->create_vector());
  SDS->diagonalize(UHF_NOs, UHF_NOONs);
  UHF_NOONs->print();
}

void NewUHF::finalize()
{
  compute_nos();
  UHF::finalize();
}

}} // psi::scf
