/*
 *@BEGIN LICENSE
 *
 * uhfno by Psi4 Developer, a plugin to:
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>

#include "newuhf.h"

INIT_PLUGIN

using namespace boost;

namespace psi{ namespace scf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "UHFNO"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
    }

    return true;
}

extern "C" 
PsiReturnType uhfno(Options& options)
{
    int print = options.get_int("PRINT");

    boost::shared_ptr<Wavefunction> oldwfn = Process::environment.wavefunction();
    boost::shared_ptr<PSIO> psio(new PSIO);
    boost::shared_ptr<Wavefunction> scf(new NewUHF(options, psio));    
    scf->compute_energy();
    Process::environment.set_wavefunction(oldwfn); // avoid segfault on exit

    return Success;
}

}} // End namespaces

