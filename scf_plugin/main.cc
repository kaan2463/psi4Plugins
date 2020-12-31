/*
 *@BEGIN LICENSE
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

#include "psi4/libplugin/plugin.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/psi4-dec.h"
#include "psi4/libmints/wavefunction.h"
#include "scf_plugin.h"

extern "C" PSI_API
int read_options(std::string name, Options &options)
{
    if (name == "SCF_PLUGIN"|| options.read_globals()) {
        options.add_int("SCF_E_CONVERGENCE", 8);
        options.add_int("SCF_D_CONVERGENCE", 6);
        options.add_int("SCF_MAXITER", 50);
        options.add_double("E_CONVERGENCE", 1.0E-10);
        options.add_double("D_CONVERGENCE", 1.0E-6);
        options.add_str("METHOD","RHF");
    }

    return true;
}

extern "C" PSI_API
SharedWavefunction scf_plugin(SharedWavefunction ref_wfn, Options& options)
{
    SharedWavefunction scf_plugin = SharedWavefunction(new SCF_PLUGIN(ref_wfn, options));
    scf_plugin->compute_energy();
    return ref_wfn;
}
