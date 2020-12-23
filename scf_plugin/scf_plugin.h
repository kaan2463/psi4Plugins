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

#ifndef scf_plugin_h
#define scf_plugin_h

#include "psi4/libmints/wavefunction.h"
#include <string>
#include "psi4/libfock/jk.h"
#include "tensors.h"


using namespace psi;
using namespace scf_plugin; 

class SCF_PLUGIN : public Wavefunction
{
    public:
        SCF_PLUGIN(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
        virtual ~SCF_PLUGIN();

        virtual double compute_energy();

    protected:
        // Functions
        void title();
        void common_init();
        
        // Variables
        SharedTensor2d Sso;
        SharedTensor2d Tso;
        SharedTensor2d Vso;
        SharedTensor2d Hso;
        SharedTensor2d Eri;

        double tol_E; 
        double tol_D; 
        int scf_maxiter;  
        int nso;
        int doccA;
        double Enuc;
};

#endif
