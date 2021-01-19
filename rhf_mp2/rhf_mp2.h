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

#ifndef mp_plugin_h
#define mp_plugin_h

#include "psi4/libmints/wavefunction.h"
#include <string>
#include "psi4/libfock/jk.h"
#include "tensors.h"

#define ABS(x) (x<0?-x:x)
#define CUT_OFF_DIAGONALIZE 1.0E10-15

using namespace psi;
using namespace psi4_plugin; 

class MP_PLUGIN : public Wavefunction
{
    public:
        MP_PLUGIN(std::shared_ptr<Wavefunction> reference_wavefunction, Options &options);
        virtual ~MP_PLUGIN();

        virtual double compute_energy();

    protected:
        // Functions
        void title();
        void common_init();

        // Variables
        double tol_E; 
        double tol_D; 
        int scf_maxiter; 
        int nso;

        SharedTensor2d fockA;
        SharedTensor2d bQiaA;

        int nfrzv;
        int nfrzc; 
        int nmo;
        int noccA;
        int nQ ;
        int noccB;
        int naoccA;
        int nvirA;
        int navirA;
        int nvirB;
        double Enuc;
        double Escf;

        std::string method;
};

#endif
