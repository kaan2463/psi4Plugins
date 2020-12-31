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

#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/basisset.h"
#include "scf_plugin.h"
#include "psi4/psi4-dec.h"
#include <cmath>
#include "psi4/libpsi4util/PsiOutStream.h"
#include <memory>

SCF_PLUGIN::SCF_PLUGIN(SharedWavefunction ref_wfn, Options &options) : Wavefunction(options)
{
	reference_wavefunction_ = ref_wfn;
	shallow_copy(ref_wfn);
	title();
	common_init();
}

SCF_PLUGIN::~SCF_PLUGIN()
{

}

void SCF_PLUGIN::common_init()
{

    scf_maxiter=options_.get_int("SCF_MAXITER");
    tol_E= options_.get_double("E_CONVERGENCE");
    tol_D=options_.get_double("D_CONVERGENCE");
    method=options_.get_str("METHOD");

    MintsHelper mints(reference_wavefunction_->basisset());
    Enuc = reference_wavefunction_->molecule()->nuclear_repulsion_energy(reference_wavefunction_->get_dipole_field_strength());


    nso = reference_wavefunction_->nso();
    noccA = reference_wavefunction_->doccpi()[0]+reference_wavefunction_->soccpi()[0];
    noccB = reference_wavefunction_->doccpi()[0];
    //S: overlap matrix
    Sso = std::make_shared<Tensor2d>("Overlap Matrix",nso,nso);
    Sso->set(mints.ao_overlap());
    //S^-1/2: 

    SharedTensor2d eigen_vectors = std::make_shared<Tensor2d>(nso,nso);
    SharedTensor1d eigen_values = std::make_shared<Tensor1d>(nso);
    Sso->diagonalize(eigen_vectors,eigen_values,CUT_OFF_DIAGONALIZE);

    SharedTensor2d powered_eigen_matrix = std::make_shared<Tensor2d>(nso,nso);

    for(int i=0;i<nso;i++){
        powered_eigen_matrix->set(i,i,1.0/sqrt(eigen_values->get(i)));
    }

    X=std::make_shared<Tensor2d>("S^-1/2",nso,nso);

    X ->back_transform(powered_eigen_matrix,eigen_vectors);


    //X = UHF_PLUGIN::power(Sso,-0.5);
    //T: kinetic energy integrals
    Tso = std::make_shared<Tensor2d>("kinetic energy integrals matrix",nso,nso);
    Tso->set(mints.ao_kinetic());
    //V : potential energy integrals
    Vso =std::make_shared<Tensor2d>("potential energy integrals matrix",nso,nso);
    Vso->set(mints.ao_potential());
    //H : one electron integral matrix
    Hso=std::make_shared<Tensor2d>("one e integral matrix",nso,nso);
    Hso->set(Tso);
    Hso->add(Vso);
    //2 electron integrals
    Eri=std::make_shared<Tensor2d>("2 electron integral matrix",nso*nso,nso*nso);
    Eri->set(mints.ao_eri());
    
}


double SCF_PLUGIN::compute_energy()
{
    printf("Method = %s\n",method.c_str());
    if(method.compare("UHF")==0){
         printf("uhf is running\n");
        return uhf();
    } else if(method.compare("RHF")==0){
         printf("rhf is running\n");
        return rhf();
    }
     printf("is unrecognized!!\n");
    return 0.0; 
}

void SCF_PLUGIN::title()
{
	outfile->Printf("SCF PLUGIN\n");
}
