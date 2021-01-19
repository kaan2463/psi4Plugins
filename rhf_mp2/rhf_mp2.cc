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
#include "psi4/libmints/vector.h"
#include "psi4/lib3index/dftensor.h"
#include "rhf_mp2.h"
#include "psi4/psi4-dec.h"
#include "tensors.h"
#include <cmath>
#include "psi4/libpsi4util/PsiOutStream.h"
#include <memory>

MP_PLUGIN::MP_PLUGIN(SharedWavefunction ref_wfn, Options &options) : Wavefunction(options)
{
	reference_wavefunction_ = ref_wfn;
	shallow_copy(ref_wfn);
	title();
	common_init();
}

MP_PLUGIN::~MP_PLUGIN()
{

}

void MP_PLUGIN::common_init()
{
    // Auxiliary basis set
    std::shared_ptr<BasisSet> auxiliary = reference_wavefunction_->get_basisset("DF_BASIS_SCF");
    nQ = auxiliary->nbf();
    noccA = reference_wavefunction_->doccpi()[0];
    nfrzv = reference_wavefunction_->frzvpi()[0];
    nfrzc = reference_wavefunction_->frzcpi()[0];
    naoccA=noccA-nfrzc; 
    nmo = reference_wavefunction_->nmo();
    nvirA = nmo - noccA;
    navirA=nvirA-nfrzv;

    fockA = SharedTensor2d(new Tensor2d("Fock Matrix", nmo, nmo));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC A (Q|IA)", nQ, naoccA, navirA));
    SharedVector epsilon_a = reference_wavefunction_->epsilon_a();
    for(int i=0;i<nmo;i++){
        fockA->set(i,i,epsilon_a->get(i));
    }

    // DFTensor A
    std::shared_ptr<DFTensor> DFA (new DFTensor(reference_wavefunction_->basisset(), auxiliary, reference_wavefunction_->Ca(), noccA, nvirA, naoccA, navirA, options_));

    // Alpha OV Block
    SharedMatrix ovA = DFA->Qov();
    bQiaA->set(ovA);
    ovA.reset();
}


double MP_PLUGIN::compute_energy()
{
    SharedTensor2d T=SharedTensor2d(new Tensor2d("Tiajb", naoccA, navirA,naoccA,navirA));
    SharedTensor2d TD=SharedTensor2d(new Tensor2d("<ij|ab>", naoccA, navirA, naoccA,navirA));
    SharedTensor2d U=SharedTensor2d(new Tensor2d("Uijab", naoccA, navirA, naoccA,navirA));

    T->gemm(true,false,bQiaA,bQiaA,1.0,0.0);
    TD->set(T);
    T->apply_denom_chem(nfrzc,naoccA,fockA);
    U->set(T);

    SharedTensor2d temp=SharedTensor2d(new Tensor2d("Ujaib", naoccA, navirA, naoccA,navirA));
    temp->sort(3214,U,1.0,0.0);
    U->scale(2.0);
    U->add(-1.0,temp);

    double mp2=U->vector_dot(TD);
    
    outfile->Printf("mp2 correction= %8.23f",mp2);


    return mp2;
}

void MP_PLUGIN::title()
{
	outfile->Printf("MP PLUGIN\n");
}
