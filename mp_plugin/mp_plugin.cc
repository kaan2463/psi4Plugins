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
#include "mp_plugin.h"
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
    std::shared_ptr<BasisSet> auxiliary = get_basisset("DF_BASIS_SCF");
    nQ = auxiliary->nbf();
    noccA = reference_wavefunction_->doccpi()[0] + reference_wavefunction_->soccpi()[0];
    noccB = reference_wavefunction_->doccpi()[0]; 
    nfrzv = reference_wavefunction_->frzvpi()[0];
    nfrzc = reference_wavefunction_->frzcpi()[0];
    naoccA=noccA-nfrzc;
    naoccB=noccB-nfrzc; 
    nmo = reference_wavefunction_->nmo();
    nvirA = nmo - noccA;
    nvirB = nmo - noccB;
    navirA=nvirA-nfrzv;
    navirB=nvirB-nfrzv;

    fockA = SharedTensor2d(new Tensor2d("fock A", nmo, nmo));
    fockB = SharedTensor2d(new Tensor2d("fock B", nmo, nmo));
    bQiaA = SharedTensor2d(new Tensor2d("DF_BASIS_CC A (Q|IA)", nQ, naoccA, navirA));
    bQiaB = SharedTensor2d(new Tensor2d("DF_BASIS_CC B (Q|IA)", nQ, naoccB, navirB));
    SharedVector epsilon_a = reference_wavefunction_->epsilon_a();
    SharedVector epsilon_b = reference_wavefunction_->epsilon_b();
    for(int i=0;i<nmo;i++){
        fockA->set(i,i,epsilon_a->get(i));
        fockB->set(i,i,epsilon_b->get(i));
    }

    // DFTensor A
    std::shared_ptr<DFTensor> DFA (new DFTensor(reference_wavefunction_->basisset(), auxiliary, reference_wavefunction_->Ca(), noccA, nvirA, naoccA, navirA, options_));
    std::shared_ptr<DFTensor> DFB (new DFTensor(reference_wavefunction_->basisset(), auxiliary, reference_wavefunction_->Cb(), noccB, nvirB, naoccB, navirB, options_));


    // Alpha OV Block
    SharedMatrix ovA = DFA->Qov();
    SharedMatrix ovB = DFB->Qov();
    bQiaA->set(ovA);
    bQiaB->set(ovB); 
    ovA.reset();
    ovB.reset();
}


double MP_PLUGIN::compute_energy()
{
    SharedTensor2d TIAJB=SharedTensor2d(new Tensor2d("TIAJB", naoccA, navirA,naoccA,navirA));
    SharedTensor2d IAJB=SharedTensor2d(new Tensor2d("(IA|JB)", naoccA, navirA, naoccA,navirA));
    SharedTensor2d Tiajb=SharedTensor2d(new Tensor2d("Tiajb",naoccB,navirB,naoccB,navirB));
    SharedTensor2d iajb=SharedTensor2d(new Tensor2d("(ia|jb)",naoccB,navirB,naoccB,navirB));
    SharedTensor2d TIAjb=SharedTensor2d(new Tensor2d("TIAjb",naoccA,navirA,naoccB,navirB));
    SharedTensor2d IAjb=SharedTensor2d(new Tensor2d("(IA|jb)",naoccA,navirA,naoccB,navirB));
    
    SharedTensor2d tempA = SharedTensor2d(new Tensor2d("temp for sort",naoccA,navirA,naoccA,navirA));
    SharedTensor2d tempB = SharedTensor2d(new Tensor2d("temp for sort",naoccB,navirB,naoccB,navirB));
    SharedTensor2d tempAB = SharedTensor2d(new Tensor2d("temp for sort",naoccA,naoccB,navirA,navirB));

    IAJB->gemm(true,false,bQiaA,bQiaA,1.0,0.0);
    tempA->sort(1432,IAJB,1.0,0.0);
    IAJB->add(-1,tempA);
    
    iajb->gemm(true,false,bQiaB,bQiaB,1.0,0.0);
    tempB->sort(1432,iajb,1.0,0.0);
    iajb->add(-1,tempB);

    TIAJB->set(IAJB);
    TIAJB->apply_denom_chem(nfrzc,naoccA,fockA);
    
    Tiajb->set(iajb);
    Tiajb->apply_denom_chem(nfrzc,naoccB,fockB);
    
    IAjb->gemm(true,false,bQiaA,bQiaB,1.0,0.0);
    tempAB->sort(1324,IAjb,1.0,0.0);
    tempAB->apply_denom_os(nfrzc, naoccA, naoccB, fockA, fockB);
    TIAjb->sort(1324,tempAB,1.0,0.0);
    
    double mp2_correction = (TIAJB->vector_dot(IAJB)+Tiajb->vector_dot(iajb))/4.0 + TIAjb->vector_dot(IAjb);
    
    outfile->Printf("MP2 CORRECTION = %8.23f \n",mp2_correction);

    return mp2_correction;
}

void MP_PLUGIN::title()
{
	outfile->Printf("MP PLUGIN\n");
}
