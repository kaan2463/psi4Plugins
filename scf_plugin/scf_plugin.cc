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
#include <cmath>
#include "psi4/libfock/jk.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include <memory>

#define ABS(x) (x<0?-x:x)
#define CUT_OFF_DIAGONALIZE 1.0E10-15

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

    MintsHelper mints(reference_wavefunction_->basisset());
    Enuc = reference_wavefunction_->molecule()->nuclear_repulsion_energy(reference_wavefunction_->get_dipole_field_strength());


    nso = reference_wavefunction_->nso();
    doccA = reference_wavefunction_->doccpi()[0]+reference_wavefunction_->soccpi()[0];
    //S: overlap matrix
    Sso = std::make_shared<Tensor2d>("Overlap Matrix",nso,nso);
    Sso->set(mints.ao_overlap());
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
    SharedTensor2d Evecs = std::make_shared<Tensor2d>("Eigenvectors of Overlap Matrix",nso,nso);
    SharedTensor1d Evals = std::make_shared<Tensor1d>("Eigenvalues of overlap matrix",nso);
    SharedTensor2d SqrtEvals = std::make_shared<Tensor2d>("Sqrt of eigenvalues",nso,nso);
    SharedTensor2d X = std::make_shared<Tensor2d>("Sso^-1/2",nso,nso);

    Sso->diagonalize(Evecs,Evals,CUT_OFF_DIAGONALIZE);
    // 1/sqrt of eigen values
    for(int i=0;i<nso;i++){
        SqrtEvals->set(i,i,1.0/sqrt(Evals->get(i)));
    }

    X->back_transform(SqrtEvals,Evecs);

    SharedTensor2d Ft = std::make_shared<Tensor2d>("Transformed F matrix",nso,nso);
    SharedTensor2d Ct = std::make_shared<Tensor2d>("Transformed C matrix",nso,nso);
    SharedTensor1d Fteval = std::make_shared<Tensor1d>("Transformed F eigen vals",nso);
    SharedTensor2d C = std::make_shared<Tensor2d>("C matrix",nso,nso);
    SharedTensor2d D = std::make_shared<Tensor2d>("Density Matrix",nso,nso);

    D->zero();

    SharedTensor2d F = std::make_shared<Tensor2d>("Fock Matrix",nso,nso);
    
    double Etotal=0.0;
    double dE = __DBL_MAX__;
    int iterNumber=0;
    while(iterNumber<50 && ABS(dE)>tol_E){
        // Fock Matrix
        F->set(Hso);
        for(int n=0;n<nso;n++){
            for(int v =0;v<nso;v++){
                double item=0.0;
                for(int p=0;p<nso;p++){
                    for(int q=0;q<nso;q++){
                        item+= D->get(p,q)*(2*Eri->get(n*nso+v,p*nso+q)-Eri->get(n*nso+p,v*nso+q)); 
                    }
                }
                F->add(n,v,item);
            }
        }

        // Electronic energy
        double Eelec = 0.0;
        for(int n=0;n<nso;n++){
            for(int v =0;v<nso;v++){
               Eelec += D->get(n,v)*(Hso->get(n,v)+F->get(n,v));
            }
        }

        double EtotalOld=Etotal;
        Etotal = Eelec + Enuc;
        
        Ft->transform(F,X);

        Ft->diagonalize(Ct,Fteval,CUT_OFF_DIAGONALIZE);

        C->gemm(false,false,X,Ct,1,0);

        SharedTensor2d Dold= std::make_shared<Tensor2d>("Old D matrix",nso,nso);

        Dold->set(D);

        for(int m=0;m<nso;m++){
            for(int v =0;v<nso;v++){
                double Dmv = 0.0;
                for(int i=0;i<doccA;i++){
                    Dmv += C->get(m,i) * C->get(v,i);
                }
                D->set(m,v, Dmv);
            }
        }

        double rms= D->rms(Dold);

        dE=Etotal-EtotalOld;
        iterNumber++;

        outfile->Printf(" %d  Energy = %8.23f %8.23f %8.23f %8.23f\n",iterNumber,Eelec,Etotal,ABS(dE),rms);

    }
    return Etotal; 
}

void SCF_PLUGIN::title()
{
	outfile->Printf("Hello world!\n");
}
