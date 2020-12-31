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
#include "psi4/psi4-dec.h"
#include "scf_plugin.h"
#include <cmath>
#include <memory>
#include "psi4/libpsi4util/PsiOutStream.h"

double SCF_PLUGIN::rhf()
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
                for(int i=0;i<noccA;i++){
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