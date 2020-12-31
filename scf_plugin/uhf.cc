
#include "scf_plugin.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

double SCF_PLUGIN::uhf()
{
    SharedTensor2d FtA = std::make_shared<Tensor2d>("Transformed Fock Matrix for A",nso,nso); 
    SharedTensor2d CtA = std::make_shared<Tensor2d>("Transformed C Matrix for A",nso,nso);
    SharedTensor2d FtB = std::make_shared<Tensor2d>("Transformed Fock Matrix for B",nso,nso); 
    SharedTensor2d CtB = std::make_shared<Tensor2d>("Transformed C Matrix for B",nso,nso);
    SharedTensor2d CA = std::make_shared<Tensor2d>("SCF eigenvector for A",nso,nso);
    SharedTensor2d CB = std::make_shared<Tensor2d>("SCF eigenvector for B",nso,nso);
    SharedTensor2d DA = std::make_shared<Tensor2d>("Density Matrix for A",nso,nso);
    SharedTensor2d DB = std::make_shared<Tensor2d>("Density Matrix for B",nso,nso);

    FtA->transform(Hso,X);
    FtB->set(FtA);
    
    SharedTensor1d temp = std::make_shared<Tensor1d>(nso);
    //diagonalize SCF
    FtA->diagonalize(CtA,temp,CUT_OFF_DIAGONALIZE);
    CtB->set(CtA);
    //multiply S^-1/2 and tronsformed C matrixes
    CA->gemm(false,false,X,CtA,1,0);
    CB->set(CA);

    for(int m=0;m<nso;m++){
        for(int v=0;v<nso;v++){
            double item_a=0.0;
            for(int i=0;i<noccA;i++){
                item_a+=CA->get(m,i)*CA->get(v,i);
            }
            DA->add(m,v,item_a);
        }
    }

    for(int m=0;m<nso;m++){
        for(int v=0;v<nso;v++){
            double item_b=0.0;
            for(int i=0;i<noccB;i++){
                item_b+=CB->get(m,i)*CB->get(v,i);
            }
            DB->add(m,v,item_b);
        }
    }

    //Fock matrix for A
    SharedTensor2d FA = std::make_shared<Tensor2d>("Fock Matrix for A",nso,nso);
    SharedTensor2d FB = std::make_shared<Tensor2d>("Fock Matrix for B",nso,nso);
    SharedTensor2d Dtot = std::make_shared<Tensor2d>("Density Matrix Total",nso,nso);
    SharedTensor2d Dtot_old = std::make_shared<Tensor2d>("Density Matrix Total Old",nso,nso);

    double dE=0.0;
    double rms=0.0;
    double e_total=0;
    double e_total_old=0.0;
    int iter_n=0;
    do{
      
        FA->set(Hso);
        FB->set(Hso);
        
        Dtot->set(DA);
        Dtot->add(DB);

        for(int m =0;m<nso;m++){
            for(int v=0;v<nso;v++){
                double item_a =0.0;
                double item_b =0.0;
                for(int l=0;l<nso;l++){
                    for(int s=0;s<nso;s++){
                        item_a+=DA->get(l,s)*(Eri->get(m*nso+v,l*nso+s)-Eri->get(m*nso+l,v*nso+s))+DB->get(l,s)*Eri->get(m*nso+v,l*nso+s);
                        item_b+=DB->get(l,s)*(Eri->get(m*nso+v,l*nso+s)-Eri->get(m*nso+l,v*nso+s))+DA->get(l,s)*Eri->get(m*nso+v,l*nso+s);
                    }   
                }
                FA->add(m,v,item_a);
                FB->add(m,v,item_b);
            }
        }

        Dtot->set(DA);
        Dtot->add(DB);

        double e_elec = 0;
        

        for(int m=0;m<nso;m++){
            for(int v=0;v<nso;v++){
                e_elec+=(Dtot->get(m,v)*Hso->get(m,v)+DA->get(m,v)*FA->get(m,v)+DB->get(m,v)*FB->get(m,v))/2;
            }
        }
        e_total=e_elec+Enuc;

        SharedTensor1d temp = std::make_shared<Tensor1d>(nso);
        FtA->transform(FA,X);
        FtB->transform(FB,X);
        
        //diagonalize SCF
        FtA->diagonalize(CtA,temp,CUT_OFF_DIAGONALIZE);
        FtB->diagonalize(CtB,temp,CUT_OFF_DIAGONALIZE);
        //multiply S^-1/2 and tronsformed C matrixes
        CA->gemm(false,false,X,CtA,1,0);
        CB->gemm(false,false,X,CtB,1,0);
        for(int m=0;m<nso;m++){
            for(int v=0;v<nso;v++){
                double item_a=0;
                for(int i=0;i<noccA;i++){
                    item_a+=CA->get(m,i)*CA->get(v,i);
                }
                DA->set(m,v,item_a);
            }
        }

        for(int m=0;m<nso;m++){
            for(int v=0;v<nso;v++){
                double item_b=0;
                for(int i=0;i<noccB;i++){
                    item_b+=CB->get(m,i)*CB->get(v,i);
                }
                DB->set(m,v,item_b);
            }
        }

        //mean square of density matrix
        rms = Dtot->rms(Dtot_old);
        Dtot_old->set(Dtot);

        dE=e_total-e_total_old;
        e_total_old=e_total;
        iter_n++;

        outfile->Printf("%d %8.23f %8.23f %8.23f %8.23f %8.23f %8.23f\n",iter_n,e_total,e_elec,dE,rms,tol_E,tol_D);
    }while(scf_maxiter>iter_n && (ABS(dE)>tol_E || rms>tol_D));

    return e_total; 
}
