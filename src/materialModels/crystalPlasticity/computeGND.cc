#include "../../../include/crystalPlasticity.h"
/* Based on the small deformation formulation*/

template <int dim>
void crystalPlasticity<dim>::computeGND(unsigned int cellID, unsigned int quadPtID, FEValues<dim>&fe_values, FullMatrix<double> &sModMat, const unsigned int &qptCt, const unsigned int &locDoft)
{
    //std::cout << "Suceessfully in the GND Function" << std::endl;
 
    // Cross product result of the slip plane normals
    std::vector<FullMatrix<double>> smTen(this->n_slip_systems, FullMatrix<double>(dim,dim));
    for(unsigned int i = 0; i < this->n_slip_systems; i++){

        smTen[i] = 0.0;

        if(dim == 3){
        smTen[i][0][1] = n_alpha[i][2];
        smTen[i][0][2] = n_alpha[i][1] * -1;
        
        smTen[i][1][0] = n_alpha[i][0] * -1;
        smTen[i][1][2] = n_alpha[i][0];

        smTen[i][2][0] = n_alpha[i][1];
        smTen[i][2][1] = n_alpha[i][2] * -1;
        }
        else if(dim == 2){
        smTen[i][0][1] = n_alpha[i][1] * -1;
        smTen[i][1][0] = n_alpha[i][0];}
    }

    // Rotation matrix of the crystal orientation
    FullMatrix<double> tempR1(dim,dim),SG(dim,dim);  // Temporary Matrices
    tempR1=0.0; SG=0.0;
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)
    rot1=rot_conv[cellID][quadPtID];
    FullMatrix<double> rotmat(dim,dim);
    rotmat=0.0;
    odfpoint(rotmat,rot1);

    // Find the GND denisties 
    for(unsigned int i = 0; i < this->n_slip_systems; i++){

        // Rotate the temporary tensor
        rotmat.mmult(tempR1, smTen[i]);
        tempR1.mTmult(SG, rotmat);

        // Get slip fraction graident
        Vector<double> slipQpts(qptCt), slipNodal(locDoft);
        for(unsigned int qd = 0; qd < qptCt; qd++){
            slipQpts[qd] = this->slipfraction_iter[cellID][qd][i];
        }
        sModMat.vmult(slipNodal, slipQpts);

        Tensor<1, dim> grad_gamma; // Gradient of slip fraction
        grad_gamma = 0.0;
            
        // Note - this is on our "single dof" per node FE, so using qptCt for dofs_per_cell
        for(unsigned int d = 0; d < locDoft; d++){
            grad_gamma += fe_values.shape_grad(d,quadPtID) * slipNodal[d];
        }
        
        

        Vector<double> curGNDVec(dim); 
        curGNDVec=0.0;
        for(unsigned int r = 0; r < dim; r++){
            for(unsigned int c = 0; c < dim; c++){
                curGNDVec[r] += SG[r][c]*grad_gamma[c];
            }
        }
        //curGNDVec[0] = SG[0][0]*grad_gamma[0] + SG[0][1]*grad_gamma[1] + SG[0][2]*grad_gamma[2];
        //curGNDVec[1] = SG[1][0]*grad_gamma[0] + SG[1][1]*grad_gamma[1] + SG[1][2]*grad_gamma[2];
        //curGNDVec[2] = SG[2][0]*grad_gamma[0] + SG[2][1]*grad_gamma[1] + SG[2][2]*grad_gamma[2];

        gndDensityPSS[cellID][quadPtID][i] = curGNDVec.l2_norm() / this->userInputs.burgVecMags[i];
        gndDensity[cellID][quadPtID] += gndDensityPSS[cellID][quadPtID][i];
        
        
        // std::cout << "Current Burgers Vector: " << this->userInputs.burgVecMags[i] << std::endl;
    }

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
