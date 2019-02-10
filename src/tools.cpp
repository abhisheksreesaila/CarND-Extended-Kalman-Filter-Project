#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
    

    VectorXd rmse(4);
    rmse <<  0,0,0,0;
    VectorXd hold(4);
    hold << 0,0,0,0;
    
    if(estimations.size() !=0)
        if(estimations.size() == ground_truth.size())
        {
            for (int i=0; i < estimations.size(); ++i) {
                VectorXd diff =  estimations[i] - ground_truth[i];
                VectorXd power =  diff.array() * diff.array();
                hold = hold + power;
            }
        }
    
    
    rmse = hold/ground_truth.size();
    rmse = rmse.array().sqrt();
    

    // return the result
    return rmse;
    
    
    
}







MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
    
    MatrixXd Hj(3,4);
    // recover state parameters
    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);
    
    float pxSquarePluspySquare = pow(px, 2) + pow(py, 2);
    float sqrt_pxSquarePluspYSquare = sqrt(pxSquarePluspySquare);
    float threeby2_pxSquarePluspySquare = pow(pxSquarePluspySquare, 3/2);
    
    // check division by zero
    if(px + py == 0)
        cout << "division by zero error";
    else
    {
        Hj (0,0) = px/sqrt_pxSquarePluspYSquare;
        Hj (0,1) = py/sqrt_pxSquarePluspYSquare;
        Hj (0,2) = 0;
        Hj (0,3) = 0;
        Hj (1,0) = -py/pxSquarePluspySquare;
        Hj (1,1) = px/pxSquarePluspySquare;
        Hj (1,2) = 0;
        Hj (1,3) = 0;
        Hj (2,0) = (py *(py*vx-px*vy))/threeby2_pxSquarePluspySquare;
        Hj (2,1) = (px *(px*vy-vx*py))/threeby2_pxSquarePluspySquare;
        Hj (2,2) = px/sqrt_pxSquarePluspYSquare;
        Hj (2,3) = py/sqrt_pxSquarePluspYSquare;
        
    }
    
    // compute the Jacobian matrix
    
    return Hj;
    
}
