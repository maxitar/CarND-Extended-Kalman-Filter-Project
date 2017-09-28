#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  int n = estimations.size();
  for (int i = 0; i < n; ++i) {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array()*diff.array();
    rmse += diff;
  }

  rmse /= n;
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian.
  */
	MatrixXd Hj(3, 4);
	Hj.fill(0);
	double px = x_state(0);
	double py = x_state(1);
	double vx = x_state(2);
	double vy = x_state(3);
	double denom = px*px + py*py;
	if (denom < 1e-8) {
		std::cout << "CalculateJacobian() - Error - Division by zero";
	}
	//compute the Jacobian matrix
	else {
          Hj(0, 0) = px / sqrt(denom);
          Hj(0, 1) = py / sqrt(denom);
          Hj(1, 0) = -py / denom;
          Hj(1, 1) = px / denom;
          Hj(2, 0) = py*(vx*py - vy*px) / pow(denom, 1.5);
          Hj(2, 1) = px*(vy*px - vx*py) / pow(denom, 1.5);
          Hj(2, 2) = Hj(0, 0);
	  Hj(2, 3) = Hj(0, 1);
	}
	return Hj;
}
