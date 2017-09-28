#include "kalman_filter.h"
#include <iostream>
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
    * predict the state
  */
  x_ = F_*x_;
  P_ = F_*P_*F_.transpose() + Q_;
}

void KalmanFilter::UpdateState(const VectorXd &y) {
  MatrixXd PHt = P_*H_.transpose();
  MatrixXd S = H_*PHt + R_;
  MatrixXd K = PHt*S.inverse();

  x_ = x_ + K*y;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  P_ = (I - K*H_)*P_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_*x_;
  VectorXd y = z - z_pred;
  UpdateState(y);
}

VectorXd toPolar(const VectorXd& x_state) {
  auto px = x_state(0);
  auto py = x_state(1);
  auto vx = x_state(2);
  auto vy = x_state(3);
  auto rho = std::sqrt(px*px + py*py);
  auto theta = std::atan2(py, px);
  auto vrho = (px*vx + py*vy) / rho;
  VectorXd x_polar(3);
  x_polar << rho, theta, vrho;
  return x_polar;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    * update the state by using Extended Kalman Filter equations
  */
  auto z_pred = toPolar(x_);
  VectorXd y = z - z_pred;
  y(1) = fmod(y(1), 2. * 3.14159265);
  UpdateState(y);
}
