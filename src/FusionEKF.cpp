#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_      = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  H_laser_ <<
	  1, 0, 0, 0,
	  0, 1, 0, 0;
  Hj_.fill(0);
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
    */
    // first measurement
    std::cout << "EKF: " << std::endl;
    ekf_.x_ = VectorXd(4);
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      auto rho = measurement_pack.raw_measurements_[0];
      auto theta = measurement_pack.raw_measurements_[1];
      ekf_.x_ << rho*std::cos(theta), rho*std::sin(theta), 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    previous_timestamp_ = measurement_pack.timestamp_;

    ekf_.F_ = MatrixXd::Identity(4, 4);
    ekf_.F_(0, 2) = 1.;
    ekf_.F_(1, 3) = 1.;

    ekf_.P_ = MatrixXd::Identity(4, 4);
    ekf_.P_(2, 2) = 1000.;
    ekf_.P_(3, 3) = 1000.;

    ekf_.Q_ = MatrixXd::Identity(4, 4);

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  double dt2 = dt*dt;
  double dt3 = dt2*dt;
  double dt4 = dt3*dt;
  ekf_.Q_(0, 0) = 0.25*dt4*noise_ax_;
  ekf_.Q_(1, 1) = 0.25*dt4*noise_ay_;
  ekf_.Q_(2, 2) = dt2*noise_ax_;
  ekf_.Q_(3, 3) = dt2*noise_ay_; 
  ekf_.Q_(0, 2) = 0.5*dt3*noise_ax_;
  ekf_.Q_(2, 0) = 0.5*dt3*noise_ax_;
  ekf_.Q_(1, 3) = 0.5*dt3*noise_ay_;
  ekf_.Q_(3, 1) = 0.5*dt3*noise_ay_;
  

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
