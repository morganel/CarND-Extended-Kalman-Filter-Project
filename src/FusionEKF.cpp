#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>       /* cos sin */

using namespace std;
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
  Hj_ = MatrixXd(3, 4);

  /**
  TODO:
    * Finish initializing the FusionEKF.
  */

  // measurement noise (assume independent by dimension)
  float e_laser = 0.01;
  float e_radar = 0.0225; // more noise for radar data

  R_laser_ << e_laser, 0.0,
		  	  0.0, e_laser;

  R_radar_ << e_radar, 0.0, 0.0,
  		  	  0.0, e_radar, 0.0,
			  0.0, 0.0, e_radar;

  // measurement matrix
  H_laser_ << 1.0, 0.0, 0.0, 0.0,
		  	  0.0, 1.0, 0.0, 0.0;

  // acceleration noise for Q
  noise_ax = 5.;
  noise_ay = 5.;


  // initialize ekf
  VectorXd x_in(4);
  x_in << 1,1,1,1;
  MatrixXd P_in(4,4);
  P_in << 10.0, 0.0, 0.0, 0.0,
		0.0, 10.0, 0.0, 0.0,
		0.0, 0.0, 1000.0, 0.0,
		0.0, 0.0, 0.0, 1000.0; // more uncertainty for velocity than position

  MatrixXd F_in(4,4);
  F_in << 1.0, 0.0, 1.0, 0.0,
  		  0.0, 1.0, 0.0, 1.0,
  		  0.0, 0.0, 1.0, 0.0,
  		  0.0, 0.0, 0.0, 1.0;

  MatrixXd H_in(4,4);
  MatrixXd R_in(2,2);
  MatrixXd Q_in(4,4);
  ekf_.Init(x_in, P_in, F_in, H_in, R_in, Q_in);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  cout << "starting process measurement" << endl;
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement

    cout << measurement_pack.sensor_type_ << endl;
    cout << MeasurementPackage::RADAR << endl;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];

      cout << rho*cos(phi) << endl;
      cout << rho*sin(phi) << endl;
      cout << rho_dot*cos(phi) << endl;
      cout << rho_dot*sin(phi) << endl;
      ekf_.x_ << rho*cos(phi), rho*sin(phi), rho_dot*cos(phi), rho_dot*sin(phi);

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_)/ 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
  			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
  			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

	  Hj_ = tools.CalculateJacobian(ekf_.x_);
	  ekf_.H_ = Hj_;
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
