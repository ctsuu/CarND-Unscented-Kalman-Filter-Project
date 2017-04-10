#include <iostream>
#include "ukf.h"
using namespace std;


/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd (5,5);
  P_ <<   0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
         -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
          0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
         -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
         -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  // cout << P_ << endl;

  n_x_ = 5;
  n_aug_ = 7;
  n_sig_ = 2*n_aug_+1;
  lambda_ = 3 - n_aug_;


  Xsig_pred_ = MatrixXd(n_x_, n_sig_);
  // Process noise standard deviation longitudinal acceleration in m/s^2

  std_a_ = 1.5;
  std_yawdd_ = 0.5;


  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  



  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;



  // Weight initialization
  weights_ = VectorXd(n_sig_);
  weights_(0) =lambda_/(lambda_+n_aug_);
  for (int i=1; i<n_sig_; i++) {  // iterate the rest sigma points
    weights_(i) = 0.5/(n_aug_+lambda_);
  }


  is_initialized_ = false;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


  
  if(! is_initialized_) {
    time_us_ = meas_package.timestamp_;
    if(meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // cout << "Laser Measurement" << endl;
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];

      // CTRV Model State Vector
      x_ << px, py, 0, 0, 0;
      // time_us_ = meas_package.timestamp_;
    }
    else if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      // cout << "Radar Measurement" << endl;
      /**
      Convert radar from polar to cartesian coordinates
      */

      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      float px = rho * cos(phi);
      float py = rho * sin(phi);

      // CTRV Model State Vector
      x_ << px, py, 0, 0, 0;
      // time_us_ = meas_package.timestamp_;

    }
    is_initialized_=true;
  }
  else {

    // Prediction
    double delta_t = (meas_package.timestamp_ - time_us_) /1000000.0;
    time_us_= meas_package.timestamp_;

    Prediction(delta_t);


    // switch between lidar and radar measurements
    if(meas_package.sensor_type_==MeasurementPackage::LASER){
      
      UpdateLidar(meas_package);
    }
    else if(meas_package.sensor_type_==MeasurementPackage::RADAR) {
      
      UpdateRadar(meas_package);
    }


  }
  

   

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // Sigma points
  // create 7x1 augmented mean vector
  VectorXd x_aug_(7);
  x_aug_.fill(0);
  x_aug_.head(n_x_) = x_;

  // cout << x_aug_ << endl;

  
  // create 7x7 augmented covariance matrix
  MatrixXd P_aug_(n_aug_, n_aug_);
  P_aug_.fill(0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = pow(std_a_, 2);
  P_aug_(6,6) = pow(std_yawdd_, 2);

  
  // create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();


  // create augmented sigma points
  MatrixXd Xsig_aug_(n_aug_, n_sig_);

  // first column
  Xsig_aug_.col(0) = x_aug_;
  for(int i = 0; i<n_aug_; i++) {
    //Xsig_aug_.col(i+1) = x_aug_ + s * L_.col(i);
    //Xsig_aug_.col(n_aug_+i+1) = x_aug_ - s * L_.col(i);
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }
  
  // cout << Xsig_aug_ << endl;
  
  // sigma point prediction

  for(int i = 0; i<n_sig_; i++) {  //iterate over sigma points

    // Sigma Points Prediction, from 7x15 matrix to 5x15 matrix
    double p_x = Xsig_aug_.col(i)(0);
    double p_y = Xsig_aug_.col(i)(1);
    double v = Xsig_aug_.col(i)(2);
    double yaw = Xsig_aug_.col(i)(3);
    double yawd = Xsig_aug_.col(i)(4);
    double nu_a = Xsig_aug_.col(i)(5);
    double nu_yawdd = Xsig_aug_.col(i)(6);

    double px_p,py_p;
    if(fabs(yawd) > 0.001){
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );

    }else{
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);

    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;


  }


  //predicted state mean

  x_ = Xsig_pred_ * weights_;

  //cout << "Predicted state mean x_" << endl;
  //cout << x_ << endl;


  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;

    //angle normalization
    //while(x_diff_(3)> M_PI) x_diff_(3) -=2.*M_PI;
    //while(x_diff_(3)<-M_PI) x_diff_(3) +=2.*M_PI;
    x_diff_(3) = x_diff_(3)-ceil((x_diff_(3)-M_PI)/(2.*M_PI))*2.*M_PI;


    P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose();

    //cout << "Predicted state covariance P_" << endl;
    //cout << P_ << endl;

  }



}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  //transform sigma points into 2x15 measurement space
  MatrixXd Zsig_ = MatrixXd(2, 15);
  VectorXd z_pred_ = VectorXd(2);

  // measurement model
  Zsig_.row(0) = Xsig_pred_.row(0);        //px
  Zsig_.row(1) = Xsig_pred_.row(1);        //py

  z_pred_ = Zsig_ * weights_;

  //measurement covariance matrix S
  MatrixXd S_ = MatrixXd(2, 2);
  S_.fill(0);

  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(5, 2);
  Tc_.fill(0);

  

  for(int i = 0; i<n_sig_; i++) {  //iterate over sigma points
    
    // measurement residual
    VectorXd tz_diff_ = Zsig_.col(i) - z_pred_;
    S_ = S_ + weights_(i)*tz_diff_*tz_diff_.transpose();
    //cout << "Update Laser Measurement S" << endl;
    //cout << S_ << endl;
    

    // state residual
    VectorXd tx_diff_ = Xsig_pred_.col(i) - x_;
    Tc_ = Tc_ + weights_(i)*tx_diff_*tz_diff_.transpose();
    //cout << "Update Laser Measurement Tc" << endl;
    //cout << Tc_ << endl;
  }

  // Laser measurement noise covariance matrix
  MatrixXd R_laser_(2,2);
  R_laser_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;
  // cout << "Create Laser Measurement Noise matrix R" << endl;
  // cout << R_laser_ << endl;


  // add measurement noise covariance matrix
  S_ = S_ + R_laser_;
  // cout << "Update Laser Measurement Noise S" << endl;
  // cout << S_ << endl;

  // calculate karman gain K
  MatrixXd K_ = Tc_*(S_.inverse());
  // cout << "calculate karman gain K" << endl;
  // cout << K_ << endl;

  // VectorXd xo = x_;

  // Actual measurement against predicted measurement
  VectorXd A_ = meas_package.raw_measurements_-z_pred_;  

  // Update state 
  x_ = x_ + K_*A_;

  // Update NIS
  NIS_laser_ = (A_.transpose())*S_.inverse()*A_;

  // Update stat covariance
  P_ = P_ - K_*S_*K_.transpose();


}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // Transform sigma points into 3x15 measurement space
  MatrixXd Zsig_ = MatrixXd(3, n_sig_);

  // Predicted measurement mean
  VectorXd z_pred_ = VectorXd(3);


  // Copy over the sigma point predictions
  for(int i = 0; i<n_sig_; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);

    double vx = cos(yaw) * v;
    double vy = sin(yaw) * v;

    // predicted measurement rho
    Zsig_(0,i) = sqrt(px*px+py*py);    //rho
    // predicted measurement phi
    Zsig_(1,i) = atan2(py,px);           //phi

    // predicted measurement rho_dot
    if(Zsig_(0,i) > 0.001){
      Zsig_(2,i) = (px*vx + py*vy)/sqrt(px*px + py*py);     //rho_dot

    }else{

      Zsig_(2,i) = 0.0;                    //rho_dot
    }


  }

  // Calculate Predicted Measurement Mean
  z_pred_ = Zsig_ * weights_;

  // Create measurement covariance matrix S
  MatrixXd S_ = MatrixXd(3, 3);
  S_.fill(0);

  // Create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(5, 3);
  Tc_.fill(0);


  for(int i = 0; i<n_sig_; i++) {  // iterate over all sigma points
    // Calculate the difference
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    // Calculate 3x3 predicted measurement covariance S_
    S_ = S_ + weights_(i)*z_diff_*z_diff_.transpose();
    
    // Calculate state residual
    VectorXd tx_diff_ = Xsig_pred_.col(i) - x_;
    
    // Calculate measurement residual
    VectorXd tz_diff_ = Zsig_.col(i) - z_pred_;
    Tc_ = Tc_ + weights_(i)*tx_diff_*tz_diff_.transpose();

  }

  //cout << "calculate Tc" << endl;
  //cout << Tc_ << endl;



  // Create Radar Measurement Noise Covariance Matrix
  MatrixXd R_radar_(3,3);
  R_radar_ << std_radr_*std_radr_,0,0,
              0,std_radphi_*std_radphi_,0,
              0,0,std_radrd_*std_radrd_;
  
  // add measurement noise covariance matrix
  S_ = S_ + R_radar_;

 
  // Calculate Karman Gain K 
  MatrixXd K_ = Tc_*(S_.inverse());
  //cout << "calculate karman gain K" << endl;
  //cout << K_ << endl;

  // Actual measurement against predicted measurement
  VectorXd A_ = meas_package.raw_measurements_-z_pred_;  

  // Update state 
  x_ = x_ + K_*A_;

  // Update NIS
  NIS_radar_ = (A_.transpose())*S_.inverse()*A_;
  //cout << " Updated NIS_radar " << endl;
  //cout << NIS_radar_ << endl;


  // Update stat covariance
  P_ = P_ - K_*S_*K_.transpose();
  //cout << " Upated state covariance P_" << endl;
  //cout << P_ << endl;


}
