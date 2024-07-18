#include <iostream>

#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_.fill(0.0);

  // initial predicted state vector
  x_pred_ = VectorXd(n_x_);
  // x_pred_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  P_.fill(0.0);

  // initial covariance matrix
  P_pred_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.2; // 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3; // 30;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;

  // set design parameter
  lambda_ = 3;

  // Xsig_ initialize
  Xsig_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_.fill(0.0);

  // Xsig_pred_ initialize
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);

  // Weights initialize
  weights_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2 * n_aug_).setConstant(1 / (2 * (lambda_ + n_aug_)));
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (!is_initialized_)
  {
    // cout << "Kalman Filter Initialization " << endl;
    // different x initialization depending on type of measurement
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {

      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];

      // convert radar from polar to cartesian coords
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx * vx + vy * vy);

      // set the state with the initial location and zero velocity
      // x_ << px, py, v, 0, 0;
      x_ << px, py, v, rho, rho_dot;
      // P_ << std_radr_ * std_radr_, 0, 0, 0, 0,
      //     0, std_radr_ * std_radr_, 0, 0, 0,
      //     0, 0, std_radrd_ * std_radrd_, 0, 0,
      //     0, 0, 0, std_radphi_ * std_radphi_, 0,
      //     0, 0, 0, 0, std_radphi_ * std_radphi_;

      P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 0.5, 0, 0,
          0, 0, 0, 0.5, 0,
          0, 0, 0, 0, 0.5;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];

      x_ << px, py, 0, 0, 0;

      // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
      //     0, std_laspy_ * std_laspy_, 0, 0, 0,
      //     0, 0, 1, 0, 0,
      //     0, 0, 0, 1, 0,
      //     0, 0, 0, 0, 1;

      P_ << 0.5, 0, 0, 0, 0,
          0, 0.5, 0, 0, 0,
          0, 0, 100, 0, 0,
          0, 0, 0, 100, 0,
          0, 0, 0, 0, 10;
    }
    // // set initial covariance matrix
    // P_ << std_laspx_ * std_laspx_, 0, 0, 0, 0,
    //     0, std_laspy_ * std_laspy_, 0, 0, 0,
    //     0, 0, std_radr_ * std_radr_, 0, 0,
    //     0, 0, 0, std_radphi_ * std_radphi_, 0,
    //     0, 0, 0, 0, std_radrd_ * std_radrd_;

    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // compute the time elapsed between the current and previous measurements
  // dt - expressed in seconds
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Process Radar measurement
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    // Process LiDAR measurement
    UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double delta_t)
{
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  GenerateAugmentedSigmaPoints();
  PredictSigmaPoints(delta_t);
  PredictMeanAndCovariance();

  // Update the state and covariance matrix with the predicted values
  x_ = x_pred_;
  P_ = P_pred_;
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  // Lidar measurement dimension is 2: px and py
  int n_z = 2;

  // Extract the measurement as a vector
  VectorXd z = VectorXd(n_z);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

  // Measurement matrix for Lidar
  MatrixXd H = MatrixXd(n_z, n_x_);
  H << 1, 0, 0, 0, 0,
      0, 1, 0, 0, 0;

  // Measurement noise covariance matrix for Lidar
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;

  // Calculate the residual
  VectorXd y = z - H * x_;

  // Calculate the measurement covariance matrix S
  MatrixXd S = H * P_ * H.transpose() + R;

  // Calculate the Kalman gain K
  MatrixXd K = P_ * H.transpose() * S.inverse();

  // Update the state mean and covariance matrix
  x_ = x_ + K * y;
  P_ = P_ - K * H * P_;

  // Print the result
  std::cout << "Updated state x: " << std::endl
            << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl
            << P_ << std::endl;
}

void UKF::NormalizeAngle(double &angle)
{
  while (angle > M_PI)
    angle -= 2.0 * M_PI;
  while (angle < -M_PI)
    angle += 2.0 * M_PI;
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // set number measurement variables
  int n_z = 3;

  // initial predicted measurement vector
  z_pred_ = VectorXd(n_z);
  z_pred_.fill(0.0);

  // initial predicted measurment covariance
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);

  // create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z, 2 * n_aug_ + 1);

  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    double px_p = Xsig_pred_(0, i);
    double py_p = Xsig_pred_(1, i);
    double v_p = Xsig_pred_(2, i);
    double psi_p = Xsig_pred_(3, i);

    double rho = sqrt(px_p * px_p + py_p * py_p);
    double phi = atan2(py_p, px_p);
    double rho_dot = (px_p * v_p * cos(psi_p) + py_p * v_p * sin(psi_p)) / rho;

    Zsig_.col(i) << rho, phi, rho_dot;
  }

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    z_pred_ += weights_(i) * Zsig_.col(i);
  }

  // calculate innovation covariance matrix S
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    NormalizeAngle(z_diff(1)); // z(1) is angle phi

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  // for (int i = 0; i < 2 * n_aug_ + 1; i++)
  // {
  //   S += weights_(i) * (Zsig_.col(i) - z_pred_) * (Zsig_.col(i) - z_pred_).transpose();
  // }
  S += R;

  // read in actual measurement z (should be k+1)
  VectorXd z = VectorXd(n_z);

  double rho = meas_package.raw_measurements_[0];
  double phi = meas_package.raw_measurements_[1];
  double rho_dot = meas_package.raw_measurements_[2];

  // set the state with the initial location and zero velocity
  z << rho, phi, rho_dot;

  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  // Calculate cross correlation matrix
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    NormalizeAngle(x_diff(3)); // Normalize the angle

    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    NormalizeAngle(z_diff(1)); // Normalize the angle

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // // calculate cross correlation matrix
  // for (int i = 0; i < 2 * n_aug_ + 1; i++)
  // {
  //   Tc += weights_(i) * (Xsig_pred_.col(i) - x_) * (Zsig_.col(i) - z_pred_).transpose();
  // }

  // calculate Kalman gain K;
  MatrixXd K = MatrixXd(n_x_, n_z);
  K = Tc * S.inverse();

  // update state mean and covariance matrix
  // Update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_;
  NormalizeAngle(z_diff(1)); // Normalize the angle

  // x_ = x_ + K * (z - z_pred_);
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  /**
   * Student part end
   */

  // print result
  std::cout << "Updated state x: " << std::endl
            << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl
            << P_ << std::endl;

  // // print result
  // std::cout << "z_pred: " << std::endl
  //           << z_pred_ << std::endl;
  // std::cout << "S: " << std::endl
  //           << S << std::endl;
}

void UKF::GenerateAugmentedSigmaPoints()
{
  // create augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);

  // create augmented state covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  // create augmented sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // create augmented covariance matrix
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
      0, std_yawdd_ * std_yawdd_;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug.bottomRightCorner(2, 2) = Q;

  // create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  // set remaining sigma points
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  Xsig_ = Xsig_aug;
}

void UKF::PredictSigmaPoints(double delta_t)
{
  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for each column from Xsig_aug for readability
    double p_x = Xsig_(0, i);
    double p_y = Xsig_(1, i);
    double v = Xsig_(2, i);
    double yaw = Xsig_(3, i);
    double yawd = Xsig_(4, i);
    double nu_a = Xsig_(5, i);
    double nu_yawdd = Xsig_(6, i);

    // predicted state values
    double px_p;
    double py_p;

    // avoid division by zero, sub in the sigma points into CTRV process model
    if (fabs(yawd) > 0.001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }

    double v_p = v; // v_dot is 0 so v hasn't changed, assumption of CTRV
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd; // constant turn rate too, so yawd const

    // add the process noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;

    yaw_p = yaw_p + 0.5 * delta_t * delta_t * nu_yawdd;
    // Normalize the yaw angle
    NormalizeAngle(yaw_p);
    yawd_p = yawd_p + delta_t * nu_yawdd;

    // write predicted sigma points into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }
  // print result
  // std::cout << "Xsig_pred = " << std::endl
  //           << Xsig_pred_ << std::endl;
}

void UKF::PredictMeanAndCovariance()
{
  // reset predicted mean and covariance to zero
  x_pred_.fill(0.0);
  P_pred_.fill(0.0);

  for (int i = 0; i < weights_.size(); i++)
  {
    x_pred_ += weights_(i) * Xsig_pred_.col(i);
  }

  // predict state covariance matrix
  // for (int i = 0; i < weights_.size(); i++)
  // {
  //   P_pred_ += weights_(i) * (Xsig_pred_.col(i) - x_pred_) * (Xsig_pred_.col(i) - x_pred_).transpose();
  // }
  for (int i = 0; i < weights_.size(); i++)
  {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;

    // Normalize the angles
    NormalizeAngle(x_diff(3)); // Normalize yaw angle
    NormalizeAngle(x_diff(4)); // Normalize yaw rate angle (if applicable)

    P_pred_ += weights_(i) * x_diff * x_diff.transpose();
  }
  // // print result
  // std::cout << "Predicted state" << std::endl;
  // std::cout << x_pred_ << std::endl;
  // std::cout << "Predicted covariance matrix" << std::endl;
  // std::cout << P_pred_ << std::endl;
}