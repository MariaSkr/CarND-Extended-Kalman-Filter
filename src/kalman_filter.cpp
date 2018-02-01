#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
//predict the state
void KalmanFilter::Predict() {

    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}
//update the state by using Kalman Filter equations
void KalmanFilter::Update(const VectorXd &z) {
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    PHt = P_ * Ht;
    MatrixXd S = H_ * PHt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  PHt * Si;
   // New state
   int x_size = x_.size();
   x_ = x_ + (K * y);
   MatrixXd I = MatrixXd::Identity(x_size, x_size);
   P_ = (I - K * H_) * P_;
}
//update the state by using Extended Kalman Filter equations
void KalmanFilter::UpdateEKF(const VectorXd &z) {

    // Get predicted location in polar coords.

    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

  // Check the division by 0.

    float eps = 0.000001;
    if (fabs(px) < eps && fabs(py) < eps) {
        px = eps;
        py = eps;
 }

    else if (fabs(px) < eps) {
     px = eps;
 }
  // Recalculate x object state to rho, theta, rho_dot coordinates

    float rho = sqrtf(powf(px, 2) + powf(py, 2));
    float theta = atan2f(py, px);
    float rho_dot = (px * vx + py * vy) / rho;
    VectorXd h = VectorXd(3);   // h(x_)
    h << rho, theta, rho_dot;
    VectorXd y = z - h;

    y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));


    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * PHt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  PHt * Si;

    // New state
    int x_size = x_.size();
    x_ = x_ + (K * y);
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}
