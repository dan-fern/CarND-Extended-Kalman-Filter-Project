#include "kalman_filter.h"


using Eigen::MatrixXd;
using Eigen::VectorXd;


// NOTE: Eigen lib does not initialize VectorXd or MatrixXd objects with zeros.


/** Constructor */
KalmanFilter::KalmanFilter( )
{

}


/** Destructor */
KalmanFilter::~KalmanFilter( )
{

}


/** Init Initializes Kalman filter
    * @param x_in Initial state
    * @param P_in Initial state covariance
    * @param F_in Transition matrix
    * @param H_in Measurement matrix
    * @param R_in Measurement covariance matrix
    * @param Q_in Process covariance matrix
*/
void KalmanFilter::Init( VectorXd &x_in,
                         MatrixXd &P_in,
                         MatrixXd &F_in,
                         MatrixXd &H_in,
                         MatrixXd &R_in,
                         MatrixXd &Q_in )
{
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
}


/** Predicts the state and state covariance using the process model
    * @param delta_T Time between k and k+1 in s
*/
void KalmanFilter::Predict( )
{

    //MatrixXd Ft = F_.transpose( );

    x_ = F_ * x_;

    P_ = F_ * P_ * F_.transpose( ) + Q_;

}


/** Updates the state by using standard Kalman Filter equations
    * @param z The measurement at k+1
*/
void KalmanFilter::Update( const VectorXd &z )
{

    VectorXd y = z - ( H_ * x_ );

    MatrixXd S = H_ * P_ * H_.transpose() + R_;

    MatrixXd K = P_ * H_.transpose( ) * S.inverse( );

    //new estimate
    x_ = x_ + ( K * y );

    MatrixXd I = MatrixXd::Identity( x_.size( ), x_.size( ) );

    P_ = ( I - K * H_ ) * P_;

}


/** Updates the state by using Extended Kalman Filter equations
    * @param z The measurement at k+1
*/
void KalmanFilter::UpdateEKF( const VectorXd &z )
{
    //state parameters
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

    float rho = sqrt( px * px + py * py );

    if( fabs( rho ) < 0.00001 )
    {
        px += 0.001;
        py += 0.001;
        rho = sqrt( px * px + py * py );
    }

    float phi = atan2( py, px );

    float rho_dot = ( px * vx + py * vy ) / rho;

    VectorXd z_pred( 3 );
    z_pred << rho, phi, rho_dot;

    VectorXd y = z - z_pred;

    //normalize angle
    while( y( 1 ) > M_PI )
    {
        y( 1 ) -= 2 * M_PI;
    }
    while( y( 1 ) < -M_PI )
    {
        y( 1 ) += 2 * M_PI;
    }

    MatrixXd S = H_ * P_ * H_.transpose( ) + R_;

    MatrixXd K = P_ * H_.transpose( ) * S.inverse( );

    //new estimate
    x_ = x_ + ( K * y );

    MatrixXd I = MatrixXd::Identity( x_.size( ), x_.size( ) );

    P_ = ( I - K * H_ ) * P_;

}
