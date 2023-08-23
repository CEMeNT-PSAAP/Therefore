#include <cmath>
#include <vector>

# define e 2.7182818

class mms{
    public:
        double A;
        double B;
        double C;
        double D;
        double F;
        double v1;
        double v2;

        double group1afCONT(double mu, double x, double t, double sigma){
            //the continuous angular flux function for group 1 //
            return ( A*cos(mu) + B*pow(x,2)*t ) ;
        }

        double group1afUNINT(){
            return ( 1.0 );
        }

        double group1sourceUNINT(double mu, double x, double t, double sigma){
            //the un-evaluated integral of the source term used to take our continuous source function into the proper democratization for SCB-TDMB//
            return ( (2*(B*sigma*t*v1+B)*pow(x,3))/(3*v1)+2*B*mu*t*pow(x,2)+2*sigma*acos(mu)*x );
        }

        std::vector<double> group1source(double x, double dx, double t, double dt, double mu, double sigma){
            //actually evaluating the integrals using the un-evaluated term and returning them as values//
            double x_m = x-dx/2;
            double x_i = x;
            double x_p = x+dx/2;

            double t_p = t+dt/2;
            double t_m = t-dt/2;

            // evaluating the integral over dx_i-1/2
            double Qa = group1sourceUNINT(mu,x_i,t_p, sigma) - group1sourceUNINT(mu,x_m,t_p,sigma);
            // evaluating the integral over dx_i+1/2
            double Qb = group1sourceUNINT(mu,x_p,t_p, sigma) - group1sourceUNINT(mu,x_i,t_p,sigma);
            // finding the time averaged integral dx_i-1/2 
            double Qc = ( Qa + (group1sourceUNINT(mu,x_i,t_m, sigma) - group1sourceUNINT(mu,x_m,t_m,sigma)) ) / 2;
            // finding the time averaged integral dx_i+1/2 
            double Qd = ( Qb +  group1sourceUNINT(mu,x_p,t_m, sigma) - group1sourceUNINT(mu,x_i,t_m,sigma) ) / 2;

            return( std::vector<double> {Qa, Qb, Qc, Qd} );
        }

        //
        // GROUP 2
        //
        
        double group2afCONT(double mu, double x, double t, double sigma){
            //the continuous angular flux function for group 2 //
            return ( C * pow(e, -F*x*t ) + D*pow(mu,2) ) ;
        }

        double group2afUNINT(){
            return ( 1.0 );
        }

        double group2sourceCONT(double mu, double x, double t, double sigma){
            //continuous source from MMS//
            return( (2*pow(e ,-F*t*x ) * (D*mu*2*sigma*v2*pow( e,F*t*x )) - C*F*x - C*(F*mu*t-sigma)*v2) / v2 );
        }

        double group2sourceUNINT(double mu, double x, double t, double sigma){
            //the un-evaluated integral (dx) of the source term used to take our continuous source function into the proper discretization for SCB-TDMB//
            //(2*((self.C*(self.F*t*x+1)*math.exp(-self.F*t*x))/(self.F*t**2)-(self.C*sigma*self.v*math.exp(-self.F*t*x))/(self.F*t)+self.C*mu*self.v*math.exp(-self.F*t*x)+self.D*mu**2*sigma*self.v*x))/self.v )
            return ( (2*(C*(F*t*x+1)*pow(e,-F*t*x))/(F*pow(t,2)) - (C*sigma*v2*pow(e,-F*t*x))/(F*t) + C*mu*v2*pow(e,-F*t*x)+D*pow(mu,2)*sigma*v2*x) / v2 );
        }

        std::vector<double> group2source(double x, double dx, double t, double dt, double mu, double sigma){
            //actually evaluating the integrals using the un-evaluated term and returning them as values//
            double x_m = x-dx/2;
            double x_i = x;
            double x_p = x+dx/2;

            double t_p = t+dt/2;
            double t_m = t-dt/2;

            double Qa = group2sourceUNINT(mu,x_i,t_p, sigma) - group2sourceUNINT(mu,x_m,t_p,sigma);
            double Qb = group2sourceUNINT(mu,x_p,t_p, sigma) - group2sourceUNINT(mu,x_i,t_p,sigma);
            double Qc = ( Qa + (group2sourceUNINT(mu,x_i,t_m, sigma) - group2sourceUNINT(mu,x_m,t_m,sigma)) ) / 2;
            double Qd = ( Qb +  group2sourceUNINT(mu,x_p,t_m, sigma) - group2sourceUNINT(mu,x_i,t_m,sigma) ) / 2;

            return( std::vector<double> {Qa, Qb, Qc, Qd} );
        }
        
        double group2af(double x, double dx, double t, double dt, double mu, double sigma){
            //The proper spatial integral and time edge values of the chosen continuous angular flux//
            return( 1.0 );
        }

};