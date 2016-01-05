#include <iomanip>
#include <iostream>
#include "TLorentzVector.h"
#include "TVector3.h"

using namespace std;

int main(){

/* UNSMEARED CM EVENTS
1 -11 0 0 0 0 0.0381544 0.0367909 -0.0196287 0.0565232 0.000511
1 11 0 0 0 0 0.00613642 0.00961213 0.000183961 0.0114168 0.000511
1 22 0 0 0 0 -0.0442908 -0.046403 0.0194448 0.06703 9.31323e-10
3
1 -11 0 0 0 0 0.0100433 -0.010453 0.00759892 0.0163749 0.000511
1 11 0 0 0 0 0.0361004 -0.00384014 0.0377222 0.0523566 0.000511
1 22 0 0 0 0 -0.0461437 0.0142932 -0.0453211 0.0662385 -1.31709e-09
*/

/* Smeared Lab Events
1 -11 0 0 0 0 2.24861029 -1.61543235 -5.05352421 5.76228923 0.000510999995
1 11 0 0 0 0 0.314349887 -0.221386091 -0.704369765 0.802474099 0.000511
1 22 0 0 0 0 1.29650502 -1.04847936 -3.03537639 3.46319853 0
3
1 -11 0 0 0 0 0.275098331 0.771824966 0.169790892 0.836792733 0.000511
1 11 0 0 0 0 1.60099036 4.6169779 0.994288045 4.986809 0.000510999999
1 22 0 0 0 0 1.32239361 4.06489203 0.791782425 4.34729703 5.96046448e-08
*/


TLorentzVector x1,x2,x3,xpi;
x1.SetXYZM(2.24861029, -1.61543235, -5.05352421 ,0.000510999995);
x2.SetXYZM(0.314349887, -0.221386091, -0.704369765, 0.000511);
x3.SetXYZM(1.29650502, -1.04847936, -3.03537639, 0);
xpi = x1+x2+x3;

double gamma,beta;
gamma = xpi.E() / xpi.M();
beta = xpi.P() / xpi.E();

double E1cm,P1cm, E2cm,P2cm,  E3cm, P3cm;
TLorentzVector x1cm,x2cm,x3cm;

TLorentzVector x1c,x2c,x3c;

x1cm.SetXYZM(0.0381544, 0.0367909, -0.0196287, 0.000511);
x2cm.SetXYZM(0.00613642, 0.00961213, 0.000183961, 0.000511);
x3cm.SetXYZM(-0.0442908, -0.046403, 0.0194448, 9.31323e-10);

gamma = x1.E()/ x1.M();
beta = x1.P() / x1.E();
E1cm = gamma*(x1.E() - beta* x1.P()) ;/// ( gamma * (1- beta*beta) );
P1cm = gamma*(x1.P() - beta* x1.E()); /// ( gamma * (1- beta*beta) );
x1c.SetXYZM(P1cm, x1.Py(),x1.Pz(),x1.M());


gamma = x2.E()/ x2.M();
beta = x2.P() / x2.E();
E2cm = gamma*(x2.E() - beta* x2.P()); /// ( gamma * (1- beta*beta) );
P2cm = gamma*(x2.P() - beta* x2.E()); /// ( gamma * (1- beta*beta) );
x2c.SetXYZM(P2cm, x2.Py(), x2.Pz(), x2.M());

beta = x3.P() / x3.E();
E3cm = (x3.E() - beta* x3.P()); /// ( gamma * (1- beta*beta) );
P3cm = (x3.P() - beta* x3.E()); /// ( gamma * (1- beta*beta) );

TVector3 test(-x1.Px(),-x1.Py(),-x1.Pz() );

TLorentzVector test4v;
test4v.SetXYZM(x1.Px(),x1.Py(),x1.Pz(),x1.M());

x1.Boost(test);

cout<<x1.Px()<<" "<<x1.Py()<<endl;
cout<<test.X()<<endl;
cout<<test4v.Px()<<endl;;
cout<<test4v.E()<<" "<<test4v.P()<<endl;


cout<<"calculated E1*="<<x1c.E()<<"  CM E1*="<<x1cm.E()<<"  calculated P1*="<<x1c.P()<<"  CM P1*="<<x1cm.P()<<endl;

cout<<"calculated E2*="<<x2c.E()<<"  CM E2*="<<x2cm.E()<<"  calculated P2*="<<x2c.P()<<"  CM P2*="<<x2cm.P()<<endl;

cout<<"calculated E3*="<<E3cm<<"  CM E3*="<<x3cm.E()<<"  calculated P3*="<<P3cm<<"  CM P3*="<<x3cm.P()<<endl;

cout<<"CM calculated P* sum = "<<P1cm+P2cm+P3cm<<endl;
}
