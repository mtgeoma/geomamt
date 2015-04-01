#include<iostream>
#include<iomanip>
using namespace std;

int main(int argc, char* argv[]){
    int nch=atoi(argv[1]);
    cout<<"N="<<(nch-2)*(nch-1)/2<<endl;

    for (size_t i=0;i<nch-2;i++) {
        for (size_t j=0;j<=i;j++) {
            cout<<'['<<setw(2)<<i*(i+1)/2+j<<']'<<i<<j;
        }
        cout<<endl;
    }
}
/*
int main() {
    complex<double> c=complex<double>(3,3);
    cout<<c.real()<<" "<<c.imag()<<endl;
    cout<<abs(c)<<endl;
    cout<<(pow(c.real(),2)+pow(c.imag(),2))<<endl;
    cout<<norm(c)<<endl;
    cout<<(180./3.1415927)*arg(c)<<endl;
}
*/
