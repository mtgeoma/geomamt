#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <iomanip>

vector<string> getwords(string line);

int main(int argc, char* argv[]) {
    if(argc != 3) {
        cerr<<argv[0]<<": transforma o arquivo *.cts de Re. Im. para Amp. fase"<<endl
            <<"Uso: "<<argv[0]<<" <arquivo de entrada> <arquivo de saida>"<<endl;
        exit(1);
    }

    ifstream cal(argv[1]);
    if(!cal) {
        cerr<<"Nao pode abrir o arquivo \""<<argv[1]<<'\"'<<endl;
        exit(1);
    }

    ofstream out(argv[2]);
    if(!out) {
        cerr<<"Nao pode abrir o arquivo \""<<argv[2]<<'\"'<<endl;
        exit(1);
    }

    string line;

    // get header
    getline(cal,line);

    vector<double> freq;
    vector< complex<double> > ex,ey,hx,hy,hz;

    // separate the calibration file
    while(getline(cal,line)) {
        vector<string> words=getwords(line);
//        for(int i=0;i<words.size();i++)
//           cout<<i<<") "<<setw(12)<<words[i];
//        cout<<endl;

        if(words.size()!=12) {
            cerr<<"error reading calibration file "<<argv[1]<<endl;
            exit(1);
        }

            freq.push_back(atof(words[0].c_str()));
            complex<double> cmp(atof(words[2].c_str()),atof(words[3].c_str()));
            ex.push_back(cmp);
            cmp=complex<double>(atof(words[4].c_str()),atof(words[5].c_str()));
            ey.push_back(cmp);
            cmp=complex<double>(atof(words[6].c_str()),atof(words[7].c_str()));
            hx.push_back(cmp);
            cmp=complex<double>(atof(words[8].c_str()),atof(words[9].c_str()));
            hy.push_back(cmp);
            cmp=complex<double>(atof(words[10].c_str()),atof(words[11].c_str()));
            hz.push_back(cmp);
    }
    cal.close();

/*  out<<setw(12)<<1./freq[freq.size()-1]
       <<setw(12)<<abs(ex[freq.size()-1])
       <<setw(12)<<arg(ex[freq.size()-1])
       <<setw(12)<<abs(ey[freq.size()-1])
       <<setw(12)<<arg(ey[freq.size()-1])
       <<setw(12)<<abs(hx[freq.size()-1])
       <<setw(12)<<arg(hx[freq.size()-1])
       <<setw(12)<<abs(hy[freq.size()-1])
       <<setw(12)<<arg(hy[freq.size()-1])
       <<setw(12)<<abs(hz[freq.size()-1])
       <<setw(12)<<arg(hz[freq.size()-1])<<endl;
*/
    for(int i=freq.size()-1;i>=0;i--) {
        out<<setw(12)<<1./freq[i]
           <<setw(12)<<abs(ex[i])
           <<setw(12)<<(180./M_PI)*arg(ex[i])
           <<setw(12)<<abs(ey[i])
           <<setw(12)<<(180./M_PI)*arg(ey[i])
           <<setw(12)<<abs(hx[i])
           <<setw(12)<<(180./M_PI)*arg(hx[i])
           <<setw(12)<<abs(hy[i])
           <<setw(12)<<(180./M_PI)*arg(hy[i])
           <<setw(12)<<abs(hz[i])
           <<setw(12)<<(180./M_PI)*arg(hz[i])<<endl;
    }
    out.close();
}

// get words separated by commas in line
vector<string> getwords(string line) {
    vector<string> words;
    while(line.find(",")!=string::npos) {
        string word=line.substr(0,line.find(","));
        words.push_back(word);
        line=line.substr(line.find(",")+1,string::npos);
    }
    return words;
}
