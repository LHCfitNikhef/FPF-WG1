//Change the tolerance for PDF profiling in all subdirectories
//  g++ --std=gnu++0x -o changeToleranceEXE changeTolerance.cxx
//  ./changeToleranceEXE 
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

bool inStreamOpen(ifstream &in, string filename) {
    in.open(filename.c_str());
    if (!in.is_open()) {
        cout << "Can't open " << filename << endl;
        return false;        
    } else return true;
}

bool outStreamOpen(fstream &out, string filename) {
    out.open(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!out.is_open()) {
        cout << "Can't open " << filename << endl;
        return false;        
    } else return true;
}

//A function to find the line containing the substring "changeme" in the
//"filepath" file. If sub == true, only the substring is changed. Otherwise the 
//line is replaced entirely by "replacement".
int reWrite(string filepath, string changeme, string replacement, bool sub) {

    stringstream sstream;
    ifstream in;
    fstream out;
    string line, lineBegin, lineEnd;
    
    if (!inStreamOpen(in,filepath)) return -1;
    if (!outStreamOpen(out,filepath+"tmp")) {in.close();  return -1;}

    while (getline(in,line)) {
        if (line.find(changeme)!=string::npos) {
            if (sub) {
                lineBegin = line;    lineEnd = line;
                lineBegin.erase(line.find(replacement));
                lineEnd.erase(0,line.find(replacement)+replacement.size());
                out << lineBegin << replacement << lineEnd << endl;
            } else out << replacement << endl;
        } else out << line << endl;
    }
    
    in.close();
    out.close();
    
    //Replace previous file with new
    sstream.str("");  sstream.clear();
    sstream << "mv " << filepath << "tmp " << filepath;
    system(sstream.str().c_str());
    
    return 1;
}

//List all subdirectory paths under dir in a returned vector of strings
vector<string> getSubDirs(string dir) {
    stringstream sstream;
    sstream << "ls -d " << dir << "/*/ > lstmp.txt";
    system(sstream.str().c_str());
    ifstream in;
    inStreamOpen(in,"lstmp.txt");
    string stmp, line;
    vector<string> ret;
    while (getline(in,line)) {
        sstream.str("");  sstream.clear();
        sstream << line << endl;
        while (sstream >> stmp) {
            stmp.pop_back();  //Remove "/"
            ret.push_back(stmp);
        }
    }
    system("rm lstmp.txt");
    return ret;
}

int main() {
    stringstream sstream;
    double tol = 3.1622777;  //sqrt(10) to correspond to PDF4LHC21
    //double tol = 5.7445626;  //sqrt(33) to correspond to EPPS21
    sstream << "  scalePdfs: "
            << tol
            << "         "
            << "# rescale (divide) all PDF eigenvectors by the given factor";
    string scalePdfsStr = sstream.str();    

    //Fetch paths to all run directories
    vector<string> dirs = getSubDirs(".");
    vector<string> subdirs;
    for (string d : dirs) for (string sd : getSubDirs(d)) subdirs.push_back(sd);

    //Replace the scalePdfs line in each parameters.yaml
    for (string sd : subdirs) {
        reWrite(sd+"/parameters.yaml","scalePdfs: ",scalePdfsStr,false);
    }
    return 1;
}
