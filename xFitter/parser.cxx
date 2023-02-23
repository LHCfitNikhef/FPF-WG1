//  g++ --std=gnu++0x -o parserEXE parser.cxx
//  ./parserEXE
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

//AUX functions to open file I/O streams
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

//AUX function to check if dir exists, make it recursively if not
void dirCheck(string dir) {

    //Deal w/ terminating '/' and save working copy of dir string
    if (dir.back()=='/') dir.erase(dir.size()-1);
    string path = dir;
    
    //Store system command state, 0 if successful
    int s=system(("ls "+path).c_str());
    if (s!=0) {
        if (path.find("/")!=string::npos) {  //See if previous dir exists
            reverse(path.begin(),path.end());
            path.erase(0,path.find("/"));
            reverse(path.begin(),path.end());
            if (system(("ls "+path).c_str())!=0) dirCheck(path);  //Recursion
        }    
        s = system(("mkdir "+dir).c_str());    
        if (s==0) cout << "-> made " << dir << endl;
    }
}

//AUX function for reliable string to double conversion
double str2d(string str) {
    stringstream sstream;
    sstream << str;
    double d;
    sstream >> d;
    return d;
}

//AUX function to read column "iRead" into "arr" from datafile "input" */
bool readCols(string input, int iRead, vector<double> &arr) {
    ifstream in;    
    if (!inStreamOpen(in, input)) return false;        
    string line, stmp;
    stringstream sstream;
    double dtmp;
    while (getline(in,line)) {
        sstream.str("");    sstream.clear();    //Reinit
        sstream << line;
        for (int i=0; i!=iRead; ++i)  sstream >> stmp;
        sstream >> dtmp;
        arr.push_back(dtmp);
    }
    in.close();    
    return true;
} //readCols

//AUX function to replace all instances of substring "oldsub" by "newsub" in
//the input "instr". Returns a new string.
string replace(string instr, string oldsub, string newsub) {
    string ret = instr;
    if (newsub.find(oldsub)!=string::npos) {
        cout << "WARNING: replace can't handle the following:" <<endl;
        cout << "    instr  = " << instr  << endl;
        cout << "    oldsub = " << oldsub << endl;
        cout << "    newsub = " << newsub << endl;
    } else {
        while (ret.find(oldsub)!=string::npos) {
            string retstart = ret;
            string retend   = ret;
            retstart.erase(retstart.find(oldsub));
            retend.erase(0,retend.find(oldsub)+oldsub.size());
            ret = retstart + newsub + retend;
        }
    }
    return ret;
}

//AUX function to split a string at whitespaces, return substrings in a vector
vector<string> stringSplit(string str) {
    vector<string> ret;
    stringstream sstream;
    sstream << str;
    string stmp;
    while (sstream >> stmp) ret.push_back(stmp);
    return ret;
}

//AUX function to write an xFitter table
//Param  mdir     Master directory above datatables and grids
//       sdir     Subdirectory under mdir where to write these tables
//       expname  Experiment name
//       nuID     Neutrino PDG ID as str, start w/ "+-" to include both nu & antinu
//       index    Dataset must be assigned a unique index
//       binF     Flags (int 1 or 0) if bin used in fit or not
//       x,Q2     Vec. of x and Q^2 values for each data point
//       xsec     Vec. of cross-sections
//       stat     Vec. of statistical uncertainties for each data point
//       uncor    Vec. of uncorrelated uncertainties
//       sStr     Names of correlated systematic uncertainty sources
//       syst     Matrix of corr. syst. unc. values
void xFtableWriter(string mdir, string sdir,
                   string expname, string nuID_in, int index, 
                   vector<int>    binF_in, 
                   vector<double> x, 
                   vector<double> Q2, 
                   vector<double> xsec, 
                   vector<double> stat, 
                   vector<double> uncor, 
                   vector<string> sStr, 
                   vector< vector<double> > syst)
{
    string nuID = replace(nuID_in,"-+","+-");
    string expID = expname + "_" + nuID;

    //LaTeX and ROOT style strings for experiment identification
    string exptag = replace(expID, "v",     "$\\nu$");
    exptag        = replace(exptag,"_14",  " $\\nu_{\\mu}$");
    exptag        = replace(exptag,"_-14", " $\\bar{\\nu}_{\\mu}$");
    exptag        = replace(exptag,"_+-14"," $\\nu_{\\mu}+\\bar{\\nu}_{\\mu}$");
    string expstr = replace(exptag,"$","");
    expstr        = replace(expstr,"\\","#");
    vector<int> binF = binF_in;
    
    //Check if uncor, stat unc & all xsec (or those w/ binflag 1) given
    bool useUncor = uncor.size() > 0;
    bool useStat  = stat.size()  > 0;
    bool xsecWriteAll = xsec.size() == binF.size();

    dirCheck(mdir+sdir);
        
    fstream out;
    outStreamOpen(out,mdir+sdir+expID+"-thexp.dat");
    out << "* " << exptag << " CC DIS pseudodata\n";
    out << "&Data\n";
    out << "   Name = '" << exptag << "'\n";
    out << "   IndexDataset = " << index << "\n"; //Give set a unique ID
    out << "   Reaction = 'CC nup'\n\n";
    out << "   NData = " << x.size() << "\n";
    int Nerr = sStr.size();
    if (useStat) ++Nerr;
    if (useUncor) ++Nerr;
    int Ncol = 1 + 2 + 1 + Nerr;
    out << "   NColumn = " << Ncol << "\n";
    out << "   ColumnType = 'Flag',2*'Bin','Sigma'";
    if (Nerr>0) out << "," << Nerr << "*'Error'";
    out << "\n";
    out << "   ColumnName = 'binFlag','x','Q2','Sigma'";
    if (useStat)  out << ",'stat'";
    if (useUncor) out << ",'uncor'";
    for (auto s : sStr) out << ",'" << s << "'";
    out << "\n";
    out << "   TheoryType = 'expression\'\n";     
    if (nuID.find("+-")!=string::npos) {
        out << "   TermType   = 'reaction','reaction'\n";
        out << "   TermName   = 'P','A'\n";
        out << "   TermSource = 'PineAPPL','PineAPPL'\n";
        out << "   TermInfo   = 'GridName="
            << mdir << "grids/"
            << replace(expID,"+-","")
            << "/nu_A_1-XSFPFCC.pineappl.lz4',\n";
        out << "                'GridName="
            << mdir << "grids/"
            << replace(expID,"+-","m")
            << "/nub_A_1-XSFPFCC.pineappl.lz4'\n";
        out << "   TheorExpr  = 'P+A'\n\n";
    } else {        
        out << "   TermType   = 'reaction'\n";
        out << "   TermName   = 'P'\n";
        out << "   TermSource = 'PineAPPL'\n";
        out << "   TermInfo   = 'GridName=";
        if (nuID.find("-")!=string::npos) {
            out << mdir << "grids/"
                << replace(expID,"-","m")
                << "/nub_A_1-XSFPFCC.pineappl.lz4'\n";
        } else {            
            out << mdir << "grids/"
                << replace(expID,"-","m")
                << "/nu_A_1-XSFPFCC.pineappl.lz4'\n";
        }
        out << "   TheorExpr  = 'P'\n\n";
    }
    out << "   Percent = " << Nerr << "*True\n";
    out << "&End\n";
    out << "&PlotDesc\n";
    out << "   PlotN = 4\n";
    out << "   PlotDefColumn = 'Q2'\n";
    out << "   PlotVarColumn = 'x'\n";
    out << "   PlotDefValue = 3., 5., 11., 110., 1100.\n";
    out << "   PlotOptions(1)  = 'Experiment:"
        << expstr << "@ExtraLabel:#nu p (CC)";
    out <<   " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}";
    out <<   " @Title:Q^{2} > 4 @Ymin:0.01@Xlog@Ylog'\n";
    out << "   PlotOptions(2)  = 'Experiment:"
        << expstr << "@ExtraLabel:#nu p (CC)";
    out <<   " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}";
    out <<   " @Title:Q^{2} > 10 @Ymin:0.01@Xlog@Ylog'\n";
    out << "   PlotOptions(3)  = 'Experiment:"
        << expstr << "@ExtraLabel:#nu p (CC)";
    out <<   " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}";
    out <<   " @Title:Q^{2} > 100 @Ymin:0.01@Xlog@Ylog'\n";
    out << "   PlotOptions(4)  = 'Experiment:"
        << expstr << "@ExtraLabel:#nu p (CC)";
    out <<   " @XTitle: x @YTitle: d^{2}#sigma/dxdQ^{2}";
    out <<   " @Title:Q^{2} > 1000 @Ymin:0.01@Xlog@Ylog'\n";
    out << "&End\n";
    out << "*binFlag          x              Q2           Sigma";
    if (useStat)  out << "            stat";
    if (useUncor) out << "           uncor";
    for (auto su : sStr) out << setfill(' ') << setw(16) << su;
    out << endl;
    int ixsec=0;
    double statCut=50.;  //Ignore bins w/ sat uncertainty above this percentage
    for (int ix=0; ix!=x.size(); ++ix) {        
        if (stat[ix]*100.>statCut) binF[ix]=0;
        out << "  " << binF[ix];
        out << setfill(' ') << setw(16) << x[ix];
        out << setfill(' ') << setw(16) << Q2[ix];
        if (xsecWriteAll || binF[ix]==1) {
            out << setfill(' ') << setw(16) << xsec[ixsec];
            ++ixsec;
        } else out << "        0.000000";
        //Turn err into %
        if (useStat)  out << setfill(' ') << setw(16) << 100.*stat[ix];
        if (useUncor) out << setfill(' ') << setw(16) << 100.*uncor[ix];
        for (int isu=0; isu!=sStr.size();++isu) {  //Systematic uncertainties
            out << setfill(' ') << setw(16) << 100.*syst[isu][ix];
        }
        out << endl; 
    }
    out.close();
}// END xFtableWriter

//AUX function to read binned_events files under ../results
//The results are saved by appending into the vecs passed by reference
void readBinnedEvents(string infile, bool useStat, bool useUncor,
                      vector<int>&    binF1,
                      vector<double>& xlo1,
                      vector<double>& xhi1,
                      vector<double>& xav1,
                      vector<double>& Q2lo1,
                      vector<double>& Q2hi1,
                      vector<double>& Q2av1,
                      vector<double>& sigma1,
                      vector<double>& stat1,
                      vector<double>& uncor1)
{
    //Read the predictions / pseudodata template
    cout << "Processing " << infile << endl;
    ifstream in;
    inStreamOpen(in,infile);
    vector<string> lines;
    string line;
    int step=0;
    while (getline(in,line)) {
        if (step>1) lines.push_back(line);  //Remove 2 lines of header
        else ++step;
    }
    for (auto l : lines) {
        vector<string> lsplit = stringSplit(l);
        binF1.push_back(1);
        xlo1.push_back(str2d(lsplit[0]));
        xhi1.push_back(str2d(lsplit[1]));
        string xavstr = lsplit[2];
        if (xavstr=="nan") xav1.push_back(0.5*(xlo1.back()+xhi1.back()));
        else xav1.push_back(str2d(xavstr));
        Q2lo1.push_back(str2d(lsplit[3]));
        Q2hi1.push_back(str2d(lsplit[4]));
        string Q2avstr = lsplit[5];
        if (Q2avstr=="nan") Q2av1.push_back(0.5*(Q2lo1.back()+Q2hi1.back()));
        else Q2av1.push_back(str2d(Q2avstr));
        string sigstr = lsplit[9];
        if (sigstr=="nan") {
            sigma1.push_back(0.);
            binF1.back() = 0;
        } else sigma1.push_back(pow(10.,str2d(sigstr)));  //Convert log10(sigma)
        string Nevtstr = lsplit[10];
        if (useStat) {
            if (Nevtstr=="0.00000") {
                stat1.push_back(0.);
                binF1.back() = 0;
            }
            else stat1.push_back(str2d(lsplit[11])/str2d(Nevtstr));
        }
        if (useUncor) {
            if (Nevtstr=="0.00000") {
                uncor1.push_back(0.);
                binF1.back() = 0;
            } else uncor1.push_back(str2d(lsplit[12])/str2d(Nevtstr));
        } 
    } //for l : lines
    in.close();
} // END readBinnedEvents

// AUX function to read th_orig from preliminary xFitter run(s)
vector<double> readPrelXsec(string prelname) {
    ifstream in;
    inStreamOpen(in,prelname);
    vector<string> lines;
    int step=0;
    string line;
    while (getline(in,line)) {
        if (step>7) lines.push_back(line);  //Remove 8 lines of header
        else ++step;
    }
    vector<double> xsec1;
    for (auto l : lines) {
        vector<string> lsplit = stringSplit(l);
        xsec1.push_back(str2d(lsplit[6]));
    }
    in.close();
    return xsec1;
}
    
// BEGIN "main"
int main() {

    //True:  write tables for preliminary xFitter run
    //False: write final tables to be used as pseudodata in fits.
    bool writePrel = false;

    string PDF="PDF4LHC21";
    string suffix = ".txt";
    vector<string> expnames = {"FASERv","FASERv2"};
    vector<string> nuIDs = {"14","-14","+-14"};
    string mdir = "datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    
    //For generating univariate Gaussian random numbers
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,1.0);
    
    //Produce subdirectories for tables
    dirCheck(mdir);
    dirCheck(mdir+"/"+PDF+"/prel/");
    //TODO to be most useful, should also untar files copied from theory repo
    //dirCheck(mdir+"grids/");
    //copy_tree("../theory/", mdir+"grids/")  //Pseudocode

    int iexp=136;
    for (auto expname : expnames) {
        
        for (auto nuID : nuIDs) {
                
            iexp += 1;
            //Flags for including uncertainties
            bool useStat  = true;
            bool useUncor = true;
            //Syst. unc are included if a non-empty array is provided
        
            //Init
            string expID = expname + "_" + nuID;
            string infile;
            vector<int> binF;
            vector<double> xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,stat,uncor;
            vector<int> binFp,binFm;
            vector<double> xlom,xhim,xavm,Q2lom,Q2him,Q2avm,sigmam,statm,uncorm;
            if (nuID.find("+-")!=string::npos) {
                //Read nu
                infile = "../results/binned_sysevents_"
                       + replace(expID,"+-","") + suffix;
                binFp.clear();
                readBinnedEvents(infile,useStat,useUncor,
                                 binFp,xlo,xhi,xav,Q2lo,Q2hi,Q2av,
                                 sigma,stat,uncor);
                //Read nub
                infile = "../results/binned_sysevents_"
                       + replace(expID,"+-","-") + suffix;
                binFm.clear();
                xlom.clear();
                xhim.clear();
                xavm.clear();
                Q2lom.clear();
                Q2him.clear();
                Q2avm.clear();
                sigmam.clear();
                statm.clear();
                uncorm.clear();
                readBinnedEvents(infile,useStat,useUncor,
                                 binFm,xlom,xhim,xavm,Q2lom,Q2him,Q2avm,
                                 sigmam,statm,uncorm);
                //Form data vecs for combined nu+nub case
                for (int i=0; i!=xlo.size(); ++i) {                    
                    binF.push_back(binFp[i]*binFm[i]);
                    sigma[i] += sigmam[i];
                    stat[ i] = sqrt( stat[i]*stat[i]  +  statm[i]*statm[i]);
                    uncor[i] = sqrt(uncor[i]*uncor[i] + uncorm[i]*uncorm[i]);
                }
            } else {
                infile = "../results/binned_sysevents_" + expID + suffix;
                readBinnedEvents(infile,useStat,useUncor,
                                 binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,
                                 sigma,stat,uncor);
            }
    
            vector<double> dvtmp;
            vector<string> svtmp;
            vector< vector<double> > systtmp;
    
            //Write xFitter tables for prel run
            if (writePrel) xFtableWriter(mdir,PDF+"/prel/",expname,nuID,
                                         iexp,binF,xav,Q2av,sigma,
                                         stat,uncor,svtmp,systtmp);
            else { //write xFitter tables for final run
                
         
                //Read xFitter output from preliminary run
                vector<double> xsec;
                if (expID.find("+-")!=string::npos) {
                    vector<double> xsecp = readPrelXsec("PDF_profiling/"
                                                        +PDF+"/prel/"
                                                        +replace(expID,"+-","")
                                                        +"/output/"
                                                        +"fittedresults.txt"
                                                        +"_set_0000");
                    vector<double> xsecm = readPrelXsec("PDF_profiling/"
                                                        +PDF+"/prel/"
                                                        +replace(expID,"+-","m")
                                                        +"/output/"
                                                        +"fittedresults.txt"
                                                        +"_set_0000");
                    //Only take xsec at indices where both files have binFlag 1
                    int ip=0;
                    int im=0;
                    for (int ibF=0; ibF!=binF.size(); ++ibF) {                        
                        if ((ip >= xsecp.size()) || (im > xsecm.size())) break;
                        if (binFp[ibF]==1 && binFm[ibF]==1) {
                            xsec.push_back(xsecp[ip]+xsecm[im]);
                            ++ip;
                            ++im;
                        } else if (binFp[ibF]) ++ip;
                        else if   (binFm[ibF]) ++im;
                    }
                } else xsec=readPrelXsec("PDF_profiling/"+PDF+"/prel/"
                                         +replace(expID,"-","m")
                                         +"/output/fittedresults.txt_set_0000");
         
                //Form pseudodata: vary results by the estimated exp. unc.
                vector<double> xsv;         //Cross sections varied 
                vector<double> xsvuncor05;  //..w/in different combinations
                vector<double> xsvstat;     //...of uncertainties
                vector<double> uncor05;     //Optimistic case: halved uncor unc.
                for (double u : uncor) uncor05.push_back(0.5*u);
                for (int ixs=0; ixs!=xsec.size(); ++ixs) {                    
                    double xs = xsec[ixs];
                    double err=0.;
                    double erruncor05=0.;
                    double errstat=0.;
                    if (useStat) {
                        double stat2 = stat[ixs]*stat[ixs];
                        err        += stat2;
                        erruncor05 += stat2;
                        errstat    += stat2;
                    }
                    if (useUncor) {
                        err        +=   uncor[ixs]*uncor[ixs];
                        erruncor05 += uncor05[ixs]*uncor05[ixs];
                    }
                    err        = sqrt(err);
                    erruncor05 = sqrt(erruncor05);
                    errstat    = sqrt(errstat);
                    double rndm = distribution(generator);
                    xsv.push_back(       xs*(1.+       err*rndm));
                    xsvuncor05.push_back(xs*(1.+erruncor05*rndm));
                    xsvstat.push_back(   xs*(1.+   errstat*rndm));
                }
        
                //Write tables for final runs
                xFtableWriter(mdir,PDF+"/",
                              expname,nuID,iexp,
                              binF,xav,Q2av,xsv,
                              stat,uncor,
                              svtmp,systtmp);
                xFtableWriter(mdir,PDF+"/noVar/",
                              expname,nuID,iexp, 
                              binF,xav,Q2av,xsec,
                              stat,uncor,
                              svtmp,systtmp);
                xFtableWriter(mdir,PDF+"/uncorHalf/",
                              expname,nuID,iexp, 
                              binF,xav,Q2av,xsvuncor05,
                              stat,uncor05,
                              svtmp,systtmp);
                xFtableWriter(mdir,PDF+"/statOnly/",
                              expname,nuID,iexp, 
                              binF,xav,Q2av,xsvstat,
                              stat,dvtmp,
                              svtmp,systtmp);
    
            } //else write xFitter tables for final run
        } //for nuID
    } //for expname

    if (writePrel) {
        cout << "Wrote tables for preliminary xFitter runs" << endl;
        cout << "These are used for obtaining predictions convoluted w/ "
             << PDF << endl;
        cout << "which then enter the final fit/profilin stage. Now run the\n"
             << "preliminary fits, then rerun this code w/ writePrel = false"
             << endl;
    } else cout << "Wrote tables for final xFitter runs" << endl;

    return 1;
} // END main