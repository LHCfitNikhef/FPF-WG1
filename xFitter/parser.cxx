//  g++ --std=gnu++0x -o parserEXE parser.cxx
//  ./parserEXE
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <math.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <vector>
#include <utility>

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
    sstream << setprecision(15) << str;
    double d=0;
    if (str!="nan") sstream >> d;
    return d;
}

//AUX function to compare bools 
bool dAgree(double a, double b, double tol=0.001) {
    return (a*(1.+tol) > b) && (a*(1.-tol) < b); 
}

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

/* AUX function to read binned_events files under ../results
 * The results are saved by appending into the vecs passed by reference
 */
vector< vector<double> > readBinnedEvents(string infile, bool useStat, bool useSyst) {

    vector<double> binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,N,stat,tmp;  //Init
    vector< vector<double> > syst;
    cout << "Processing " << infile;
    ifstream in;
    inStreamOpen(in,infile);
    vector<string> lines;
    string line;
    int step=0;
    while (getline(in,line)) {
        if (step>1) lines.push_back(line);  //Omit 2 lines of header
        else ++step;
    }
    for (auto l : lines) {
        
        vector<string> lsplit = stringSplit(l);
        lsplit.pop_back();  //Discard #MC samples at line end
        
        //Read bin limits & flags
        binF.push_back(1.);
        xlo.push_back(str2d(lsplit[0]));
        xhi.push_back(str2d(lsplit[1]));
        xav.push_back(str2d(lsplit[2]));
        Q2lo.push_back(str2d(lsplit[3]));
        Q2hi.push_back(str2d(lsplit[4]));
        Q2av.push_back(str2d(lsplit[5]));

        //Read x-sec
        string sigstr = lsplit[9];
        if (sigstr=="nan") {
            sigma.push_back(0.);
            binF.back() = 0.;
        } else sigma.push_back(str2d(sigstr));
        //} else sigma.push_back(pow(10.,str2d(sigstr)));  //ALT^ convert log10(sigma)

        //Read number of events N and stat unc's effect on N, find stat unc in %
        string Nevtstr = lsplit[10];
        N.push_back(str2d(Nevtstr));
        if (useStat && lsplit.size()>11) {
            if (Nevtstr.find("0.00000")!=string::npos) {
                stat.push_back(0.);
                binF.back() = 0.;
            }
            else stat.push_back(100.0*str2d(lsplit[11])/str2d(Nevtstr)); //Turn into %
        }
        
        //Init matrix of syst unc sources
        vector<double> empty;
        int syststart = 13;  //Read syst unc percentages starting on this column
        if (useSyst && syst.size()==0) {
            for (int i=syststart; i<lsplit.size(); ++i) syst.push_back(empty);
        }
        for (int i=0; useSyst && i+syststart<lsplit.size(); ++i) {
            if (Nevtstr=="0.00000") {
                binF.back() = 0;
                syst[i].push_back(0.);
            } else syst[i].push_back(100.0*str2d(lsplit[i+syststart])); //Turn into %
        }
        
    }  //for l : lines
    in.close();
    cout << "...done" << endl;

    vector< vector<double> > ret = {binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,N};
    if (useStat) ret.push_back(stat);
    for (auto s : syst) ret.push_back(s); //Syst uncs
    
    return ret;    
    
} // END readBinnedEvents

//AUX function to avoid boiler plate
string exptagConstructor(string expID, string origin, bool sub) {
    string lep = sub ? "e" : "\\mu";
    string exptag = replace(expID, "_sub", "");
    exptag        = replace(expID, origin, replace(origin,"_"," "));
    exptag        = replace(exptag,"Rv", "R$\\nu$");
    exptag        = replace(exptag,"_nochargediscrimination",
                                   " $\\nu_{"+lep+"}+\\bar{\\nu}_{"+lep+"}$");
    exptag        = replace(exptag,"_nub", " $\\bar{\\nu}_{"+lep+"}$");
    exptag        = replace(exptag,"_nu",  " $\\nu_{"+lep+"}$");
    return exptag;
}

/* AUX function to write an xFitter table
 * Param  mdir       Master directory above datatables and grids
 *        sdir       Subdirectory under mdir where to write these tables
 *        expname    Experiment name
 *        nuID       "nu" for neutrinos, "nub" for antineutrinos
 *        origin     E.g. "" or "_charm"
 *        sub        True if based on results/*sub i.e. electrons, else muons
 *        index      Dataset must be assigned a unique index
 *        binF_in    Flags 1/0 if bin used in fit; "in"=doubles, changed to int
 *        x,Q2       Vec. of x and Q^2 values for each data point
 *        xsec       Vec. of cross-sections
 *        stat       Vec. of stat. uncertainties for each data point (in %)
 *        uncor      Vec. of uncorrelated uncertainties (in %)
 *        systnames  Names of correlated systematic uncertainty sources
 *        syst       Matrix of corr. syst. unc. values (in %)
 * Returns false if bin mismatch found & tables not written, otherwise true.
 */
bool xFtableWriter(string mdir, 
                   string sdir,
                   string expname, 
                   string nuID,
                   string origin,
                   bool sub,
                   int index, 
                   vector<double> binF_in, 
                   vector<double> x, 
                   vector<double> Q2, 
                   vector<double> xsec, 
                   vector<double> stat, 
                   vector<double> uncor, 
                   vector<string> systnames, 
                   vector< vector<double> > syst)
{
    string expID = expname + origin + (sub ? "_sub_" : "_") + nuID;
    bool prel = sdir.find("prel")!=string::npos;  //Check if writing prel table
    
    //LaTeX and ROOT style strings for experiment identification
    
    string exptag = exptagConstructor(expID,origin,sub);
    string expstr = replace(replace(exptag,"$",""),"\\","#");
    vector<int> binF;
    for (double b : binF_in) binF.push_back((int)b);
    
    //Check if uncor, stat unc given
    bool useStat  = stat.size()  > 0;
    bool useUncor = uncor.size() > 0;
                         
    //Check if all xsec given, or only those w/ binflag 1 in which case bin matching is required 
    bool xsecWriteAll = xsec.size() == binF.size();
    if (!xsecWriteAll) {
        int bFsum=0;
        for (int b : binF) bFsum+=b;
        if (bFsum!=xsec.size()) {  //Ensure xsec vec size equals #[expected active bins]
            cout << "ERROR binFlag sum=" << bFsum
                 << " != xsec size=" << xsec.size() << " in " << expID << endl;
            return false;
        }
    }

    dirCheck(mdir+sdir);
        
    fstream out;
    string thexp = "-thexp.dat";
    string outname = mdir+sdir+expID+thexp;
    outStreamOpen(out,outname);
    out << "* " << exptag << " CC DIS pseudodata\n";
    out << "&Data\n";
    out << "   Name = '" << replace(exptag,"_optimistic"," optimistic")
        << "'\n";
    out << "   IndexDataset = " << index << "\n"; //Give set a unique ID
    out << "   Reaction = 'CC nup'\n\n";
    out << "   NData = " << x.size() << "\n";
    int Nerr = systnames.size();
    if (useStat)  ++Nerr;
    if (useUncor) ++Nerr;
    int Ncol = 1 + 2 + 1 + Nerr;  //Flag + bins + sigma + err.s
    out << "   NColumn = " << Ncol << "\n";
    out << "   ColumnType = 'Flag',2*'Bin','Sigma'";
    if (Nerr>0) out << "," << Nerr << "*'Error'";
    out << "\n";
    out << "   ColumnName = 'binFlag','x','Q2','Sigma'";
    if (useStat)  out << ",'stat'";
    if (useUncor) out << ",'uncor'";
    for (auto s : systnames) out << ",'" << s << "'";
    out << "\n";
    out << "   TheoryType = 'expression\'\n";     
    //TODO theory path must account for sub once grids available
    //Theory files for nochargediscrimination case are still divided, sum them
    string chorigin = origin == "_charm" ? origin : "";
    if (nuID.find("nochargediscrimination")!=string::npos) {
        out << "   TermType   = 'reaction','reaction'\n";
        out << "   TermName   = 'P','A'\n";
        out << "   TermSource = 'PineAPPL','PineAPPL'\n";
        out << "   TermInfo   = 'GridName="
            << mdir << "grids/grids-"
            << replace(expID,"_optimistic","")
            << "-a1/nu_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4',\n";
        out << "                'GridName="
            << mdir << "grids/grids-"
            << replace(expID,"_optimistic","")
            << "-a1/nub_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4'\n";
        out << "   TheorExpr  = 'P+A'\n\n";
    } else {
        out << "   TermType   = 'reaction'\n";
        out << "   TermName   = 'P'\n";
        out << "   TermSource = 'PineAPPL'\n";
        out << "   TermInfo   = 'GridName=";
        out << mdir << "grids/grids-"
            << replace(expID,"_optimistic","")
            << "-a1/"+nuID+"_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4'\n";
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
    out << "*binFlag" << setfill(' ') << setw(19) << "x"
                      << setfill(' ') << setw(24) << "Q2"
                      << setfill(' ') << setw(24) << "Sigma";
    if (useStat)  out << setfill(' ') << setw(24) << "stat";
    if (useUncor) out << setfill(' ') << setw(24) << "uncor";
    for (auto su : systnames) out << setfill(' ') << setw(24) << su;
    out << endl;
    int ixsec=0;
    double uncCut=50.;  //In final tables, ignore bins w/ any unc above this %
    for (int ix=0; ix!=x.size(); ++ix) {
        out << "  ";
        int bt = binF[ix];  //N.B. don't modify binF[ix], consistency* required!
        //if (!prel && stat[ix] > uncCut) bt=0;
        //for (int isu=0; isu!=syst.size(); ++isu) if (!prel && syst[isu][ix]>uncCut) bt=0;
        if (xsec[ixsec] < 0) bt=0;  //Don't include bin if x-sec negative
        out << bt;
        out << setfill(' ') << setw(24) << setprecision(15) << x[ix];
        out << setfill(' ') << setw(24) << setprecision(15) << Q2[ix];
        if (xsecWriteAll || binF[ix]==1) {  //what matters here is the original*
            out << setfill(' ') << setw(24) << xsec[ixsec];
            ++ixsec;
        } else out << "        0.000000";
        //Turn err into %
        if (useStat)  out << setfill(' ') << setw(24) << stat[ix];
        if (useUncor) out << setfill(' ') << setw(24) << uncor[ix];
        for (int isu=0; isu!=syst.size(); ++isu) {  //Syst. unc.
            out << setfill(' ') << setw(24) << syst[isu][ix];
        }
        out << endl; 

    } //Loop ix (over x vector size)
    
    out.close();
    cout << "Wrote " << outname << endl;
    return true;
}// END xFtableWriter

/* AUX function to read a previously produced .dat.
 * Returns the table contents in a matrix w/ column index first
 * Also reads/overwrites the column names into colNames
 */
vector< vector<double> > xFtableReader(vector<string>& colNames,
                                       string mdir, 
                                       string sdir,
                                       string expname, 
                                       string nuID, 
                                       string origin,
                                       bool sub)
{
    string line;
    stringstream sstream;
    ifstream in;
    double dtmp;
    vector<double> vtmp;
    vector< vector<double> > Mtmp;
    string expID = expname + origin + (sub ? "_sub_" : "_") + nuID;
    string infile = mdir + sdir + expID + "-thexp.dat";

    //Open xFitter .dat file for reading
    inStreamOpen(in,infile);
    colNames.clear();

    //Read header
    while (getline(in,line) && line.find("*binFlag")==string::npos) {
        //Find list of column names, trim wrapping apostrophes & save
        if (line.find("ColumnName")!=string::npos &&
            line.find("=")!=string::npos            )
        {
            line.erase(0,line.find("=")+1);
            sstream.clear();  sstream.str("");
            sstream << line;
            string stmp;
            while(sstream.good()) {
                getline(sstream,stmp,',');
                if (stmp.find("'")!=string::npos) {
                    stmp.erase(0,stmp.find("'")+1);
                    stmp.erase(stmp.find("'"));
                }
                colNames.push_back(stmp);  //Save via par passed as ref
            }
        }
    }

    //Read data table
    while (getline(in,line)) {
        vtmp.clear();
        sstream.clear();  sstream.str("");
        sstream << line;
        while (sstream >> dtmp) vtmp.push_back(dtmp);
        Mtmp.push_back(vtmp);
    }
    
    in.close();

    //Check that data was obtained
    if (Mtmp.size()==0) {
        cout << "ERROR in xFtableReader processing " << infile;
        return Mtmp;
    }

    //Transpose Mtmp s.t. Mret[i] yields ith column of the read table format
    vector< vector<double> > Mret;
    for (int i=0; i!=Mtmp[0].size(); ++i) {
        vtmp.clear();
        for (int j=0; j!=Mtmp.size(); ++j) vtmp.push_back(Mtmp[j][i]);
        Mret.push_back(vtmp);
    }
    return Mret;
    
} // END xFtableReader

// AUX function to read th_orig from preliminary xFitter run(s)
vector< vector<double> > readPrelXsec(string prelname) {

    vector<double> xbin, Q2bin, xsec;

    //Read output of preliminary xFitter run
    ifstream in;
    inStreamOpen(in,prelname);
    vector<string> lines;
    int step=0;
    string line;

    //Remove header
    while (getline(in,line) && line.find("th orig")==string::npos);

    //Read data
    while (getline(in,line)) lines.push_back(line);
    
    //Store relevant columns for return
    for (auto l : lines) {
        vector<string> lsplit = stringSplit(l);
        xbin.push_back( str2d(lsplit[0]));
        Q2bin.push_back(str2d(lsplit[1]));
        xsec.push_back( str2d(lsplit[6]));
    }
    in.close();
    return {xbin, Q2bin, xsec};
}

/* Compute full correlation matrix based on stat and syst uncertainties
 * Param  N      Total number of events per bin (left for reference)
 *        stat   Vector including statistical uncertainties
 *        systs  Matrix of systematic unc., each entry is a vector corresponding
 *               to a single source of systematic uncertainty
 *        fred   Reduction factor for systematic uncertainties
 * Returns the correlation matrix
 */
vector<vector<double>> correlations(vector<double> N,
                                    vector<double> stat,
                                    vector<vector<double>> systs,
                                    double fred)
{
    //Init
    vector<double> tmp;
    for (int i=0; i!=N.size(); ++i) tmp.push_back(0.0);
    vector< vector<double> > cov;
    for (auto itmp : tmp) cov.push_back(tmp);
    double fred2 = fred*fred;
    
    //Check if stat unc contribution can be added
    bool addStat = stat.size()==N.size();
    if (!addStat) cout << "WARNING can't add stat contribution to corr" << endl;
    
    //Find covariance matrices first, turn later into correlation matrix
    for (int i=0; i!=N.size(); ++i) {
        for (int j=0; j!=N.size(); ++j) {
            //Assume uncorrelated stat unc -- technically leads to diag(1,..,1)
            //for the relative covariance, left here for future extendability
            if (addStat && i==j) cov[i][j]+=stat[i]*stat[j];
            //Fully correlated sst unc
            for (vector<double> syst : systs) {                
                if (syst.size()==N.size()) cov[i][j] += fred2*syst[i]*syst[j];
                else cout << "WARNING correlations: syst size mismatch" << endl;
            }
            //Only relative cov. needed, so #evts and %-conversion would cancel
            //cov[i][j] *= N[i]*N[j];  //Left for reference
        }
    }
    
    //Construct correlation matrix from relative covariance
    vector< vector<double> > corr = cov;
    for (int i=0; i!=N.size(); ++i) {
        for (int j=0; j!=N.size(); ++j) {
            if (!dAgree(cov[i][i],0) && !dAgree(cov[j][j],0)) {
                corr[i][j] /= (sqrt(cov[i][i])*sqrt(cov[j][j]));
            } else corr[i][j] = 0;
        }
    }
        
    return corr;
} // END correlations

/* Write covariance matrix table in xFitter format
 * Param  expname  Experiment name string e.g. "FASERv2","FASERv"...
 *        nuID     Consider neutrinos ("nu") antineutrinos ("nub") or 
 *                 both ("nochargediscrimination")
 *        origin   Write tables from "_inclusive" or "_charm" production data
 *        sub      True if based on results/*sub i.e. electrons, else muons
 *        nameAdd  Add identifier tag to filename, e.g. "onlyEl"
 * Returns true if table wrote successfully, false in case of errors
 */
bool writeCorr(string expname, 
               string nuID, 
               string origin,
               bool sub,
               string nameAdd,
               vector<double> xav,
               vector<double> Q2av,
               vector< vector<double> > corr)
{               
    string expID  = expname + origin + (sub ? "_sub_" : "_") + nuID;
    string exptag = exptagConstructor(expID,origin,sub);

    //Check dimension agreement
    if (xav.size()!=Q2av.size()) {
        cout << "ERROR writeCorr xav & Q2av size mismatch in " << exptag << endl;
        return false;
    }
    if (xav.size()!=corr.size()) {
        cout << "ERROR writeCorr xav & corr size mismatch in " << exptag << endl;
        return false;
    }
    for (int i=0; i!=corr.size(); ++i) {
        if (xav.size()!=corr[i].size()) {
            cout << "ERROR writeCorr xav & corr[" << i << "] size mismatch in "
                 << exptag << endl;
            return false;
        }
    }

    string datadir = "./datafiles/lhc/fpf/neutrinoDIS/pseudodata/";    
    string outname = datadir + expID + nameAdd + ".corr";
    fstream out;
    outStreamOpen(out,outname);
    out << "! Full correlation matrix" << endl;
    out << "&StatCorr"           << endl;
    out << "  Name1 = '"
        << replace(exptag,"_optimistic"," optimistic")
        << "'"   << endl;
    out << "  Name2 = '"
        << replace(exptag,"_optimistic"," optimistic")
        << "'\n" << endl;
    out << "  NIdColumns1 = 2"          << endl;
    out << "  NIdColumns2 = 2\n"        << endl;
    out << "  IdColumns1 = 'x', 'Q2'"   << endl;
    out << "  IdColumns2 = 'x', 'Q2'\n" << endl;
    out << "  NCorr = " << xav.size()*xav.size() << "\n" << endl;
    out << "  MatrixType = 'Full correlation matrix'" << endl;
    out << "&End" << endl;
    for (int i=0; i!=xav.size(); ++i) {
        for (int j=0; j!=xav.size(); ++j) {
            out <<setfill(' ')<<setw(20)<<setprecision(15)<< xav[ i];
            out <<setfill(' ')<<setw(20)<<setprecision(15)<< Q2av[i];
            out <<setfill(' ')<<setw(20)<<setprecision(15)<< xav[ j];
            out <<setfill(' ')<<setw(20)<<setprecision(15)<< Q2av[j];
            out <<setfill(' ')<<setw(25)<<setprecision(15)<< corr[i][j] <<endl;
        }
    }
    out.close();
    cout << "Wrote " << outname << endl;
    
    return true;
} // END writeCorr
    
/* Write xFitter tables for preliminary run
 * Param  expname  Experiment name string e.g. "FASERv2","FASERv"...
 *        nuID     Consider neutrinos ("nu") antineutrinos ("nub") or 
 *                 both ("nochargediscrimination")
 *        origin   Write tables from "_inclusive" or "_charm" production data
 *        sub      Use tables from the results/*_sub directories?
 *        iexp     xFitter datatable index, all datasets used in a single fit
 *                 require unique indices for identification
 *        fred     Corr. tables produced with and without this reduction factor
 * Returns exit status of xFtableWriter: false if bin mismatch found & tables 
 *                                      not written, else true.
 */
bool writeDatPrel(string expname, 
                  string nuID, 
                  string origin,
                  bool sub,
                  int iexp, 
                  bool useStat, 
                  bool useSyst,
                  double fred)
{
    string expID = expname + origin + (sub ? "_sub_" : "_") + nuID;
    string suffix = ".txt";
    string mdir = "datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    dirCheck(mdir);
    string origindir = (origin.find("_charm")!=string::npos ? "CHARM"
                                                            : "INCLUSIVE");
    origindir += (sub ? "_sub/" : "/");
    string binevtbase = "../results/"
                      + origindir
                      + expname
                      + "/clipped_nan/clipped_nan_binned_sysevents_";

    //Make sure pineappl grids exist
    string thpath = "../theory/grids/";
    string gridsub = "grids-" + expID + "-a1";  //Subdir name in git
    string chorigin = "";
    if (origin.find("_charm")!=string::npos) chorigin = origin;
    vector<string> grids;
    //TODO paths may account for mu or e once new grids become available
    grids.push_back("nu_A_1-" +expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4");
    grids.push_back("nub_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4");
    dirCheck(thpath+gridsub);
    string gitln = "https://github.com/juanrojochacon/FPF-WG1/raw/main/";
    bool dloadth = false;
    for (string grid : grids) {
        if (system(("ls "+thpath+gridsub+"/"+grid).c_str())==512) {
            cout << "Theory for " << endl;
            cout << "    " << thpath << gridsub << "/" << grid << endl;
            cout << "    not found, attempting to download" << endl;
            dloadth = true;
        }
    }
    if (dloadth) {
        system(("wget " + gitln + "theory/grids/" + gridsub + ".tar").c_str());
        system(("tar -xf " + gridsub + ".tar").c_str());
        for (string grid : grids) {
            system(("mv grids/"+grid+" " + thpath+gridsub+"/"+grid).c_str());
        }
        system("rm -rf grids");
        system(("rm " + gridsub + ".tar").c_str());
    }
    //Put links to main theory/grids subdirs under xfitter datafiles
    string thpath4ln =  "../../../../../../theory/grids/";
    if (system(("ls " + mdir + "grids/").c_str())==512) {
        system(("ln -s " + thpath4ln + " " + mdir + "grids").c_str());
    }
    
    //Init
    string infile;
    vector<double> binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,N,stat;
    vector<double> binFp,binFm,xlom,xhim,xavm,Q2lom,Q2him,Q2avm,sigmam,Nm,statm;

    //Include unc.s already in prel tables if available, handy to read all input
    //info from prel table in the next stage, no need to revisit binned_events
    vector<string> systnames;
    vector< vector<double> > syst;

    //Read binned_events, no more need to sum tables for nochargediscrimination
    infile = binevtbase + expID + suffix;
    vector< vector<double> > BE = readBinnedEvents(infile,useStat,useSyst);
    binF  = BE[0];
    xlo   = BE[1];
    xhi   = BE[2];
    xav   = BE[3];
    Q2lo  = BE[4];
    Q2hi  = BE[5];
    Q2av  = BE[6];
    sigma = BE[7];
    N     = BE[8];
    if (useStat) stat = BE[9];
    for (int i=(useStat ? 10 : 9); i<BE.size(); ++i) syst.push_back(BE[i]);
    
    //Compute unc as quadratic sum of all systs here
    //REDUNDANT since syst unc must be called uncor anyway (see below) in 
    //current approach, including a non-empty uncor would just account for the 
    //same sources twice, but this code is left for reference.
    vector<double> uncor;
    //for (int row=0; useSyst && row!=syst[0].size(); ++row) {
    //    double dtmp=0;
    //    for (int col=0; col!=syst.size(); ++col) {
    //        dtmp += syst[col][row]*syst[col][row];
    //    }
    //    uncor.push_back(sqrt(dtmp));
    //}
    
    //Construct generic syst unc names
    stringstream sstream;
    for (int i=0; i!=syst.size(); ++i) {
         sstream.clear();  sstream.str("");
         //When using full correlation matrix or systematic correlations, 
         //systematic uncertainty columns must be called 'uncor const'
         sstream << "uncor const";
         systnames.push_back(sstream.str());
     }

    //Correlations with all syst unc, no reduction factor
    vector<vector<double>> corr = correlations(N,stat,syst,1.0);
    if (!writeCorr(expname,nuID,origin,sub,"",xav,Q2av,corr)) return false;
    //With reduction factor
    vector<vector<double>> corrfred = correlations(N,stat,syst,fred);
    if (!writeCorr(expname,nuID,origin,sub,
                   "_systVar05",
                   xav,Q2av,corrfred)) return false;

    //Correlations with only lepton energy in syst unc, no reduction factor
    vector< vector<double> > systEl;
    systEl.push_back(syst[1]);  //El is 2nd syst unc in binned_events files
    vector<vector<double>> corrEl = correlations(N,stat,systEl,1.0);
    if (!writeCorr(expname,nuID,origin,sub,
                   "_onlyEl",
                   xav,Q2av,corrEl)) return false;
    //With reduction factor
    vector<vector<double>> corrElfred = correlations(N,stat,systEl,fred);
    if (!writeCorr(expname,nuID,origin,sub,
                   "_systVar05_onlyEl",
                   xav,Q2av,corrElfred)) return false;
    
    //Write datatable for preliminary xFitter run
    return xFtableWriter(mdir,"prel/",expname,nuID,origin,sub,iexp,
                         binF,xav,Q2av,sigma,stat,uncor,systnames,syst);
    
} // END writeDatPrel

bool writeDatFinal(string PDF, 
                   string expname, 
                   string nuID, 
                   string origin,
                   bool sub,
                   int iexp,
                   double fred)
{
    string mdir = "datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    
    //For generating univariate Gaussian random numbers
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,1.0);
    
    //Init
    string expID = expname + origin + (sub ? "_sub_" : "_") + nuID;
    string infile;
    vector<double> dvtmp;
    vector<string> svtmp, systnames, systnamesEl;
    vector< vector<double> > systtmp, systunc, systuncEl;
    double fred2 = fred*fred;

    //Read xsec values from previous xFitter run
    vector< vector<double> > prelXS;  //To contain xbins, Q2bins, th orig
    string prelname = "PDF_profiling/"+PDF+"/prel/"
                    + replace(replace(expID,"-","m"),"+","p")
                    + "/output/fittedresults.txt_set_0000";
    prelXS = readPrelXsec(prelname);
    if (prelXS.size()!=3) {
        cout << "ERROR prelXS size!=3 for " << prelname << endl;
        return false;
    } else if (prelXS[0].size()==0) {
        cout << "ERROR could not read preliminary x-sec from: " << endl;
        cout << "  " << prelname << endl;
        return false;
    }
    vector<double> xb=prelXS[0],Q2b=prelXS[1],xsec=prelXS[2];

    //Read input data to the previous run^, required for 
    //matching tables as xFitter output skips unused bins
    vector<string> cols;
    vector< vector<double> > Mdat = xFtableReader(cols, mdir, "prel/",
                                                  expname,nuID,origin,sub);
    if (Mdat.size()==0) return false;
    map<string,int> colMap;
    for (int i=0; i!=cols.size(); ++i) colMap[cols[i]] = i;
    vector<double> binF = Mdat[colMap["binFlag"]];
    vector<double> xav  = Mdat[colMap["x"]];
    vector<double> Q2av = Mdat[colMap["Q2"]];
    
    //Check if stat unc given
    vector<double> stat;
    bool useStat = find(cols.begin(),cols.end(),"stat")!=cols.end();
    if (useStat) stat = Mdat[colMap["stat"]];

    //Check if uncor given (expect quadratic sum of all syst unc)
    vector<double> uncor;
    bool useUncor = find(cols.begin(),cols.end(),"uncor")!=cols.end();
    if (useUncor) uncor = Mdat[colMap["uncor"]];
    
    //Assume syst unc are listed after stat and uncor, if those are given
    int isyst = colMap["Sigma"] + 1 + (int)useStat + (int)useUncor;
    for (int i=isyst; i<cols.size(); ++i) {
        systnames.push_back(cols[i]);
        systunc.push_back(Mdat[i]);
        //Assume lepton E is 2nd listed unc in binned events
        if (i==isyst+1) {
            systnamesEl.push_back(cols[i]);
            systuncEl.push_back(Mdat[i]);
        }
    }
    bool useSyst = systnames.size() > 0;

    //Bin-match xsec vector
    int i=0;
    vector<double> xsec1;
    for (int j=0; j!=Mdat[0].size(); ++j) {
        if (   ((int)Mdat[0][j])==1
            && i<xb.size()
            && dAgree(xav[j],  xb[i])
            && dAgree(Q2av[j], Q2b[i]))
        {
            xsec1.push_back(xsec[i]);
            ++i;
        } else {
            binF[j]=0;
            xsec1.push_back(0.0);
        }
    }
    xsec.clear();
    xsec=xsec1;
    
    //Form pseudodata: vary cross sections by the estimated exp unc
    //and different combinations of uncertainties. In forming the
    //pseudodata, correlated uncertainties are expected to be more
    //constraining than uncorrelated, so the syst unc are taken as
    //uncorrelated in the pseudodata variations even if correlations
    //are used for the theory computation in the fit. To estimate 
    //how much better the correlated data would constrain the PDFs, 
    //a reduction factor fred (e.g. 0.5) is applied to the correlated case
    vector<double> xsv, xsvstat, xsvuncfred, xsvuncElfred;
    for (int ixs=0; ixs!=xsec.size(); ++ixs) {                    

        //Generate N = (#syst.unc.s + #stat.unc.) random numbers
        vector<double> rndm;
        for (int i=0; i!=1+systunc.size(); ++i) {
            rndm.push_back(distribution(generator));
        }

        //Two options to form pseudodata variations:
        //(1) varied = original*(1 + r*delta),
        //    where: r = univariate gaussian random number,
        //           delta=sqrt(stat^2 + sum_i syst_i^2)
        //(2) varied = original*(1 + r*delta_stat + sum_i r_i*delta_syst_i ),
        //    where: r¸ r_i all different univ. gauss. random numbers,
        //           delta_stat   =      stat unc (no quadratic sum!)
        //           delta_syst_i = i:th syst unc (no quadratic sum!)
        
        //For obtaining quadratic sums of all stat and syst uncs
        double stat2 = (useStat ? pow(stat[ixs],2) : 0.);
        double syst2=0;
        for (int i=0; i!=systunc.size(); ++i) syst2 += pow(systunc[i][ixs],2);
        //For possibility to consider E_lepton as only syst unc src:
        double systEl2=0;
        if (systuncEl.size()!=0) systEl2 = pow(systuncEl[0][ixs],2);

        //Option 1:
        //[Random gaussian]*[total uncertainty from quad. sum of stat & syst]:
        //Factors of 0.01 due to Div by 100 since unc.s given in %
        double errstat         = 0.01*rndm[0]*sqrt(stat2);
        double errstatsyst     = 0.01*rndm[0]*sqrt(stat2 + syst2);
        double errstatsystfred   = 0.01*rndm[0]*sqrt(stat2 + syst2*fred2);
        double errstatsystElfred = 0.01*rndm[0]*sqrt(stat2 + systEl2*fred2);
        
        //Option 2: fetch sums of various combinations of stat & syst uncs
        double errsyst = 0.0;
        double errsystfred = 0.0;
        for (int i=0; i!=systunc.size(); ++i) {
            //Alternatively, individual factors for separate uncs, no quad sum:
            //Factors of 0.01 due to Div by 100 since unc.s given in %
            //Already account for rndm to assign them individually per source
            errsyst     += 0.01*systunc[i][ixs]*rndm[i+1];
            errsystfred += 0.01*systunc[i][ixs]*rndm[i+1]*fred;
        }
        double errsystfredEl = 0.01*sqrt(systEl2)*rndm[useStat ? 2 : 1]*fred;
        
        //Vary x-sec to obtain pseuodata 
        double xs = xsec[ixs];
        xsvstat.push_back(xs*(1.0 + errstat));  //Stat only, same for 1 & 2 opts
        //Option (1): only a single random number, quadratic sums
        bool useQuadSum = true;
        if (useQuadSum) {
            xsv.push_back(         xs*(1.0 +     errstatsyst));
            xsvuncfred.push_back(  xs*(1.0 +   errstatsystfred));
            xsvuncElfred.push_back(xs*(1.0 + errstatsystElfred));
        } else {  //Opt. 2: separate random numbers for all sources, no quad sum
            xsv.push_back(         xs*(1.0 + errstat + errsyst      ));
            xsvuncfred.push_back(  xs*(1.0 + errstat + errsystfred  ));
            xsvuncElfred.push_back(xs*(1.0 + errstat + errsystfredEl));
        }
        
    }
    
    /* Write tables for final runs */

    string sdir=""; //Init subdirectory name (of the form "PDF/uncLevel")
    
    //Including all syst unc
    sdir = PDF+"/syst/";
    if (!xFtableWriter(mdir,sdir,
                       expname,nuID,origin,sub,iexp,
                       binF,xav,Q2av,xsv,                 //3 bins, x-sec
                       stat,dvtmp,                        //stat, uncor
                       systnames,systunc)) return false;  //systnames, systmatrix                       
                       
    //Use all syst unc but decrease them in pseudodata variation by factor fred
    sdir = PDF+"/systVar05/";
    if (!xFtableWriter(mdir,sdir,
                       expname,nuID,origin,sub,iexp,
                       binF,xav,Q2av,xsvuncfred,          //3 bins, x-sec
                       stat,dvtmp,                        //stat, uncor
                       systnames,systunc)) return false;  //systnames, systmatrix

    //Include only lepton energy uncertainty in syst
    sdir = PDF+"/systVar05El/";
    if (!xFtableWriter(mdir,sdir,
                       expname,nuID,origin,sub,iexp,
                       binF,xav,Q2av,xsvuncElfred,            //3 bins, x-sec
                       stat,dvtmp,                            //stat, uncor
                       systnamesEl,systuncEl)) return false;  //systnames, systmatrix
                       
    //Stat unc only, no syst; neither as unc or in pseudodata variation
    if (!xFtableWriter(mdir,PDF+"/statOnly/",
                       expname,nuID,origin,sub,iexp, 
                       binF,xav,Q2av,xsvstat,         //3 bins, x-sec
                       stat,dvtmp,                    //stat, uncor
                       svtmp,systtmp)) return false;  //systnames, systmatrix

    //Syst unc considered fully uncorrelated
    if (!xFtableWriter(mdir,PDF+"/uncor/",
                       expname,nuID,origin,sub,iexp, 
                       binF,xav,Q2av,xsv,             //3 bins, x-sec
                       stat,uncor,                    //stat, uncor
                       svtmp,systtmp)) return false;  //systnames, systmatrix

    ////No variation in pseudodata, for checking chi2=0 case
    //if (!xFtableWriter(mdir,PDF+"/noVar/",
    //                   expname,nuID,origin,sub,iexp, 
    //                   binF,xav,Q2av,xsec,            //3 bins, x-sec
    //                   stat,uncor,                    //stat, uncor
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix
    
    ////Stat unc only, as uncor
    //if (!xFtableWriter(mdir,PDF+"/statAsUncor/",
    //                   expname,nuID,origin,sub,iexp, 
    //                   binF,xav,Q2av,xsvstat,         //3 bins, x-sec
    //                   dvtmp,stat,                    //stat, uncor; flip intended
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix

    ////Stat unc as uncor, syst[0] as stat
    //if (!xFtableWriter(mdir,PDF+"/systAsStat/",
    //                   expname,nuID,origin,sub,iexp, 
    //                   binF,xav,Q2av,xsvstat,         //3 bins, x-sec
    //                   systunc[0],stat,               //stat, uncor; flip intended
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix

    return true;
} // END writeDatFinal

int main() {
    
    //BEGIN user input
    vector<string> PDFs = {"PDF4LHC21","EPPS21nlo_CT18Anlo_W184"};
    vector<string> expnames = {"FASERv2","FASERv"}; //,"AdvSND","FLArE10","FLArE100"}; //,"SND"
    vector<string> nuIDs = {"nu","nub","nochargediscrimination"};
    vector<string> origins = {"_inclusive","_charm"};

    //True:  write tables for preliminary xFitter run
    //False: write final tables to be used as pseudodata in fits.
    bool prel = false;
    
    //Reduction factor for systematic uncertainties, e.g. 0.5 for systVar05
    //Tables will be produced without this i.e. at fred=1.0, and with the fred
    //set here.
    double fred = 0.5;

    //Flags for including uncertainties:
    //Turn off manually if binned_events don't contain these or not to be used
    bool useStat = true;  //Do binned_events tables contain stat unc?
    bool useSyst = true;  //              -||-              syst unc?
    //END user input

    //Write tables
    int iexp=137;
    for (string expname : expnames) {
        for (string nuID : nuIDs) {
            for (string origin : origins) {
                for (bool sub : {false}) { //TODO sub can only be true once sub th tables available
                    if (prel) {
                        if (!writeDatPrel(expname,nuID,origin,sub,iexp,useStat,useSyst,fred)) return -1;
                    } else {
                        for (string PDF : PDFs) {
                            if (!writeDatFinal(PDF,expname,nuID,origin,sub,iexp,fred)) return -1;
                        }
                    }
                    ++iexp;
                }
            }
        }
    }        
    
    //Print instructions before exiting
    if (prel) {
        cout << "Wrote tables for preliminary xFitter runs" << endl;
        cout << "used for obtaining predictions convoluted w/ " << endl;
        for (string PDF : PDFs) cout << "  " << PDF << endl;
        cout << "which then enter the final fit/profiling stage. Now run the\n"
             << "preliminary fits, then rerun this code w/ prel = false"
             << endl;
    } else cout << "Wrote tables for final xFitter runs" << endl;

    return 1;
}
