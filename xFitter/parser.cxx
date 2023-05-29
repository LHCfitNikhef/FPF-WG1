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

//AUX function to read binned_events files under ../results
//The results are saved by appending into the vecs passed by reference
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
        if (step>1) lines.push_back(line);  //Remove 2 lines of header
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
        string Nevtstr = lsplit[10];
        N.push_back(str2d(Nevtstr));

        //Read stat unc
        if (useStat && lsplit.size()>11) {
            if (Nevtstr.find("0.00000")!=string::npos) {
                stat.push_back(0.);
                binF.back() = 0.;
            }
            else stat.push_back(str2d(lsplit[11])/str2d(Nevtstr));
        }

        //Init matrix of syst unc sources
        vector<double> empty;
        if (useSyst && syst.size()==0) for (int i=12; i<lsplit.size(); ++i) syst.push_back(empty);
        //Read syst unc
        for (int i=0; useSyst && i+12<lsplit.size(); ++i) {
            if (Nevtstr=="0.00000") {
                binF.back() = 0;
                syst[i].push_back(0.);
            } else syst[i].push_back(str2d(lsplit[i+12]));
        }
        
    }  //for l : lines
    in.close();
    cout << "...done" << endl;

    vector< vector<double> > ret = {binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,N,stat};
    for (auto s : syst) ret.push_back(s); //Syst uncs
    
    return ret;    
    
} // END readBinnedEvents

//AUX function to write an xFitter table
//Param  mdir       Master directory above datatables and grids
//       sdir       Subdirectory under mdir where to write these tables
//       expname    Experiment name
//       nuID       "nu" for neutrinos, "nub" for antineutrinos
//       origin     E.g. "" or "_charm"
//       index      Dataset must be assigned a unique index
//       binF_in    Flags 1/0 if bin used in fit; "in"=doubles, changed to int
//       x,Q2       Vec. of x and Q^2 values for each data point
//       xsec       Vec. of cross-sections
//       stat       Vec. of statistical uncertainties for each data point (in %)
//       uncor      Vec. of uncorrelated uncertainties (in %)
//       systnames  Names of correlated systematic uncertainty sources
//       syst       Matrix of corr. syst. unc. values (in %)
//Returns false if bin mismatch found & tables not written, otherwise true.
bool xFtableWriter(string mdir, string sdir,
                   string expname, 
                   string nuID,
                   string origin,
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
    string expID = expname + origin + "_" + nuID;

    //LaTeX and ROOT style strings for experiment identification
    string exptag = replace(expID, origin, replace(origin,"_"," "));
    exptag        = replace(exptag,"Rv", "R$\\nu$");
    exptag        = replace(exptag,"_nochargediscrimination",
                                   " $\\nu_{\\mu}+\\bar{\\nu}_{\\mu}$");
    exptag        = replace(exptag,"_nub", " $\\bar{\\nu}_{\\mu}$");
    exptag        = replace(exptag,"_nu",  " $\\nu_{\\mu}$");
    string expstr = replace(exptag,"$","");
    expstr        = replace(expstr,"\\","#");    
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
    double uncCut=50.;  //Ignore bins w/ any unc above this %
    for (int ix=0; ix!=x.size(); ++ix) {
        out << "  ";
        int bt = binF[ix];  //N.B. don't modify binF[ix], consistency* required!
        if      (stat[ix] > uncCut) bt=0;
        else if (xsec[ixsec] < 0  ) bt=0;         
        for (int isu=0; isu!=syst.size(); ++isu) if (syst[isu][ix]>uncCut) bt=0;
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

//AUX function to read a previously produced .dat.
//Returns the table contents in a matrix w/ column index first
//Also reads/overwrites the column names into colNames
vector< vector<double> > xFtableReader(vector<string>& colNames,
                                       string mdir, 
                                       string sdir,
                                       string expname, 
                                       string nuID, 
                                       string origin)
{
    string line;
    stringstream sstream;
    ifstream in;
    double dtmp;
    vector<double> vtmp;
    vector< vector<double> > Mtmp;
    string expID = expname + origin + "_" + nuID;
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
    
//Write xFitter tables for preliminary run
//Returns exit status of xFtableWriter: false if bin mismatch found & tables not written, else true.
bool writeDatPrel(string PDF, 
                  string expname, 
                  string nuID, 
                  string origin,
                  int iexp, 
                  bool useStat, 
                  bool useSyst)
{
    string expID = expname + origin + "_" + nuID;
    string suffix = ".txt";
    string mdir = "datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    string origindir = (origin=="_charm" ? "CHARM/" : "INCLUSIVE/");
    string binevtbase = "../results/"
                      + origindir
                      + expname
                      + "/clipped_nan/clipped_nan_binned_sysevents_";

    //Links to grids
    string gridsub = "grids-" + expID + "-a1";
    dirCheck(mdir+"grids/");  //Also check if grids exist / to be downloaded / ln -s
    dirCheck(mdir+"grids/"+gridsub);
    string thpath = "../theory/grids/";
    if (system(("ls "+thpath+gridsub).c_str())==512) {
        system(("tar -xf "+thpath+gridsub+".tar").c_str());
        system(("mv grids "+thpath+gridsub).c_str());
    }
    string thpath4ln =  "../../../../../../../../theory/grids/" + gridsub;
    string chorigin = origin == "_charm" ? origin : ""; 
    for (string gd : { "nu_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4",
                      "nub_A_1-"+expname+origin+"-XSFPFCC"+chorigin+".pineappl.lz4"})
    {
        if (system(("ls "+mdir+"grids/"+gridsub+"/"+gd).c_str())==512) {
            system(("ln -s "+thpath4ln+"/"+gd
                        +" "+mdir+"grids/"+gridsub+"/"+gd).c_str());
        }
    }
    
    //Init
    string infile;
    vector<double> binF,xlo,xhi,xav,Q2lo,Q2hi,Q2av,sigma,N,stat,uncor;
    vector<double> binFp,binFm,xlom,xhim,xavm,Q2lom,Q2him,Q2avm,sigmam,Nm,statm;

    //Include unc.s already in prel tables if available, handy to read all input
    //info from prel table in the next stage, no need to revisit binned_events
    vector<string> systnames;
    vector< vector<double> > syst, systm;

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
    int isyst = (useStat ? 10 : 9);
    if (useStat) stat = BE[9];
    for (int i=isyst; i<BE.size(); ++i) syst.push_back(BE[i]);
    for (int j=0; j!=syst.size(); ++j) {
        for (int i=0; i!=syst[j].size(); ++i) {
            syst[j][i] *= (N[i]>0 ? 1./N[i] : 0.);
        }
    }

    //Compute unc as quadratic sum of all systs here
    for (int row=0; useSyst && row!=syst[0].size(); ++row) {
        double dtmp=0;
        for (int col=0; col!=syst.size(); ++col) {
            dtmp += syst[col][row]*syst[col][row];
        }
        uncor.push_back(sqrt(dtmp));
    }
    
    //Construct generic syst unc names
    stringstream sstream;
    for (int i=0; i!=syst.size(); ++i) {
         sstream.clear();  sstream.str("");
         sstream << "syst";
         if (i>0) sstream << i;
         sstream << ":C";           //Use covariance matrix
         if (i>0) sstream << ":A";  //Additive unc.s, if many given
         systnames.push_back(sstream.str());
     }

    //Turn unc.s into %. N.B. only done for prel tables, final tables read the
    //unc.s from prel tables where they're hence already in %
    for (int i=0; i!=uncor.size(); ++i) {
        uncor[i] *= 100.;
        stat[i] *= 100.;
        for (int j=0; j!=syst.size(); ++j) syst[j][i] *= 100.;
    }

    //Write prel table
    return xFtableWriter(mdir,PDF+"/prel/",expname,nuID,origin,iexp,
                         binF,xav,Q2av,sigma,stat,uncor,systnames,syst);

} // END writeDatPrel


bool writeDatFinal(string PDF, 
                   string expname, 
                   string nuID, 
                   string origin,
                   int iexp)
{
    string mdir = "datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    
    //For generating univariate Gaussian random numbers
    default_random_engine generator;
    normal_distribution<double> distribution(0.0,1.0);
    
    //Init
    string expID = expname + origin + "_" + nuID;
    string infile;
    vector<double> dvtmp;
    vector<string> svtmp, systnames;
    vector< vector<double> > systtmp, systunc;

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
    vector< vector<double> > Mdat = xFtableReader(cols, mdir,
                                                  PDF+"/prel/",
                                                  expname,nuID,origin);
    if (Mdat.size()==0) {
        cout <<"ERROR in xFtableReader at "<<PDF+"/prel/"+expID+origin<< endl;
        return false;
    }
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
    //a correction factor f = 0.5 is applied to the correlated case
    vector<double> xsv, xsvstat, xsvunc05;
    for (int ixs=0; ixs!=xsec.size(); ++ixs) {                    
        
        //Quadratic sums of various combinations of stat & syst uncs
        double syst2=0;
        for (int i=0; i!=systunc.size(); ++i) syst2 += pow(systunc[i][ixs],2);
        double stat2 = (useStat ? pow(stat[ixs],2) : 0.);
        stat2 *= 0.0001;  syst2 *= 0.0001;  //Div by 100^2 since given in %
        double errstat       = sqrt(stat2);
        double errstatsyst   = sqrt(stat2 + syst2);
        double errstatsyst05 = sqrt(stat2 + syst2*0.25);  //0.25=0.5^2
        double rndm = distribution(generator);
        
        //Vary x-sec to obtain pseuodata 
        double xs = xsec[ixs];
        xsv.push_back(      xs*(1. +   errstatsyst*rndm));
        xsvstat.push_back(  xs*(1. +       errstat*rndm));
        xsvunc05.push_back( xs*(1. + errstatsyst05*rndm));
    }
        
    /* Write tables for final runs */

    //Including syst unc
    if (!xFtableWriter(mdir,PDF+"/syst/",
                       expname,nuID,origin,iexp,
                       binF,xav,Q2av,xsv,                 //3 bins, x-sec
                       stat,dvtmp,                        //stat, uncor
                       systnames,systunc)) return false;  //systnames, systmatrix
    
    //0.5x less syst var in pseudodata variation, otherwise full syst unc
    if (!xFtableWriter(mdir,PDF+"/systVar05/",
                       expname,nuID,origin,iexp,
                       binF,xav,Q2av,xsvunc05,            //3 bins, x-sec
                       stat,dvtmp,                        //stat, uncor
                       systnames,systunc)) return false;  //systnames, systmatrix
    
    //Stat unc only, no syst; neither as unc or in pseudodata variation
    if (!xFtableWriter(mdir,PDF+"/statOnly/",
                       expname,nuID,origin,iexp, 
                       binF,xav,Q2av,xsvstat,         //3 bins, x-sec
                       stat,dvtmp,                    //stat, uncor
                       svtmp,systtmp)) return false;  //systnames, systmatrix

    //Syst unc considered fully uncorrelated
    if (!xFtableWriter(mdir,PDF+"/uncor/",
                       expname,nuID,origin,iexp, 
                       binF,xav,Q2av,xsv,             //3 bins, x-sec
                       stat,uncor,                    //stat, uncor
                       svtmp,systtmp)) return false;  //systnames, systmatrix

    ////No variation in pseudodata, for checking chi2=0 case
    //if (!xFtableWriter(mdir,PDF+"/noVar/",
    //                   expname,nuID,origin,iexp, 
    //                   binF,xav,Q2av,xsec,            //3 bins, x-sec
    //                   stat,uncor,                    //stat, uncor
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix
    
    ////Stat unc only, as uncor
    //if (!xFtableWriter(mdir,PDF+"/statAsUncor/",
    //                   expname,nuID,origin,iexp, 
    //                   binF,xav,Q2av,xsvstat,         //3 bins, x-sec
    //                   dvtmp,stat,                    //stat, uncor; flip intended
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix

    ////Stat unc as uncor, syst[0] as stat
    //if (!xFtableWriter(mdir,PDF+"/systAsStat/",
    //                   expname,nuID,origin,iexp, 
    //                   binF,xav,Q2av,xsvstat,         //3 bins, x-sec
    //                   systunc[0],stat,               //stat, uncor; flip intended
    //                   svtmp,systtmp)) return false;  //systnames, systmatrix

    return true;
} // END writeDatFinal


bool writeCov(string PDF, 
              string expname, 
              string nuID, 
              string origin)
{
    string datadir = "./datafiles/lhc/fpf/neutrinoDIS/pseudodata/";
    ifstream in;
    fstream out;
    string line;
    stringstream sstream;
    pair<double,double> bintmp;
    int i1,i2;  //Index tmps
    double xlo,xhi,xav,Q2lo,Q2hi,Q2av,cov;

    //Read and write covariance matrix tables
    int ind1,ind2;
    
    string expID  = expname + origin + "_" + nuID;
    string exptag = replace(expID, origin, replace(origin,"_"," "));
    exptag        = replace(exptag,"Rv", "R$\\nu$");
    exptag        = replace(exptag,"_nochargediscrimination",
                                   " $\\nu_{\\mu}+\\bar{\\nu}_{\\mu}$");
    exptag        = replace(exptag,"_nub", " $\\bar{\\nu}_{\\mu}$");
    exptag        = replace(exptag,"_nu",  " $\\nu_{\\mu}$");

    
    //Deduce cov. mat bin map fom binned events
    map< int, pair<double,double> > binmap;
    string suffix = ".txt";
    string origindir = (origin=="_charm" ? "CHARM/" : "INCLUSIVE/");
    string infile = "../results/"
                  + origindir
                  + expname
                  + "/clipped_nan/clipped_nan_binned_sysevents_"
                  + expID
                  + suffix;
    if (!inStreamOpen(in,infile)) return false;
    for (int i=0; i!=2; ++i) getline(in,line);  //Skip 2 lines of header
    i1=0;  //Init
    while (getline(in,line)) {
        sstream.clear();  sstream.str("");
        sstream << setprecision(15) << line;
        sstream >> xlo >> xhi >> xav >> Q2lo >> Q2hi >> Q2av;
        bintmp.first  = xav;
        bintmp.second = Q2av;
        binmap[i1] = bintmp;            
        ++i1;
    }
    in.close();

    //Read covariance matrix entries
    string covmatfile = "../results/"
                      + origindir
                      + expname
                      + "/clipped_nan/clipped_nan_covariance_"
                      + expID
                      + suffix;
    if (!inStreamOpen(in,covmatfile)) return false;
    vector<int> i1s,i2s;
    vector<double> covs;
    while (getline(in,line)) {
        sstream.clear();  sstream.str("");
        sstream << line;
        sstream >> i1 >> i2 >> cov;
        i1s.push_back(i1);
        i2s.push_back(i2);
        covs.push_back(cov);
    }
    in.close();        

    //Write covariance matrix table in xFitter format            
    string outname = datadir + PDF + "/" + expID + ".cov";
    outStreamOpen(out,outname);
    out << "! Covariance matrix" << endl;
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
    out << "  NCorr = " << covs.size()  << "\n" << endl;
    out << "  MatrixType = 'Systematic covariance matrix'" << endl;            
    out << "&End" << endl;
    //Matrix values
    for (int i=0; i!=covs.size(); ++i) {
        out << setfill(' ') << setw(20) << setprecision(15) << binmap[i1s[i]].first;
        out << setfill(' ') << setw(20) << setprecision(15) << binmap[i1s[i]].second;
        out << setfill(' ') << setw(20) << setprecision(15) << binmap[i2s[i]].first;
        out << setfill(' ') << setw(20) << setprecision(15) << binmap[i2s[i]].second;
        out << setfill(' ') << setw(20) << setprecision(15) << covs[i] << endl;
    }
    out.close();
    cout << "Wrote " << outname << endl;
    

    return true;
} // END writeCov

int main() {
    
    //BEGIN user input
    vector<string> PDFs = {"PDF4LHC21","EPPS21nlo_CT18Anlo_W184"};
    vector<string> expnames = {"FASERv2"};
    vector<string> nuIDs = {"nu","nub","nochargediscrimination"};
    vector<string> origins = {"_inclusive","_charm"};

    //True:  write tables for preliminary xFitter run
    //False: write final tables to be used as pseudodata in fits.
    bool prel = false;

    //Flags for including uncertainties -- turn off if binned_events don't contain these
    bool useStat = true;  //Do binned_events tables contain stat unc?
    bool useSyst = true;  //              -||-              syst unc?
    //END user input

    //Write tables
    int iexp=137;
    for (string PDF : PDFs) {
        for (string expname : expnames) {
            for (string nuID : nuIDs) {
                for (string origin : origins) {
                    if (prel) {
                        if (!writeDatPrel(PDF,expname,nuID,origin,iexp,useStat,useSyst)) return -1;
                    } else {
                        if (!writeDatFinal(PDF,expname,nuID,origin,iexp)) return -1;
                    }
                    if (!writeCov(PDF,expname,nuID,origin)) return -1;
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
