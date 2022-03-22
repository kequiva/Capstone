/*******************************************************************************
Program to calculate cosmological distances in standard Lambda cosmology
Copyright (C) 2003-2021  Joshua Kempner

Version 2.1.0

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

Please send any bug fixes, enhancements, or useful comments by email to
josh@kempner.net.
*******************************************************************************/


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <map>
#include <cstring>
#include <limits>

#include "cosmo.h"

using namespace std;

string version = "2.1.0";

void help()
{
  cout << "Usage: cosmic [options]\n"
       << "Options:\n"
       << "   h=value      - value of H nought (default = 71)\n"
       << "   m=value      - value of Omega matter (default = 0.27)\n"
       << "   l=value      - value of Omega Lambda (default = 0.73)\n"
       << "   z=value      - single redshift for quick mode\n"
       << "   quiet=yes    - suppress copyright message\n"
       << "   prompt=no    - don't prompt for cosmological parameters\n"
       << "   html=yes     - output formatted in HTML"
       << "   batch=file   - run in batch mode using redshifts in \"file\"\n"
       << "   outfile=file - output batch mode results to \"file\"\n"
       << "   help=yes     - print this message\n"
       << "   version=yes  - print the version number of cosmic\n";
  exit(0);
}

void printVersion()
{
  cout << "cosmic version " << version << endl;
  exit(0);
}

void printCopyleft()
{
  cout << "cosmic version " << version << ", Copyright (C) 2003-2007 Joshua Kempner\n"
       << "cosmic comes with ABSOLUTELY NO WARRANTY; for details\n"
       << "see the accompanying license.  This is free software,\n"
       << "and you are welcome to redistribute it under certain\n"
       << "conditions; see the bundled license for details. Invoke\n"
       << "this program with \"-quiet\" or \"quiet=yes\" to suppress\n"
       << "this message.\n"
       << endl;
}

void splitArg(const string& text, string& key, string& value)
{
  int n = text.length();
  int start, stop;
  string separator = "=";

  start = text.find_first_not_of(separator);
  stop = text.find_first_of(separator, start);
  if ((stop < 0) || (stop > n)) stop = n;
  key = text.substr(start, stop - start);
  if (stop != n)
  {
    start = text.find_first_not_of(separator, stop+1);
    value = text.substr(start, n - start);
    // strip off enclosing quotes around value, if there are any
    if ('"' == value[0] || '\'' == value[0])
      value.erase(0, 1);
    n = value.length() - 1;
    if ('"' == value[n] || '\'' == value[n])
	    value.erase(n, 1);
  }
  else
  {
    if ("-" != key.substr(0, 1))
	    return;
    while ("-" == key.substr(0, 1))
	    key.erase(0, 1); // strip off leading "-"
    if ("no" == key.substr(0, 2))
	  {
	    key.erase(0, 2); // strip off leading "no"
	    value = "no";
	  }
    else
	    value = "yes";
  }
}

void processArgs(int argc, char** argv, map<string, bool>& bflags,
		 map<string, string>& sflags, map<string, double>& fflags)
{
  int badArgs = 0;
  for (int i = 1; i < argc; ++i)
  {
    // split the key=value pair
    string key, value;
    splitArg(string(argv[i]), key, value );
    if (!value.length())
	  {
	    cerr << "incomplete argument: " << argv[i] << "\n";
	    ++badArgs;
	  }
    else if (bflags.find(key) != bflags.end())
	  {
	    if ('y' == value[0] || 'Y' == value[0])
	      bflags[key] = true;
	    else if ('n' == value[0] || 'N' == value[0])
	      bflags[key] = false;
      else
      {
        cerr << "invalid value for argument '" << key << "'\n";
        ++badArgs;
      }
	  }
    else if (sflags.find(key) != sflags.end())
	    sflags[key] = value;
    else if (fflags.find(key) != fflags.end())
	  {
	    if (!isNumeric(value))
	    {
	      cerr << "invalid value for argument '" << key << "'\n";
	      ++badArgs;
	    }
	    else
	      fflags[key] = atof(value.c_str());
	  }
    else
	  {
	    cerr << "unknown argument: " << argv[i] << "\n";
	    ++badArgs;
	  }
  }
  if (badArgs)
    cerr << endl;
}

int main(int argc, char** argv)
{
    // default values for arguments
    map<string, bool> bflags;
    map<string, string> sflags;
    map<string, double> fflags;
    bflags["help"] = false;
    bflags["quiet"] = false;
    bflags["prompt"] = true;
    bflags["html"] = false;
    bflags["version"] = false;
    sflags["batch"] = "";
    sflags["outfile"] = "cosmic.out";
    fflags["h"] = 71;
    fflags["m"] = 0.27;
    fflags["l"] = 0.73;
    fflags["z"] = -1;
    
    // process arguments
    processArgs(argc, argv, bflags, sflags, fflags);
    
    // print help message if requested
    if (bflags["help"])
        help();
    
    // print version number if requested
    if (bflags["version"])
        printVersion();
    
    // print copyright notice unless quiet=yes
    if (!bflags["quiet"])
        printCopyleft();
    
    // instantiate a cosmology
    Cosmo* c = new Cosmo(fflags["h"], fflags["m"], fflags["l"]);
    if (bflags["prompt"])
        c->getCosmologyFromUser(); // prompt the user for the cosmological parameters
    
    double z = 0;
    if (fflags["z"] != -1)
    {
        c->setRedshift(fflags["z"]);
        if (bflags["html"])
            c->printAsHtml();
        else
            c->printLong();
    }
    else if (!sflags["batch"].length())
    {
        string temp;
        while ((cout << "redshift (ctrl-D to quit): ") && (cin >> temp))
        {
            if (!isNumeric(temp.c_str()))
            {
                cerr << "Redshift must be numeric\n";
                continue;
            }                
            z = atof(temp.c_str());
            if (z >= 0)
            {
                c->setRedshift(z);
                cout << "\n"; // extra blank line for readability
                if (bflags["html"])
                    c->printAsHtml();
                else
                    c->printLong();
                cout << endl;
            }
            else
            {
                cerr << "  The redshift must be a number > 0." << endl;
                cin.clear();
                cin.ignore(numeric_limits<int>::max(), '\n');
            }
        }
        cout << endl;
    }
    else
    {
        ifstream inFile(sflags["batch"].c_str());
        if (!inFile)
        {
            cerr << "Error opening batch file: " << sflags["batch"] << endl;
            return 1;
        }
        
        ofstream outFile(sflags["outfile"].c_str());
        if (!outFile)
        {
            cerr << "Error opening output file: " << sflags["outfile"] << endl;
            return 1;
        }
        
        // short message to the user
        cout << "Running in batch mode. Output will be in " << sflags["outfile"]
            << endl;
        
        // loop throught the batch file and output the redshifts to cosmic.out
        int line = 0;
        char p[100];
        c->printShortHeader(outFile); // print a couple header lines
        while (inFile.getline(p, 100))
        {
            ++line;
            z = atof(p);
            if (!z)
            {
                cerr << "Non-numeric redshift found in batch file on line" << line
                << "\nExiting with no further output" << endl;
                return 1;
            }
            
            c->setRedshift(z);
            c->printShort(outFile);
        }
        inFile.close();
        outFile.close();
    }
    
    delete c;
    return 0;
}
