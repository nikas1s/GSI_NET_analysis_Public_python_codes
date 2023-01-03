// This example program calculates the rate for each
// reaction in the Reaclib file at varius T in GK.
// It takes advantage of the classes and codes of Reaclib Reader by
//Scott Warren (warren@nscl.msu.edu)

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include "reaclibReader/NuclideSearch.h"
#include "reaclibReader/reaclibReader.h"
using std::cout;
using namespace std;

int main() {
	//loop over the T that you want to get the reaction rate at
	for (int t9 = 1; t9<=8; t9 += 1){
		char nucFile[] = "example_data/nuclides.xml";
		const char reaclibFile[] = "example_data/07_2021_AME_16";
		//int t9 = 3; // Temperature = 1 GK

		// Load the nuclide list
		NuclideSearch *nuclideSearch = NuclideSearch::newInstance();
		nuclideSearch->loadNuclideFile(nucFile);

		// Load a Reaclib file
		ReaclibReader reader;
		reader.LoadREACLIB(reaclibFile);

		// Loop through the rates
		for(int i = 0; i < reader.getRateSize(); i += 1) {
			ReaclibRate *rate = reader.getRate(i);

			// Loop through the versions of that rate
			for(int j = 0; j < rate->getVersionSize(); j += 1) {
				ReaclibVer *ver = rate->getVer(j);

				// Loop through the sets (fit parameters) of this version to calculate the
				// reaction rate (rr).
				double rr = 0;
				for(int k = 0; k < ver->getSetSize(); k += 1) {
					ReaclibSet *set = ver->getSet(k);

					// Calculate one term of the form:
					// exp(a0 + a1/T9 + a2/T9^(1/3) + a3*T9^(1/3) + a4*T9 + a5*T9^(5/3) + a6*ln(T9))
					double tmp = 0;
					tmp += set->getA0();
					tmp += set->getA1() / t9;
					tmp += set->getA2() / pow(t9, 1.0f/3.0f);
					tmp += set->getA3() * pow(t9, 1.0f/3.0f);
					tmp += set->getA4() * t9;
					tmp += set->getA5() * pow(t9, 5.0f/3.0f);
					tmp += set->getA6() * log(t9);
					tmp = exp(tmp);
					rr += tmp;
				}
				string S(rate->getEqn());
				string S2=S.substr(0,3);
				string S3=S.substr(4,5);
				string S4=S.substr(4,9);
				string S5=S.substr(0,6);
				string S6=S.substr(4,15);
				string he("he4");
				string p("p +");
				string n_e("-> n +");
				string n("n +");
				string plus("+");
				string arrow("->");
				string multiple_n("he4 + he4");
				if (S4.find(n_e) != string::npos){
				    // Show the reaction rate for that version of the reaction
				  	printf("%.3e %d %s\n", rr, t9,rate->getEqn() );}
				}
				}
}
}
