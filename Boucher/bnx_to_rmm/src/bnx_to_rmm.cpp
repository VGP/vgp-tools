// reading a text file
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
using namespace std;

vector <float> reff;


int main (int argc, char* argv[]) {
 // read_ref();

 if(argc<3){
    cout<<"Usage: bnxcnv <bnx_filename> <rmm_filename> <enzyme_restriction_pattern>"<<endl;
    return 1;
 }

  int min_size=0;
  string line;
  ifstream bnxfile (argv[1]);


  string fname(argv[2]);

  string enz(argv[3]);

//  ofstream valfile ((fname+".val").c_str());
//  ofstream cometfile ((fname+".com").c_str());
  ofstream rmmfile ((fname).c_str());


  vector < vector <int> > rmap_values;

  vector < vector <int> > pos_values;

  vector < vector <float> > intensity_values;
  vector < vector <float> > snr_values;


  int count=0;
  long tot_frag_size=0;
  long tot_frags=0;

  long frag_size=0;
  int number_frags=0;

  if (bnxfile.is_open())
  {
    while ( getline (bnxfile,line) )
    {
      istringstream iss(line);
      string word;
      iss >> word;
      if(word.compare("1")!=0){
		continue;
      }
      else{
            vector <int> rmap;

            vector <int> pos;
            int old=0;
            while(iss >> word){
                int temp;
                istringstream (word) >> temp;
                pos.push_back(temp);
                if(old==0)
                    {old=temp; continue;}
                else{
                    rmap.push_back(temp-old);
                }

                tot_frag_size+=temp-old;
                tot_frags++;
                old=temp;

            }

            if(rmap.size()<min_size)
                continue;


            string snr_line;
            getline(bnxfile,snr_line);

            string temst;

            istringstream iss_snr(snr_line);

            float temfl;

            iss_snr>>temst;

            vector<float> snr;

            while(iss_snr>>temfl){
                snr.push_back(temfl);
            }



            string int_line;
            getline(bnxfile,int_line);

            istringstream iss_int(int_line);

            iss_int>>temst;

            vector<float> intensity;

            while(iss_int>>temfl){
                intensity.push_back(temfl);
            }

            rmap_values.push_back(rmap);
            pos_values.push_back(pos);
            intensity_values.push_back(intensity);
            snr_values.push_back(snr);
            count++;

      }
    }
    bnxfile.close();
  }

  else cout << "Unable to open file";


  rmmfile<<"r\t"<<rmap_values.size()<<"\t 1 \t"<<enz.size()<<"\t"<<enz<<endl;


  for(int i=0;i<pos_values.size();i++){

        rmmfile<<"R \t"<<pos_values[i].back()<<"\t"<<pos_values[i].size()<<"\t";

        for(int j=0;j<pos_values[i].size();j++){
            rmmfile<<pos_values[i][j]<<"\t";
        }

        rmmfile<<endl;

        rmmfile<<"I \t"<<intensity_values[i].size()<<"\t";

        for(int j=0;j<intensity_values[i].size();j++){
            rmmfile<<intensity_values[i][j]<<"\t";
        }
        rmmfile<<endl;

        rmmfile<<"N \t"<<snr_values[i].size()<<"\t";

        for(int j=0;j<snr_values[i].size();j++){
            rmmfile<<snr_values[i][j]<<"\t";
        }
        rmmfile<<endl;


//        cometfile<<"Rmap_"<<i;
//        valfile<<"Rmap_"<<i<<"\n";
//        valfile<<enz<<"\t"<<enz<<"\t";
//
//        for(int j=0;j<rmap_values[i].size();j++){
//                valfile<<" "<< (float)rmap_values[i][j]/1000;
//                cometfile<<" "<< (float)(float)rmap_values[i][j]/1000;
//
//        }
//
//        valfile<<"\n \n";
//        cometfile<<"\n";

  }






  cout<<"Total number of Rmaps: "<<count<<endl;
//  cout<<"Average fragment size: "<<(float)frag_size/number_frags;;
  return 0;
}
