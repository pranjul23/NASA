#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <time.h>
#include <vector>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

using namespace std;

int main(){

	ifstream file;
	ofstream result;
	result.open("Testing.txt",ios::trunc);

	vector<string> commands;
	commands.push_back("v"); //0
	commands.push_back("9"); //1
	commands.push_back("3"); //2
	commands.push_back("8"); //3
	commands.push_back("2"); //4
	commands.push_back("4"); //5
	commands.push_back("6"); //6
	commands.push_back("7"); //7
	commands.push_back("1"); //8
	commands.push_back("0"); //9
	commands.push_back("Enter"); //10
	commands.push_back("<Space>"); //11
	commands.push_back("["); //12
	commands.push_back("]"); //13

	string buffer, timestr, commandstr;
	double timeprev , timecurr;
	int hr, min, sec, msec;
	bool found;
	size_t commandID;

	//==== get specific files out of the directory ================
	const string target_path( "." );
	const boost::regex my_filter( "t.+txt" );
	vector< string > all_matching_files;
	boost::filesystem::directory_iterator end_itr;

	for( boost::filesystem::directory_iterator i( target_path ); i != end_itr; ++i )
	{
		// Skip if not a file
		if( !boost::filesystem::is_regular_file( i->status() ) ) continue;

		boost::smatch what;

		// Skip if no match
		if( !boost::regex_match(  i->path().filename().string(), what, my_filter ) ) continue;

		// File matches, store it
		all_matching_files.push_back(  i->path().filename().string() );
	}
	//==============================================================


	cout << "Processing " << all_matching_files.size() << " files\n";

	//==== Process the files and extract the commands ==============
	for(size_t f=0; f<all_matching_files.size(); f++){

		file.open(all_matching_files[f].c_str());
		timeprev = -1;

		cout << "file: " << all_matching_files[f].c_str() << "\n";

		while(getline(file, buffer)){

			size_t pos = buffer.find(">");

			if (pos == string::npos)
				continue;

			//======== find command =========
			commandstr = buffer.substr(pos+2, buffer.length()-(pos+1));

			found = false;
			for(size_t i=0; i<commands.size(); i++){
				if(commandstr.compare(commands[i]) == 0){
					found = true;
					commandID = i;
					break;
				}
			}

			if(!found) continue;

			//======== find timestamp ========
			timestr = buffer.substr(11, pos-12);

			istringstream (timestr.substr(0,2)) >> hr;
			istringstream (timestr.substr(3,2)) >> min;
			istringstream (timestr.substr(6,2)) >> sec;
			istringstream (timestr.substr(9,timestr.length()-9)) >> msec;

			timecurr = hr*3600.0 + min*60.0 + sec + msec/1000.0;

			//skip this command if the time distance between two commands is too small
			if(timecurr-timeprev < 0.3) continue;

			result << commandID << " ";

			timeprev = timecurr;
		}

		result << "\n";
		file.close();
	}
	//=============================================================
}







