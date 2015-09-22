/** @file asciifile.cpp
 * cpp file for asciifile.hpp.
 * Written by Shayan Hoshyari
 * As a part of DF_1d
 */

#include "asciifile.hpp"


void AsciiFile::open(const std::string &str){
	FuncBegin();

	fn = str;
	fs.open(fn.c_str(), std::fstream::in);
	if (!fs){
		Error::mess << str << " could not be openned.";
		ERRSET();
	}
	ln = 0;
	
	FuncEnd();
}

void AsciiFile::close(){
	FuncBegin();
	if(fs.is_open()) fs.close();
	FuncEnd();
}

AsciiFile::~AsciiFile(){
	FuncBegin();
	close();
	FuncEnd();
}

bool AsciiFile::next(){
	FuncBegin();
	char c_t = ' ';
	std::string str_t;
	bool flag;
	
	do{
		//check eof
		if (fs.eof())  return false;
		//read a line
		std::getline(fs, str_t);
		ss.clear();
		ss.str(str_t);
		ln++;
		//read a character
		ss >> c_t;
		//chek if it was a comment
		if ( (c_t == '#') || (!ss) )	flag = true;
		else {
			flag = false;
			ss.unget();
		}
  }while(flag);	
	return true;
	
	FuncEnd();
}

bool AsciiFile::find(const std::string &findee){
	FuncBegin();
	//bool flag=true;
	std::string str_t;
	
	do{
		if (!next()) return false;
		ss >> str_t;
 		if (str_t.compare(findee) == 0){
			//flag = false;
			//next();
			break;
		}
	}while(true/*flag*/);	
	return true;
	
	FuncEnd();
}

void AsciiFile::goto_beg(){
	FuncBegin();

	if(!fs.is_open()){
		Error::mess << "Can not go to beginning of not openned file";
		ERRSET();
	}
	fs.seekg(0, std::fstream::beg);
	fs.clear();
	ln = 0;

	FuncEnd();
}

void AsciiFile::enext(){
	FuncBegin();

	if(!next()) {
		Error::mess << fn << " ended unexpectedly.";
		ERRSET();
	}

	FuncEnd();
}
void AsciiFile::operator() () {
	FuncBegin();
	enext();
	FuncEnd();
}

void AsciiFile::efind(const std::string& findee){
	FuncBegin();

	if(!find(findee)) {
		Error::mess << findee << " could not be found in " << fn ;
		ERRSET();
	}

	FuncEnd();
}
	
void AsciiFile::operator() (const std::string &name){
	FuncBegin();

	std::string a_name;
	operator()(a_name, name + " String Identifier");
	if (a_name.compare(name) != 0){
		Error::mess << fn << " at line " << ln << ". Expected to read \"" 
					<< name << "\" but found \"" << a_name << "\" instead";
		ERRSET();
	}

	FuncEnd();
}
	
