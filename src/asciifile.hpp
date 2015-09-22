/** @file asciifile.hpp
 * Header file for AsciiFile class.
 * Written by Shayan Hoshyari
 * As a part of DF_1d
 * @ingroup dr_module
 */

#ifndef ASCIIFILE_HPP
#define ASCIIFILE_HPP

#include <fstream>
#include "error.hpp"

/** @ingroup dr_module
 *  @brief Helps reading an ascii file.
 *
 * Usage:
 * \li Open file. X.open("address")
 * \li Go to Beginning: X.gotobeg()
 * \li Finding: X.find("str")
 * \li NextLine: X()
 * \li Reading: X(var_to_read, "name to give warning")
 * \li Close: X.close()
 */
class AsciiFile{
public:
	/** Stores file name.
	 */
	std::string fn;
	
	/** Stores the current line that has been read.
	 */
	std::stringstream ss;
	
	/** File stream that reads the file.
	 */
	std::fstream fs;
	
	/** Line number.
	 */
	int ln;
	
	/** Calls close.
	 */
	~AsciiFile();
	
	/** Reads the next line and puts it in ss.
	 * If reached EOF, false will be returned.
	 */
	bool next();
	
	/** Goes to next line with error checking */
	void enext();
	/** Same as enext */
	void operator() () ;
	
	/** Finds a specific string and reads its next line.
	 *	If the string is found, the line after it will be read.\n
	 *  If not EOF will be reached and false will be returned.\n
	 *	@param findee the string to look for.
	 */
	bool find(const std::string& findee);
	
		/** Find with error checking.
	 *	@param findee the string to look for.
	 */
	void efind(const std::string& findee);
	
	/** Opens a file.
	 * @param str file address.
	 */
	void open(const std::string& str);
	
	/** Closes a currently open file.
	 */
	void close();
	
	/** If the file is open goes to beginning of it.
	 */
	void goto_beg();
	
	/** Reads a parameter from the stream and generates error if faild.
			@param param the parameter that should be read.
			@param name name of the parameter
	 */
	template<class T>
	void operator() (T &param, const std::string& name ){
		FuncBegin();
		ss >> param;
		if (ss.fail()){
			Error::mess << fn << " at line " << ln << ". Failed to read " << name ;
			ERRSET();
		}
		FuncEnd();
	}

	/** Reads a string identifier from the line and checks it.
			@param name the string identifier
	 */
	void operator() (const std::string& name);
};

#endif /*ASCIIFILE_HPP*/
